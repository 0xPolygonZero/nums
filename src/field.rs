//! Utilities for finite field math.

use crate::util::{bits_usize, is_quadratic_residue};
use alloc::vec;
use alloc::vec::Vec;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::One;

#[must_use]
pub(crate) fn all_sqrts_mod_prime_power(n: usize, p: usize, k: usize) -> Vec<usize> {
    if let Some(root) = sqrt_mod_prime_power(n, p, k) {
        // TODO: Need to double check this logic, are we missing some roots?
        if p == 2 {
            vec![root]
        } else {
            vec![root, p.pow(k as u32) - root]
        }
    } else {
        vec![]
    }
}

/// If `n` is a quadratic residue mod `p^k`, return one of its roots.
#[must_use]
pub(crate) fn sqrt_mod_prime_power(n: usize, p: usize, k: usize) -> Option<usize> {
    let p_pow = p.pow(k as u32);
    let n = n % p_pow;
    assert_ne!(k, 0);
    if k == 1 {
        return zn_sqrt(n % p, p);
    }
    if n == 0 {
        return Some(0);
    }

    if let Some(x) = zn_sqrt(n % p, p) {
        if x == 0 {
            assert!(n.is_multiple_of(&p));
            let p_squared = p * p;
            if n.is_multiple_of(&p_squared) {
                // x^2 = n (mod p^k)
                // (s p)^2 = t p^2 (mod p^k)
                // s^2 = t (mod p^{k-2})
                sqrt_mod_prime_power(n / p_squared, p, k - 2).map(|s| s * p)
            } else {
                None
            }
        } else {
            // See https://mathoverflow.net/a/223806/128564
            let q = p.pow(k as u32);
            let r = q / p;
            let e = (q - 2 * r + 1) / 2;
            let res = zn_pow(x, r, q) * zn_pow(n, e, q);
            Some(res % q)
        }
    } else {
        None
    }
}

/// If `n` is a quadratic residue mod `p`, return one of its roots.
#[must_use]
pub(crate) fn zn_sqrt(n: usize, p: usize) -> Option<usize> {
    assert!(n < p);
    let n_biguint = BigUint::from(n);
    let p_biguint = BigUint::from(p);
    if n == 0 {
        Some(0)
    } else if p == 2 {
        // Tonelli–Shanks doesn't work for p = 2.
        // In GF(2), 0^2 = 0 and 1^2 = 1.
        Some(n)
    } else if is_quadratic_residue(&n_biguint, &p_biguint) {
        // Apply the Tonelli–Shanks algorithm.
        // Write p - 1 = q 2^s, for some odd q.
        let s = (p - 1).trailing_zeros() as usize;
        let q = (p - 1) >> s;

        // Search for a quadratic non-residue z.
        let z = (2..p)
            .find(|&z| !is_quadratic_residue(&BigUint::from(z), &p_biguint))
            .expect("No quadratic non-residues?");

        let mut m = s;
        let mut c = zn_pow(z, q, p);
        let mut t = zn_pow(n, q, p);
        let mut r = zn_pow(n, (q + 1) / 2, p);

        while !t.is_one() {
            let i = zn_squares_until_one(t, p);

            let mut b = c;
            for _ in 0..m - i - 1 {
                b *= b;
                b %= p;
            }

            m = i;
            c = zn_mul(b, b, p);
            t = zn_mul(t, c, p);
            r = zn_mul(r, b, p);
        }
        assert_ne!(r, 0);
        Some(r)
    } else {
        None
    }
}

/// Computes `x - y` in `Z_n`.
#[must_use]
#[inline]
pub(crate) fn zn_sub(x: usize, y: usize, p: usize) -> usize {
    debug_assert!(x < p);
    debug_assert!(y < p);
    let (diff, underflow) = x.overflowing_sub(y);
    if underflow {
        diff.wrapping_add(p)
    } else {
        diff
    }
}

/// Computes `x * y` in `Z_n`.
#[must_use]
#[inline]
pub(crate) fn zn_mul(x: usize, y: usize, n: usize) -> usize {
    debug_assert!(x < n);
    debug_assert!(y < n);
    (x as u128 * y as u128 % n as u128) as usize
}

/// Computes `x^2` in `Z_n`.
#[must_use]
#[inline]
pub(crate) fn zn_square(x: usize, n: usize) -> usize {
    debug_assert!(x < n);
    zn_mul(x, x, n)
}

/// Computes `x^y` in `Z_n`.
#[must_use]
pub(crate) fn zn_pow(x: usize, y: usize, n: usize) -> usize {
    let mut current = x;
    let mut product = 1;
    for j in 0..bits_usize(y) {
        if (y >> j & 1) != 0 {
            product *= current;
            product %= n;
        }
        current *= current;
        current %= n;
    }
    product
}

/// Given an element `x` of `Z_n`, return the number of squares needed to reach 1.
#[must_use]
#[inline]
pub(crate) fn zn_squares_until_one(mut x: usize, n: usize) -> usize {
    assert_ne!(x, 0);
    assert!(x < n);
    let mut count = 0;
    while !x.is_one() {
        x = zn_square(x, n);
        count += 1;
    }
    count
}

#[cfg(test)]
mod tests {
    use crate::field::{all_sqrts_mod_prime_power, sqrt_mod_prime_power, zn_sqrt, zn_square};
    use crate::BitVec;
    use alloc::vec;

    #[test]
    fn test_sqrt_mod_prime() {
        const P: usize = 31;

        let mut squares = BitVec::new(P);
        for x in 0..P {
            squares.set(x * x % P, true);
        }

        for x in 0..P {
            let opt_root = zn_sqrt(x, P);
            if let Some(root) = opt_root {
                assert_eq!(root * root % P, x);
            } else {
                assert!(!squares.get(x));
            }
        }
    }

    #[test]
    fn test_sqrt_mod_prime_power() {
        const P: usize = 5;

        for power in 1..=5 {
            let modulus = P.pow(power as u32);
            let mut squares = BitVec::new(modulus);
            for x in 0..modulus {
                squares.set(x * x % modulus, true);
            }

            for n in 0..modulus {
                let opt_root = sqrt_mod_prime_power(n, P, power);
                if let Some(root) = opt_root {
                    assert_eq!(zn_square(root, modulus), n);
                } else {
                    assert!(!squares.get(n));
                }
            }
        }
    }

    #[test]
    fn test_all_sqrts_mod_prime_power() {
        let mut sqrts = all_sqrts_mod_prime_power(100, 3, 3);
        sqrts.sort();
        assert_eq!(sqrts, vec![10, 17]);
    }
}
