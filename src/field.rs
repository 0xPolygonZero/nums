//! Utilities for finite field math.

use crate::util::{bits_usize, is_quadratic_residue};
use num_bigint::BigUint;
use num_traits::One;

/// If `n` is a quadratic residue mod `p^k`, return one of its roots.
#[must_use]
pub(crate) fn gf_usize_sqrt(n: usize, p: usize, k: usize) -> Option<usize> {
    if k == 1 {
        return fp_usize_sqrt(n, p);
    }

    // TODO: Hensel lifting for the recursive case.
    todo!()
}

/// If `n` is a quadratic residue mod `p`, return one of its roots.
#[must_use]
pub(crate) fn fp_usize_sqrt(n: usize, p: usize) -> Option<usize> {
    let n_biguint = BigUint::from(n);
    let p_biguint = BigUint::from(p);
    if n == 0 {
        Some(n)
    } else if p == 2 {
        // Tonelli–Shanks doesn't work for p = 2.
        // In GF(2), 0^2 = 0 and 1^2 = 1.
        assert!(n < 2);
        return Some(n);
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
        let mut c = fp_usize_modmul(z, q, p);
        let mut t = fp_usize_modmul(n, q, p);
        let mut r = fp_usize_modmul(n, (q + 1) / 2, p);

        while !t.is_one() {
            let i = fp_usize_squares_until_one(t, p);

            let mut b = c;
            for _ in 0..m - i - 1 {
                b *= b;
                b %= p;
            }

            m = i;
            c = b * b;
            c %= p;
            t *= c;
            t %= p;
            r *= b;
            r %= p;
        }
        Some(r)
    } else {
        None
    }
}

/// Computes `x - y` in a prime field `F_p`.
#[must_use]
#[inline]
pub(crate) fn fp_usize_sub(x: usize, y: usize, p: usize) -> usize {
    debug_assert!(x < p);
    debug_assert!(y < p);
    let (diff, underflow) = x.overflowing_sub(y);
    if underflow {
        diff.wrapping_add(p)
    } else {
        diff
    }
}

/// Computes `x^y` if `F_p`.
#[must_use]
pub(crate) fn fp_usize_modmul(x: usize, y: usize, p: usize) -> usize {
    let mut current = x;
    let mut product = 1;
    for j in 0..bits_usize(y) {
        if (y >> j & 1) != 0 {
            product *= current;
            product %= p;
        }
        current *= current;
        current %= p;
    }
    product
}

/// Given an element `x` of `F_p`, return the number of squares needed to reach 1.
#[must_use]
#[inline]
pub(crate) fn fp_usize_squares_until_one(mut x: usize, p: usize) -> usize {
    assert_ne!(x, 0);
    let mut count = 0;
    while !x.is_one() {
        x = x * x % p;
        count += 1;
    }
    count
}

#[cfg(test)]
mod tests {
    use crate::field::fp_usize_sqrt;
    use crate::BitVec;

    #[test]
    fn test_sqrt() {
        const P: usize = 31;

        let mut squares = BitVec::new(P);
        for x in 0..P {
            squares.set(x * x % P, true);
        }

        for x in 0..P {
            let opt_root = fp_usize_sqrt(x, P);
            if let Some(root) = opt_root {
                assert_eq!(root * root % P, x);
            } else {
                assert!(!squares.get(x));
            }
        }
    }
}
