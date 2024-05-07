use crate::bitvec::BitVec;
use alloc::vec::Vec;
use num_bigint::BigUint;
use num_traits::{CheckedSub, One, Zero};

#[must_use]
pub(crate) fn is_square(n: &BigUint) -> bool {
    &n.sqrt().pow(2) == n
}

/// Determines whether `n` is a quadratic residue modulo a prime `p`.
#[must_use]
pub(crate) fn is_quadratic_residue(n: &BigUint, p: &BigUint) -> bool {
    if p == &BigUint::from(2u8) {
        return true;
    }
    let n = n % p;
    if n.is_zero() {
        return true;
    }

    // We apply Euler's criterion.
    let p_minus_one = p - BigUint::one();
    let exponent = &p_minus_one >> 1;
    let power = n.modpow(&exponent, p);
    if power.is_one() {
        true
    } else if power == p_minus_one {
        false
    } else {
        panic!("expected 1 or -1, got {}", power);
    }
}

// pub(crate) fn sqrts(_n: &BigUint) -> Option<(BigUint, BigUint)> {
//     todo!()
// }

/// Computes `ceil(a / b)`. Assumes `a + b` does not overflow.
#[must_use]
pub const fn ceil_div_usize(a: usize, b: usize) -> usize {
    (a + b - 1) / b
}

#[must_use]
pub(crate) fn distance(x: &BigUint, y: &BigUint) -> BigUint {
    x.checked_sub(y).unwrap_or_else(|| y - x)
}

#[must_use]
pub(crate) fn transpose(mat: Vec<BitVec>) -> Vec<BitVec> {
    // Naive transpose for now. If this somehow becomes a bottleneck, we should use a more
    // cache-friendly algorithm.

    let input_major_len = mat.len();
    let input_minor_len = mat[0].len();
    mat.iter()
        .for_each(|vec| assert_eq!(vec.len(), input_minor_len));

    (0..input_minor_len)
        .map(|result_major_idx| {
            let mut vec = BitVec::new(input_major_len);
            #[allow(clippy::needless_range_loop)]
            for result_minor_idx in 0..input_major_len {
                vec.set(
                    result_minor_idx,
                    mat[result_minor_idx].get(result_major_idx),
                );
            }
            vec
        })
        .collect()
}

#[must_use]
pub(crate) fn biguint_to_f64(n: &BigUint) -> f64 {
    n.to_str_radix(10).parse::<f64>().unwrap()
}

#[cfg(test)]
mod tests {
    use crate::bitvec::BitVec;
    use crate::util::transpose;
    use alloc::vec;

    #[test]
    fn test_transpose() {
        assert_eq!(
            transpose(vec![
                BitVec::from([true, false]),
                BitVec::from([true, true]),
                BitVec::from([false, false]),
            ]),
            vec![
                BitVec::from([true, true, false]),
                BitVec::from([false, true, false]),
            ]
        );
    }
}
