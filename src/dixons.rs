use alloc::vec;
use alloc::vec::Vec;

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Zero};
use tracing::{event, Level};

use crate::bitvec::BitVec;
use crate::gaussian_elimination::gaussian_elimination;
use crate::nullspace::nullspace_member;
use crate::util::{biguint_to_f64, distance, is_quadratic_residue, is_square, transpose};
use crate::{CompositeSplitter, SieveOfEratosthenes};

const MARGIN: usize = 10;

/// Dixon's random squares method for factorization.
#[derive(Copy, Clone, Debug)]
pub struct DixonsRandomSquares;

impl CompositeSplitter for DixonsRandomSquares {
    fn divisor(&self, n: &BigUint) -> BigUint {
        event!(Level::DEBUG, "Splitting {}", n);
        let q = |x: &BigUint| x * x - n;

        let n_float = biguint_to_f64(n);
        let ln_n = n_float.ln();
        // Optimal bound on primes in our base is around exp(sqrt(log(n) * log(log(n))) / 2).
        let max_base_prime = ((ln_n * ln_n.ln()).sqrt() / 2.0).exp();
        // Make sure our base isn't tiny.
        let max_base_prime = (max_base_prime as usize).max(300);
        event!(Level::DEBUG, "Using base with max of {}", max_base_prime);

        let base = SieveOfEratosthenes::generate(max_base_prime);
        let base: Vec<_> = base
            .into_iter()
            .filter(|&p| is_quadratic_residue(n, &BigUint::from(p)))
            .collect();
        event!(Level::DEBUG, "{} primes in base", base.len());

        let mut x = n.sqrt();
        let mut xs_with_sparse_y = vec![];
        let mut sparse_y_parity_vecs = vec![];

        let is_nontrivial_divisor = |d: &BigUint| !d.is_one() && d != n && n.is_multiple_of(d);

        while xs_with_sparse_y.len() < base.len() + MARGIN {
            x.inc();
            let y = q(&x);

            let x_gcd = x.gcd(n);
            let y_gcd = y.gcd(n);
            if is_nontrivial_divisor(&x_gcd) {
                event!(Level::DEBUG, "Found common factor with x={}", x);
                return x_gcd;
            }
            if is_nontrivial_divisor(&y_gcd) {
                event!(Level::DEBUG, "Found common factor with y={}", y);
                return y_gcd;
            }

            if let Some(vec) = exponent_parity_vec(y.clone(), &base) {
                event!(
                    Level::DEBUG,
                    "Found smooth number ({} of {}): {}^2 - n = {}",
                    xs_with_sparse_y.len(),
                    base.len() + MARGIN,
                    &x,
                    &y
                );
                xs_with_sparse_y.push(x.clone());
                sparse_y_parity_vecs.push(vec);
            }
        }

        // Transpose it to a matrix with n rows and n + MARGIN columns.
        let mut sparse_y_parity_vecs_t = transpose(sparse_y_parity_vecs);

        // Use Gaussian elimination to transform the matrix into row echelon form.
        gaussian_elimination(&mut sparse_y_parity_vecs_t);

        // Find a mask such that the mask applied to sparse_y_parity_vecs yields zero.
        // Thus, the mask applied to sparse y's yields something with even exponents
        // on each prime in our base.
        let mut solution_index = BigUint::one();
        loop {
            let selection = nullspace_member(&sparse_y_parity_vecs_t, &solution_index)
                .expect("No more solutions to try; need to expand margin");
            event!(Level::DEBUG, "Selection: {:?}", &selection);

            let mut prod_selected_xs = BigUint::one();
            let mut prod_selected_ys = BigUint::one();
            for (i, x) in xs_with_sparse_y.iter().enumerate() {
                if selection.get(i) {
                    prod_selected_xs *= x;
                    prod_selected_ys *= q(x);
                }
            }
            assert!(is_square(&prod_selected_ys), "RHS is not a square");
            let a = prod_selected_xs % n;
            let b = prod_selected_ys.sqrt() % n;
            assert_eq!(a.pow(2) % n, b.pow(2) % n, "Expected a^2 = b^2 mod n");

            // Now we have lhs^2 ≡ rhs^2 (mod n), so if lhs != ±rhs,
            // then gcd(lhs - rhs, n) is a nontrivial factor of n.
            if a != b && a != (n - &b) % n {
                let gcd = distance(&a, &b).gcd(n);
                assert!(!gcd.is_one());
                assert_ne!(&gcd, n);
                return gcd;
            } else {
                event!(
                    Level::DEBUG,
                    "Found solution with {} = ±{} (mod n); retrying",
                    a,
                    b
                );
                solution_index.inc();
            }
        }
    }
}

/// If `x` is smooth with respect to `base`, return the parity of its exponents.
fn exponent_vec(mut x: BigUint, base: &[usize]) -> Option<Vec<usize>> {
    let mut vec = Vec::with_capacity(base.len());
    for &b in base.iter() {
        let num_factors = division_loop(&mut x, b);
        vec.push(num_factors);
    }
    if x.is_one() {
        Some(vec)
    } else {
        None
    }
}

/// If `x` is smooth with respect to `base`, return its exponents.
fn exponent_parity_vec(x: BigUint, base: &[usize]) -> Option<BitVec> {
    exponent_vec(x, base).map(|vec| vec.into_iter().map(|exp| exp.is_odd()).collect())
}

/// Divide n by divisor repeatedly for as long as it's divisible; return the number of divisions.
// TODO: Return parity only?
fn division_loop(n: &mut BigUint, divisor: usize) -> usize {
    let d = BigUint::from(divisor);
    let mut count = 0;
    loop {
        let (quotient, remainder) = n.div_rem(&d);
        if remainder.is_zero() {
            *n = quotient;
            count += 1;
        } else {
            return count;
        }
    }
}

#[cfg(test)]
mod tests {
    use core::str::FromStr;

    use num_bigint::BigUint;
    use tracing::level_filters::LevelFilter;
    use tracing_forest::ForestLayer;
    use tracing_subscriber::layer::SubscriberExt;
    use tracing_subscriber::util::SubscriberInitExt;
    use tracing_subscriber::{EnvFilter, Registry};

    use crate::bitvec::BitVec;
    use crate::dixons::exponent_parity_vec;
    use crate::{CompositeSplitter, DixonsRandomSquares};

    #[test]
    fn test_exponent_parity_vec() {
        assert_eq!(
            exponent_parity_vec(BigUint::from(24u8), &[2, 3]),
            Some(BitVec::from([true, true]))
        );
        assert_eq!(
            exponent_parity_vec(BigUint::from(48u8), &[2, 3]),
            Some(BitVec::from([false, true]))
        );
        assert_eq!(exponent_parity_vec(BigUint::from(30u8), &[2, 3]), None);
        assert_eq!(exponent_parity_vec(BigUint::from(7u8), &[2, 3]), None);
    }

    #[test]
    fn factor_124833796738414950008916111503921() {
        let env_filter = EnvFilter::builder()
            .with_default_directive(LevelFilter::DEBUG.into())
            .from_env_lossy();

        Registry::default()
            .with(env_filter)
            .with(ForestLayer::default())
            .init();

        let p = BigUint::from_str("1081891").unwrap();
        let q = BigUint::from_str("115384818561587951104978331").unwrap();
        let n = &p * &q;
        assert_eq!(DixonsRandomSquares.split(&n), (p, q));
    }
}
