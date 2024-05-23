use alloc::string::ToString;
use alloc::vec;
use alloc::vec::Vec;
use core::str::FromStr;

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use tracing::{event, instrument, Level};

use crate::bitvec::BitVec;
use crate::field::{all_sqrts_mod_prime_power, zn_sub};
use crate::gaussian_elimination::gaussian_elimination;
use crate::nullspace::nullspace_member;
use crate::util::{biguint_to_f64, distance, is_quadratic_residue, is_square, transpose};
use crate::{CompositeSplitter, SieveOfEratosthenes};

const MARGIN: usize = 20;

/// The quadratic sieve factorization algorithm.
#[derive(Copy, Clone, Debug)]
pub struct QuadraticSieve;

impl CompositeSplitter for QuadraticSieve {
    fn divisor(&self, n: &BigUint) -> BigUint {
        event!(Level::INFO, "Splitting {}", n);
        let q = |x: &BigUint| x * x - n;

        let n_float = biguint_to_f64(n);
        let ln_n = n_float.ln();
        // Optimal bound on primes in our base is around exp(sqrt(log(n) * log(log(n))) / 2).
        let max_base_prime = ((ln_n * ln_n.ln()).sqrt() / 2.0).exp();
        // Make sure our base isn't tiny.
        let max_base_prime = (max_base_prime as usize).max(300);
        event!(Level::INFO, "Using base with max of {}", max_base_prime);

        let base = SieveOfEratosthenes::generate(max_base_prime);
        let base: Vec<_> = base
            .into_iter()
            .filter(|&p| is_quadratic_residue(n, &BigUint::from(p)))
            .collect();
        event!(Level::INFO, "{} primes in base", base.len());

        // If any p in our base divides n, just return that divisor.
        // Normally we wouldn't bother with such edge cases, but the sieve code doesn't like p | n,
        // as it would mean there's only one square root mod p.
        for &p in &base {
            let p_biguint = BigUint::from(p);
            if n.is_multiple_of(&p_biguint) {
                return p_biguint;
            }
        }

        let min_smooth_ys = base.len() + MARGIN;
        let mut sieve = Sieve::new(n, base.clone(), max_base_prime);
        let mut sieve_iteration = 0usize;
        while sieve.xs_with_smooth_ys.len() < min_smooth_ys {
            event!(
                Level::INFO,
                "Sieving {}, found {} of {}",
                sieve_iteration,
                sieve.xs_with_smooth_ys.len(),
                min_smooth_ys
            );
            // TODO: To improve it for tiny inputs, maybe do some multiple of base.len()?
            sieve.expand_by(1 << 21);
            sieve_iteration += 1;
        }

        let xs_with_smooth_ys = sieve.xs_with_smooth_ys;
        let smooth_ys: Vec<_> = xs_with_smooth_ys.iter().map(q).collect();

        let sparse_y_parity_vecs: Vec<_> = smooth_ys
            .into_iter()
            .map(|y| exponent_parity_vec(y, &base).expect("Sieved y is not actually smooth?"))
            .collect();

        // TODO: Bit-reverse each vec in sparse_y_parity_vecs? Might improve GE performance.

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
            event!(Level::INFO, "Selection: {:?}", &selection);

            let mut prod_selected_xs = BigUint::one();
            let mut prod_selected_ys = BigUint::one();
            for (i, x) in xs_with_smooth_ys.iter().enumerate() {
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
                    Level::INFO,
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

struct Sieve {
    n: BigUint,
    min_candidate: BigUint,
    base: Vec<usize>,

    /// For each base element or power thereof, the next numbers to be sieved with it, encoded as
    /// offsets relative to `min_candidate`. For odd primes, there are two such offsets, each
    /// corresponding to one of two square roots mod the prime.
    factors: Vec<SieveFactorInfo>,

    /// For each candidate, contains whatever smooth factors we've encountered so far.
    /// The `i`th entry is for `min_candidate + i`.
    y_smooth_factors: Vec<f64>,

    xs_with_smooth_ys: Vec<BigUint>,
}

struct SieveFactorInfo {
    prime_bits: f64,
    prime_power: usize,
    /// The next numbers to be sieved with this factor, encoded as offsets relative to
    /// `min_candidate`. There is one entry for each square root of `n` mod `prime_power`.
    next_offsets: Vec<usize>,
}

impl Sieve {
    fn new(n: &BigUint, base: Vec<usize>, max_divisor: usize) -> Self {
        let min_candidate = n.sqrt() + BigUint::one();

        let mut base_progress = vec![];
        for &p in &base {
            let p_bits = f64::from(p as u32).log2();
            let powers = (1u32..)
                .map(|exp| (exp, p.pow(exp)))
                .take_while(|&(_, power)| power <= max_divisor);
            for (exp, power) in powers {
                let min_candidate_reduced = (&min_candidate % power).to_usize().unwrap();
                let n_reduced = (n % power).to_usize().unwrap();
                let square_roots: Vec<_> = all_sqrts_mod_prime_power(n_reduced, p, exp as usize)
                    .into_iter()
                    .map(|root| zn_sub(root, min_candidate_reduced, power))
                    .collect();
                if !square_roots.is_empty() {
                    base_progress.push(SieveFactorInfo {
                        prime_bits: p_bits,
                        prime_power: power,
                        next_offsets: square_roots,
                    });
                }
            }
        }

        Self {
            n: n.clone(),
            min_candidate,
            base,
            factors: base_progress,
            y_smooth_factors: vec![],
            xs_with_smooth_ys: vec![],
        }
    }

    fn expand_by(&mut self, amount: usize) {
        self.expand_to(self.y_smooth_factors.len() + amount);
    }

    /// Expands the sieve to the given size, i.e. the maximum offset relative to `min_candidate`.
    #[instrument(skip(self))]
    fn expand_to(&mut self, to: usize) {
        let from = self.y_smooth_factors.len();
        assert!(to >= from);

        self.y_smooth_factors.resize(to, 0f64);

        let min_candidate_f64 = f64::from_str(&self.min_candidate.to_string()).unwrap();
        let n_f64 = f64::from_str(&self.n.to_string()).unwrap();

        for ref mut factor in self.factors.iter_mut() {
            for offset in &mut factor.next_offsets {
                while *offset < to {
                    let x_f64 = min_candidate_f64 + *offset as f64;
                    let y_f64 = x_f64 * x_f64 - n_f64;
                    let target_bits = y_f64.log2();
                    let current_factor_bits = &mut self.y_smooth_factors[*offset];
                    *current_factor_bits += factor.prime_bits;
                    if *current_factor_bits + 0.5 >= target_bits {
                        let x = &self.min_candidate + BigUint::from(*offset);
                        let y = &x * &x - &self.n;
                        assert!(exponent_vec(y, &self.base).is_some(), "false positive");
                        self.xs_with_smooth_ys.push(x);
                    }
                    *offset += factor.prime_power;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec;
    use core::str::FromStr;

    use num_bigint::BigUint;
    use tracing::level_filters::LevelFilter;
    use tracing_forest::ForestLayer;
    use tracing_subscriber::layer::SubscriberExt;
    use tracing_subscriber::util::SubscriberInitExt;
    use tracing_subscriber::{EnvFilter, Registry};

    use crate::bitvec::BitVec;
    use crate::quadratic_sieve::{exponent_parity_vec, Sieve};
    use crate::{CompositeSplitter, QuadraticSieve};

    #[test]
    fn test_sieve() {
        let n = BigUint::from(100u8);
        let mut sieve = Sieve::new(&n, vec![3, 17], 30);
        sieve.expand_to(20);
        assert_eq!(&sieve.min_candidate, &BigUint::from(11u8));
        assert_eq!(
            &sieve.y_smooth_factors,
            &[
                // We start from x=11.
                3f64.log2(),
                0f64,
                3f64.log2(),
                3f64.log2(),
                0f64,
                3f64.log2(),
                3f64.log2() + 3f64.log2() + 3f64.log2(),
                0f64,
                3f64.log2() + 3f64.log2(),
                3f64.log2(),
                0f64,
                3f64.log2(),
                3f64.log2(),
                17f64.log2(),
                3f64.log2(),
                3f64.log2() + 3f64.log2(),
                17f64.log2(),
                3f64.log2() + 3f64.log2(),
                3f64.log2(),
                0f64,
            ]
        );
    }

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
            .with_default_directive(LevelFilter::INFO.into())
            .from_env_lossy();

        Registry::default()
            .with(env_filter)
            .with(ForestLayer::default())
            .init();

        let p = BigUint::from_str("1081891").unwrap();
        let q = BigUint::from_str("115384818561587951104978331").unwrap();
        let n = &p * &q;
        assert_eq!(QuadraticSieve.split(&n), (p, q));
    }
}
