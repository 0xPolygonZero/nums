use alloc::string::ToString;
use alloc::vec;
use alloc::vec::Vec;
use core::str::FromStr;

use num_bigint::BigUint;
use num_integer::{Integer, Roots};
use num_traits::{One, ToPrimitive, Zero};
use tracing::{event, instrument, Level};

use crate::bitvec::BitVec;
use crate::field::{all_sqrts_mod_prime_power, zn_sqrt, zn_sub};
use crate::gaussian_elimination::gaussian_elimination;
use crate::nullspace::nullspace_member;
use crate::util::{biguint_to_f64, distance, is_quadratic_residue, is_square, transpose};
use crate::{CompositeSplitter, FactoredInteger, SieveOfEratosthenes};

/// The number of extra smooth numbers to find before stopping the smooth number search and
/// proceeding to the nullspace search. This extra margin ensures a large nullspace, making it very
/// unlikely that all solutions are trivial ones.
const EXTRA_SMOOTH_YS: usize = 20;

/// The quadratic sieve factorization algorithm.
#[derive(Copy, Clone, Debug)]
pub struct QuadraticSieve;

impl CompositeSplitter for QuadraticSieve {
    fn divisor(&self, n: &BigUint) -> BigUint {
        event!(Level::INFO, "Splitting {}", n);
        let f = |x: &BigUint| x * x - n;

        let n_float = biguint_to_f64(n);
        let ln_n = n_float.ln();

        // Optimal bound on primes in our base is around exp(sqrt(log(n) * log(log(n))) / 2)?
        // In practice, we multiply by 2 based on experiments.
        let max_base_prime = 2.0 * ((ln_n * ln_n.ln()).sqrt() / 2.0).exp();
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

        let mut sieve = Sieve::new(n, base.clone());
        let mut sieve_iteration = 0usize;
        loop {
            let primes_with_some_odd_exp = sieve.base_counts.iter().filter(|&&c| c > 0).count();
            let min_smooth_ys = primes_with_some_odd_exp + EXTRA_SMOOTH_YS;
            if sieve.xs_with_smooth_ys.len() >= min_smooth_ys {
                break;
            }
            event!(
                Level::INFO,
                "Sieving {}, found {} of {} (may grow to {})",
                sieve_iteration,
                sieve.xs_with_smooth_ys.len(),
                min_smooth_ys,
                base.len() + EXTRA_SMOOTH_YS
            );
            // TODO: To improve it for tiny inputs, maybe do some multiple of base.len()?
            sieve.expand_by(1 << 21);
            sieve_iteration += 1;
        }
        event!(Level::INFO, "base_counts {:?}", &sieve.base_counts);

        let xs_with_smooth_ys = sieve.xs_with_smooth_ys;
        let smooth_ys: Vec<_> = xs_with_smooth_ys.iter().map(f).collect();

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
            assert_eq!(selection.len(), xs_with_smooth_ys.len());

            let mut prod_selected_xs = BigUint::one();
            let mut prod_selected_ys = BigUint::one();
            for (i, x) in xs_with_smooth_ys.iter().enumerate() {
                if selection.get(i) {
                    prod_selected_xs *= x;
                    prod_selected_ys *= f(x);
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
    /// How often each base prime occurs in smooth ys with an odd exponent.
    base_counts: Vec<usize>,

    all_odd_power_divisors: Vec<DivisorInfo>,

    /// Information about the divisors that we're actually sieving. In a perfectly precise sieve,
    /// this would just be the primes in our factor base and all powers thereof, but in our
    /// implementation it's a bit different, e.g. we don't sieve small divisors.
    sieve_divisors: Vec<SieveFactorInfo>,

    /// For each candidate, contains whatever smooth factors we've encountered so far.
    /// The `i`th entry is for `min_candidate + i`.
    y_smooth_bits: Vec<u8>,

    xs_with_smooth_ys: Vec<BigUint>,
}

#[derive(Debug)]
struct SieveFactorInfo {
    p: usize,
    p_bits: u8,
    /// The next numbers to be sieved with this factor, encoded as offsets relative to
    /// `min_candidate`. There is one entry for each square root of `n` mod `p`.
    next_offsets: [usize; 2],
}

#[derive(Debug)]
struct DivisorInfo {
    p: usize,
    power: usize,
    /// The square roots of `n`, minus `min_candidate`, mod `p^k`.
    // todo: always 2, use array
    offsets: Vec<usize>,
}

impl Sieve {
    fn new(n: &BigUint, base: Vec<usize>) -> Self {
        let min_candidate = n.sqrt() + BigUint::one();
        let min_sieve_divisor = base.last().unwrap().sqrt().sqrt().max(7);

        let mut all_odd_power_divisors = vec![];
        let mut sieve_divisors = vec![];
        assert_eq!(base[0], 2);
        for &p in &base[1..] {
            let p_bits = f64::from(p as u32).log2() as u8 + 1;
            let min_candidate_reduced = (&min_candidate % p).to_usize().unwrap();
            let n_reduced = (n % p).to_usize().unwrap();
            if let Some(sqrt) = zn_sqrt(n_reduced, p) {
                let sqrts = [sqrt, p - sqrt];

                if p >= min_sieve_divisor {
                    let next_offsets = sqrts.map(|root| zn_sub(root, min_candidate_reduced, p));
                    sieve_divisors.push(SieveFactorInfo {
                        p,
                        p_bits,
                        next_offsets,
                    });
                }

                let divisors = (1u32..)
                    .map(|exp| (exp, p.checked_pow(exp)))
                    .take_while(|&(_, power)| power.is_some())
                    .map(|(exp, power)| {
                        let power = power.unwrap();
                        let min_candidate_reduced = (&min_candidate % power).to_usize().unwrap();
                        let n_reduced = (n % power).to_usize().unwrap();
                        let power_sqrts = all_sqrts_mod_prime_power(n_reduced, p, exp as usize);
                        if exp == 1 {
                            assert_eq!(power_sqrts, sqrts.to_vec());
                        }
                        DivisorInfo {
                            p,
                            power,
                            offsets: power_sqrts
                                .into_iter()
                                .map(|root| zn_sub(root, min_candidate_reduced, power))
                                .collect(),
                        }
                    });
                all_odd_power_divisors.extend(divisors);
            }
        }

        let base_counts = vec![0; base.len()];
        Self {
            n: n.clone(),
            min_candidate,
            base,
            base_counts,
            all_odd_power_divisors,
            sieve_divisors,
            y_smooth_bits: vec![],
            xs_with_smooth_ys: vec![],
        }
    }

    fn expand_by(&mut self, amount: usize) {
        self.expand_to(self.y_smooth_bits.len() + amount);
    }

    /// Expands the sieve to the given size, i.e. the maximum offset relative to `min_candidate`.
    #[instrument(skip(self))]
    fn expand_to(&mut self, to: usize) {
        let from = self.y_smooth_bits.len();
        assert!(to >= from);

        self.y_smooth_bits.resize(to, 0);

        let min_candidate_f64 = f64::from_str(&self.min_candidate.to_string()).unwrap();

        // See the comments on target_bits below. We add 1 to account for multiplication by 2. We
        // subtract 0.5 to create a buffer for floating point errors. An incomplete y (whose factors
        // have not all been discovered) will be short by at least one bit (ignoring the looseness
        // of our target), so 0.5 is the midpoint in a sense.
        let base_target_bits = min_candidate_f64.log2() as u8 + 1;

        for factor_idx in 0..self.sieve_divisors.len() {
            let factor = &self.sieve_divisors[factor_idx];
            let mut offsets = factor.next_offsets;
            for offset in &mut offsets {
                while *offset < to {
                    // We want some rough lower bound on y that's easy to compute. Note that
                    // y = (min_candidate + offset)^2 - n
                    //   = min_candidate^2 + offset^2 + 2 min_candidate off - n
                    //   = (min_candidate^2 - n) + offset^2 + 2 min_candidate offset
                    // The dominant term here should generally be 2 min_candidate offset; we will
                    // treat this as our target. So we have
                    //      y > 2 min_candidate offset
                    // log(y) > 1 + log(min_candidate) + log(offset)
                    let target_bits = base_target_bits + (*offset as f64).log2() as u8;
                    self.y_smooth_bits[*offset] += factor.p_bits;
                    if self.y_smooth_bits[*offset] * 3 >= target_bits * 2 {
                        let x = &self.min_candidate + BigUint::from(*offset);
                        // TODO: could optimize as
                        //     y = min_candidate^2 - n + offset (2 min_candidate + offset)
                        // with min_candidate^2 - n precomputed
                        // May be unnecessary if we do this in log space.
                        let y = &x * &x - &self.n;

                        let y_smooth_part =
                            self.y_smooth_part(*offset, y.trailing_zeros().unwrap());

                        if y_smooth_part == y {
                            let exponents = exponent_vec(y.clone(), &self.base)
                                .expect("already passed fast smoothness check");
                            self.xs_with_smooth_ys.push(x);

                            for (i, exp) in exponents.iter().enumerate() {
                                self.base_counts[i] += exp % 2;
                            }

                            let factored = FactoredInteger::from_ordered_factor_counts(
                                self.base
                                    .iter()
                                    .zip(exponents)
                                    .filter(|(_p, exp)| *exp != 0)
                                    .map(|(p, exp)| (BigUint::from(*p), exp))
                                    .collect(),
                            );
                            event!(Level::DEBUG, "Found smooth y = {:?} = {:?}", &y, factored);

                            // Reset to make sure we don't revisit this y later.
                            self.y_smooth_bits[*offset] = 0;
                        } else {
                            event!(Level::DEBUG, "False positive y = {:?}", &y);
                        }
                    }
                    *offset += factor.p;
                }
            }
            self.sieve_divisors[factor_idx].next_offsets = offsets;
        }
    }

    fn y_smooth_part(&self, offset: usize, factors_of_two: u64) -> BigUint {
        // TODO: Do in log space, with f64 or as much precision as needed,
        // so there's no bignum math. Don't necessarily need 100% accuracy yet, just
        // fewer false positives & ideally no false negatives.
        let mut y_smooth_part = BigUint::one() << factors_of_two;
        for other_factor in &self.all_odd_power_divisors {
            let diff_0 = other_factor.offsets[0] as i64 - offset as i64;
            let diff_1 = other_factor.offsets[1] as i64 - offset as i64;
            let divides = diff_0.is_multiple_of(&(other_factor.power as i64))
                || diff_1.is_multiple_of(&(other_factor.power as i64));
            if divides {
                y_smooth_part *= other_factor.p;
            }
        }
        y_smooth_part
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
    use crate::qsieve::{exponent_parity_vec, Sieve};
    use crate::{CompositeSplitter, QuadraticSieve};

    #[test]
    fn test_sieve() {
        let n = BigUint::from(100u8);
        let mut sieve = Sieve::new(&n, vec![3, 17]);
        sieve.expand_to(20);
        assert_eq!(&sieve.min_candidate, &BigUint::from(11u8));
        assert_eq!(
            &sieve.y_smooth_bits,
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
