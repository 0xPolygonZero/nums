use crate::traits::PrimalityTest;
use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};

/// The Miller-Rabin primality test.
///
/// This is a probabilistic test, but can be configured to have a cryptographically small error
/// margin.
#[derive(Copy, Clone, Debug)]
pub struct MillerRabin {
    pub error_bits: usize,
}

impl PrimalityTest for MillerRabin {
    fn is_prime(&self, n: &BigUint) -> bool {
        if n.is_zero() || n.is_one() {
            return false;
        }

        let one = BigUint::one();
        let two = BigUint::from(2u8);
        let mut rng = rand::thread_rng();

        // Decompose n - 1 into 2^s d by splitting off factors of two.
        let n_minus_1 = n - &one;
        let s = n_minus_1.trailing_zeros().unwrap();
        let d = &n_minus_1 >> s;

        let num_bases = (self.error_bits + 1) / 2;
        'base_search: for _ in 0..num_bases {
            let base = rng.gen_biguint_range(&one, n);

            // Let x = 2^d mod n.
            let x = base.modpow(&d, n);
            if x.is_one() {
                continue; // pass
            }

            // Use iterative squaring to check if x^(2^r) == n - 1 for some r < s.
            let mut square = x;
            if square == n_minus_1 {
                continue; // pass with r = 0
            }
            for _r in 1..s {
                square = square.modpow(&two, n);
                if square == n_minus_1 {
                    continue 'base_search; // pass
                }
            }
            return false; // failed; composite detected
        }
        true
    }
}
