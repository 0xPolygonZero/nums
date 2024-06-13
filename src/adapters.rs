use crate::traits::{CompositeSplitter, Factorizer, PrimalityTest};
use crate::FactoredInteger;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use tracing::{event, Level};

pub struct PrimalityTestFromFactorizer<F: Factorizer> {
    pub factorizer: F,
}

impl<F: Factorizer> PrimalityTest for PrimalityTestFromFactorizer<F> {
    fn is_prime(&self, n: &BigUint) -> bool {
        let factors = self.factorizer.prime_factors(n).to_vec();
        match factors.len() {
            0 => {
                assert!(n.is_one());
                false
            }
            1 => {
                assert_eq!(&factors[0], n);
                true
            }
            _ => false,
        }
    }
}

pub struct FactorizerFromSplitter<PT, CS>
where
    PT: PrimalityTest,
    CS: CompositeSplitter,
{
    pub primality_test: PT,
    pub composite_splitter: CS,
}

impl<PT, CS> Factorizer for FactorizerFromSplitter<PT, CS>
where
    PT: PrimalityTest,
    CS: CompositeSplitter,
{
    fn prime_factors(&self, n: &BigUint) -> FactoredInteger {
        assert!(!n.is_zero());

        if n.is_one() {
            return FactoredInteger::one();
        }

        if self.primality_test.is_prime(n) {
            return FactoredInteger::single(n.clone());
        }

        let (a, b) = self.composite_splitter.split(n);

        event!(Level::INFO, "Found split: {} = {} * {}", n, &a, &b);
        assert!(!a.is_one());
        assert!(!b.is_one());
        assert_ne!(&a, n);
        assert_ne!(&b, n);

        self.prime_factors(&a) * self.prime_factors(&b)
    }
}
