use alloc::vec;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Zero};

use crate::traits::PrimalityTest;
use crate::{FactoredInteger, Factorizer};

/// The trial division method for primality testing or factorization.
pub struct TrialDivision;

impl PrimalityTest for TrialDivision {
    fn is_prime(&self, n: &BigUint) -> bool {
        if n.is_zero() || n.is_one() {
            return false;
        }

        let mut divisor = BigUint::from(2u8);
        while &divisor * &divisor <= *n {
            if n.is_multiple_of(&divisor) {
                return false;
            }
            divisor.inc();
        }
        true
    }
}

impl Factorizer for TrialDivision {
    fn prime_factors(&self, n: &BigUint) -> FactoredInteger {
        assert!(!n.is_zero());

        let mut n = n.clone();
        let mut factors = vec![];
        let mut divisor = BigUint::from(2u8);

        while &divisor * &divisor <= n {
            let (quotient, remainder) = n.div_rem(&divisor);
            if remainder.is_zero() {
                factors.push(divisor.clone());
                n = quotient;
            } else {
                divisor.inc();
            }
        }

        if !n.is_one() {
            factors.push(n);
        }

        FactoredInteger::from_ordered_factors(factors)
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec;
    use num_bigint::BigUint;
    use num_traits::{One, Zero};

    use crate::{FactoredInteger, Factorizer, PrimalityTest, TrialDivision};

    #[test]
    fn is_prime() {
        assert!(!TrialDivision.is_prime(&BigUint::zero()));
        assert!(!TrialDivision.is_prime(&BigUint::one()));
        assert!(TrialDivision.is_prime(&BigUint::from(2u8)));
        assert!(TrialDivision.is_prime(&BigUint::from(3u8)));
        assert!(!TrialDivision.is_prime(&BigUint::from(4u8)));
        assert!(TrialDivision.is_prime(&BigUint::from(5u8)));
        assert!(!TrialDivision.is_prime(&BigUint::from(6u8)));
        assert!(TrialDivision.is_prime(&BigUint::from(7u8)));
        assert!(!TrialDivision.is_prime(&BigUint::from(8u8)));
    }

    #[test]
    fn factor() {
        assert_eq!(
            TrialDivision.prime_factors(&BigUint::one()),
            FactoredInteger::one()
        );
        assert_eq!(
            TrialDivision.prime_factors(&BigUint::from(2u8)),
            FactoredInteger::single(BigUint::from(2u8))
        );
        assert_eq!(
            TrialDivision.prime_factors(&BigUint::from(3u8)),
            FactoredInteger::single(BigUint::from(3u8))
        );
        assert_eq!(
            TrialDivision.prime_factors(&BigUint::from(4u8)),
            FactoredInteger::from_factors(vec![BigUint::from(2u8), BigUint::from(2u8)])
        );
        assert_eq!(
            TrialDivision.prime_factors(&BigUint::from(5u8)),
            FactoredInteger::single(BigUint::from(5u8))
        );
        assert_eq!(
            TrialDivision.prime_factors(&BigUint::from(6u8)),
            FactoredInteger::from_factors(vec![BigUint::from(2u8), BigUint::from(3u8)])
        );
    }
}
