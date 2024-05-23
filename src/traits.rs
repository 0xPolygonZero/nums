use alloc::vec;
use alloc::vec::Vec;
use core::fmt::{Debug, Display, Formatter};
use core::iter;
use core::ops::Mul;
use num_bigint::BigUint;
use num_traits::One;

pub trait PrimalityTest {
    fn is_prime(&self, n: &BigUint) -> bool;
}

pub trait CompositeSplitter {
    /// Undefined behavior if `n` is prime.
    fn divisor(&self, n: &BigUint) -> BigUint;

    fn split(&self, n: &BigUint) -> (BigUint, BigUint) {
        let d1 = self.divisor(n);
        let d2 = n / &d1;
        if d1 < d2 {
            (d1, d2)
        } else {
            (d2, d1)
        }
    }
}

/// Factors numbers.
pub trait Factorizer {
    /// Returns the prime factorization of the given number `n`.
    fn prime_factors(&self, n: &BigUint) -> FactoredInteger;
}

#[derive(Eq, PartialEq, Hash)]
pub struct FactoredInteger {
    /// A list of `(prime, count)` pairs, ordered starting with the smallest prime.
    /// There should be no zero counts.
    factor_counts: Vec<(BigUint, usize)>,
}

impl FactoredInteger {
    /// The value 1, which has no prime factors.
    pub fn one() -> Self {
        Self {
            factor_counts: vec![],
        }
    }

    /// A number with a single factor, which is assumed to be prime.
    pub fn single(factor: BigUint) -> Self {
        Self {
            factor_counts: vec![(factor, 1)],
        }
    }

    pub fn from_factors(mut factors: Vec<BigUint>) -> Self {
        factors.sort();
        Self::from_ordered_factors(factors)
    }

    pub fn from_ordered_factors(factors: Vec<BigUint>) -> Self {
        let mut factors_iter = factors.into_iter().peekable();
        let mut factor_counts = vec![];
        while let Some(f) = factors_iter.next() {
            let mut count = 1;
            while factors_iter.peek() == Some(&f) {
                count += 1;
                factors_iter.next();
            }
            factor_counts.push((f, count));
        }
        Self { factor_counts }
    }

    pub fn from_ordered_factor_counts(factor_counts: Vec<(BigUint, usize)>) -> Self {
        Self { factor_counts }
    }

    pub fn to_biguint(&self) -> BigUint {
        self.to_vec()
            .into_iter()
            .fold(BigUint::one(), |acc, x| acc * x)
    }

    pub fn to_vec(&self) -> Vec<BigUint> {
        self.factor_counts
            .iter()
            .flat_map(|(factor, count)| iter::repeat(factor.clone()).take(*count))
            .collect()
    }

    pub fn num_distinct_prime_factors(&self) -> usize {
        self.factor_counts.len()
    }

    pub fn is_one(&self) -> bool {
        self.factor_counts.is_empty()
    }
}

impl Mul for FactoredInteger {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut combined_vec = self.to_vec();
        combined_vec.extend(rhs.to_vec());
        Self::from_factors(combined_vec)
    }
}

impl Display for FactoredInteger {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        if self.is_one() {
            return write!(f, "1");
        }

        let mut first = true;
        for (factor, count) in &self.factor_counts {
            if first {
                first = false;
            } else {
                write!(f, " * ")?;
            }

            if *count == 1 {
                write!(f, "{}", factor)?;
            } else {
                write!(f, "{}^{}", factor, count)?;
            }
        }

        Ok(())
    }
}

impl Debug for FactoredInteger {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self)
    }
}
