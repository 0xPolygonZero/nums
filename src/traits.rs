use alloc::vec;
use alloc::vec::Vec;
use num_bigint::BigUint;

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

pub trait Factorizer {
    /// Returns an unsorted list of prime factors of n.
    /// Primes will appear multiple times equal to the highest power of the prime which divides n.
    fn prime_factors(&self, n: &BigUint) -> Vec<BigUint>;

    fn prime_factor_counts(&self, n: &BigUint) -> Vec<(BigUint, usize)> {
        let factors = self.prime_factors(n);
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
        factor_counts
    }
}
