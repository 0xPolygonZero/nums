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
        let d = self.divisor(n);
        (n / &d, d)
    }
}

pub trait Factorizer {
    /// Returns an unsorted list of prime factors of n.
    /// Primes will appear multiple times equal to the highest power of the prime which divides n.
    fn prime_factors(&self, n: &BigUint) -> Vec<BigUint>;

    fn factor_counts(&self, n: &BigUint) -> Vec<(BigUint, usize)> {
        let mut factors = self.prime_factors(n);
        factors.sort();
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
