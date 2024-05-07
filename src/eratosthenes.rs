use crate::bitvec::BitVec;
use alloc::vec;
use alloc::vec::Vec;

pub struct SieveOfEratosthenes;

impl SieveOfEratosthenes {
    /// Generate all primes up to a given `limit`, exclusive.
    pub fn generate(limit: usize) -> Vec<usize> {
        if limit < 2 {
            return vec![];
        }

        let mut is_nonprime = BitVec::new(limit);
        is_nonprime.set(0, true);
        is_nonprime.set(1, true);

        for n in 2..limit {
            if is_nonprime.get(n) {
                continue;
            }
            for multiple in (n * 2..limit).step_by(n) {
                is_nonprime.set(multiple, true);
            }
        }

        (0..limit).filter(|&n| !is_nonprime.get(n)).collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::SieveOfEratosthenes;
    use alloc::vec;

    #[test]
    fn generate() {
        assert_eq!(SieveOfEratosthenes::generate(0), vec![]);
        assert_eq!(SieveOfEratosthenes::generate(1), vec![]);
        assert_eq!(SieveOfEratosthenes::generate(2), vec![]);
        assert_eq!(SieveOfEratosthenes::generate(3), vec![2]);
        assert_eq!(SieveOfEratosthenes::generate(4), vec![2, 3]);
        assert_eq!(SieveOfEratosthenes::generate(5), vec![2, 3]);
        assert_eq!(SieveOfEratosthenes::generate(6), vec![2, 3, 5]);
        assert_eq!(SieveOfEratosthenes::generate(7), vec![2, 3, 5]);
        assert_eq!(SieveOfEratosthenes::generate(8), vec![2, 3, 5, 7]);
        assert_eq!(
            SieveOfEratosthenes::generate(100),
            vec![
                2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
                83, 89, 97
            ]
        );
    }
}
