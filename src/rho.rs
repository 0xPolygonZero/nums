use crate::util::distance;
use crate::CompositeSplitter;
use num_bigint::BigUint;
use num_integer::{gcd, Integer};
use num_traits::One;

#[derive(Copy, Clone, Debug)]
pub struct PollardRho;

impl CompositeSplitter for PollardRho {
    fn divisor(&self, n: &BigUint) -> BigUint {
        if n.is_even() {
            return BigUint::from(2u8);
        }

        let mut addend = BigUint::from(1u8);
        loop {
            if let Some(d) = rho_attempt(n, &addend) {
                return d;
            }
            addend.inc();
            if &addend == n {
                panic!("oop");
            }
        }
    }
}

fn rho_attempt(n: &BigUint, addend: &BigUint) -> Option<BigUint> {
    let start = BigUint::from(2u8);
    let mut x = start.clone();
    let mut y = start.clone();

    let g = |x: &BigUint| (x * x + addend) % n;

    loop {
        x = g(&x);
        y = g(&g(&y));
        let d = gcd(distance(&x, &y), n.clone());

        if &d == n {
            return None; // failure
        }
        if !d.is_one() {
            return Some(d);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{FactoredInteger, Factorizer, FactorizerFromSplitter, MillerRabin, PollardRho};
    use alloc::vec;
    use core::str::FromStr;
    use num_bigint::BigUint;

    #[test]
    fn baby_bear_ext5() {
        // This n is the order of the multiplicative group of BabyBear^5.
        let n = BigUint::from_str("33075446146858977625031769923874103810955673600").unwrap();
        let primality_test = MillerRabin { error_bits: 128 };
        let factor_counts = FactorizerFromSplitter {
            primality_test,
            composite_splitter: PollardRho,
        }
        .prime_factors(&n);

        assert_eq!(
            factor_counts,
            FactoredInteger::from_ordered_factor_counts(vec![
                (BigUint::from(2u8), 27),
                (BigUint::from(3u8), 1),
                (BigUint::from(5u8), 2),
                (BigUint::from(26321u16), 1),
                (BigUint::from(1081891u32), 1),
                (BigUint::from(115384818561587951104978331u128), 1),
            ])
        );
    }
}
