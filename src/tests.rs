use crate::{
    Factorizer, FactorizerFromSplitter, MillerRabin, PollardRho, PrimalityTest, TrialDivision,
};
use num_bigint::BigUint;

#[test]
fn primality_test_consistency() {
    for n in 0u8..100 {
        let n = BigUint::from(n);
        let res_trial_division = TrialDivision.is_prime(&n);
        let res_miller_tabin = MillerRabin { error_bits: 128 }.is_prime(&n);
        assert_eq!(
            res_miller_tabin, res_trial_division,
            "inconsistent for {}",
            n
        );
    }
}

#[test]
fn factor_test_consistency() {
    for n in 1u8..100 {
        let n = BigUint::from(n);

        let mut res_trial_division = TrialDivision.factors(&n);

        let primality_test = MillerRabin { error_bits: 128 };
        let mut res_pollard_rho = FactorizerFromSplitter {
            primality_test,
            composite_splitter: PollardRho,
        }
        .factors(&n);

        res_trial_division.sort();
        res_pollard_rho.sort();

        assert_eq!(
            res_pollard_rho, res_trial_division,
            "inconsistent for {}",
            n
        );
    }
}
