use crate::{
    Factorizer, FactorizerFromSplitter, MillerRabin, PollardRho, PrimalityTest, QuadraticSieve,
    TrialDivision,
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
    let primality_test = MillerRabin { error_bits: 128 };

    let rho_factorizer = FactorizerFromSplitter {
        primality_test,
        composite_splitter: PollardRho,
    };

    let qs_factorizer = FactorizerFromSplitter {
        primality_test,
        composite_splitter: QuadraticSieve,
    };

    for n in 1u8..100 {
        let n = BigUint::from(n);

        let res_trial_division = TrialDivision.prime_factors(&n);

        let res_rho = rho_factorizer.prime_factors(&n);
        let res_qs = qs_factorizer.prime_factors(&n);

        assert_eq!(res_rho, res_trial_division, "inconsistent for {}", n);
        assert_eq!(res_qs, res_trial_division, "inconsistent for {}", n);
    }
}
