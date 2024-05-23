use std::env;
use std::str::FromStr;

use num_bigint::BigUint;
use tracing::level_filters::LevelFilter;
use tracing_forest::ForestLayer;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, Registry};

use nums::{Factorizer, FactorizerFromSplitter, MillerRabin, PollardRho, QuadraticSieve};

fn main() {
    let env_filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .from_env_lossy();

    Registry::default()
        .with(env_filter)
        .with(ForestLayer::default())
        .init();

    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        println!("Usage: {} <algorithm> <number>", args[0]);
        return;
    }
    let algorithm = args[1].as_str();
    let n = BigUint::from_str(&args[2]).expect("Failed to parse number");

    let primality_test = MillerRabin { error_bits: 128 };
    let prime_factors = match algorithm {
        "qsieve" => FactorizerFromSplitter {
            primality_test,
            composite_splitter: QuadraticSieve,
        }
        .prime_factors(&n),
        "rho" => FactorizerFromSplitter {
            primality_test,
            composite_splitter: PollardRho,
        }
        .prime_factors(&n),
        _ => panic!("Unrecognized algorithm: {}", algorithm),
    };
    println!("Factors: {:?}", prime_factors);
}
