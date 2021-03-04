use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};
use rs_ec_perf::*;

/// Create a new testset for a particular RS encoding.
macro_rules! instanciate_test {
	($name:literal, $mp:ident) => {
		pub mod $mp {
			/// number of shards we want
			/// equal to number of validators
			/// unrelated to payload size
			const VALIDATOR_COUNT: usize = 2000;
			/// max payload size is 10_000_000
			/// this allows for quicker iterations with smaller
			/// payload sizes.
			const PAYLOAD_SIZE_CUTOFF: usize = 53;

			use super::super::$mp::{encode, reconstruct};
			use super::super::{roundtrip, BYTES, SMALL_RNG_SEED};
			use crate::drop_random_max;
			use criterion::{black_box, Criterion};



			#[test]
			fn criterion_encode_integrity() {
				roundtrip(encode, reconstruct, black_box(&BYTES[..PAYLOAD_SIZE_CUTOFF]), VALIDATOR_COUNT);
			}

			pub fn bench_encode(crit: &mut Criterion) {
				crit.bench_function(concat!($name, " encode"), |b| {
					b.iter(|| {
						let _ = encode(black_box(&BYTES[..PAYLOAD_SIZE_CUTOFF]), VALIDATOR_COUNT);
					})
				});
			}

			#[test]
			fn criterion_roundtrip_integrity() {
				roundtrip(encode, reconstruct, black_box(&BYTES[..PAYLOAD_SIZE_CUTOFF]), VALIDATOR_COUNT);
			}

			pub fn bench_roundtrip(crit: &mut Criterion) {
				crit.bench_function(concat!($name, " roundtrip"), |b| {
					b.iter(|| {
						roundtrip(encode, reconstruct, black_box(&BYTES[..PAYLOAD_SIZE_CUTOFF]), VALIDATOR_COUNT);
					})
				});
			}


			pub fn bench_reconstruct(crit: &mut Criterion) {
				use rand::{rngs::SmallRng, SeedableRng};

				crit.bench_function(concat!($name, " decode"), |b| {
					let encoded = encode(&BYTES[..PAYLOAD_SIZE_CUTOFF], VALIDATOR_COUNT);
					let shards = encoded.clone().into_iter().map(Some).collect::<Vec<_>>();

					let mut rng = SmallRng::from_seed(SMALL_RNG_SEED);

					b.iter(|| {
						let mut shards2 = shards.clone();
						drop_random_max(&mut shards2, VALIDATOR_COUNT, VALIDATOR_COUNT / 3, &mut rng);
						let _ = reconstruct(black_box(dbg!(shards2)), VALIDATOR_COUNT);
					})
				});
			}
		}
	};
}

pub mod tests {
	instanciate_test!("novel poly basis", novel_poly_basis);
	instanciate_test!("status quo", status_quo);
}

fn adjusted_criterion() -> Criterion {
	let crit = Criterion::default()
		.sample_size(10)
		.warm_up_time(Duration::from_secs(1))
		.measurement_time(Duration::from_secs(60));
	crit
}

criterion_group!(name = acc_novel_poly_basis; config = adjusted_criterion(); targets =  tests::novel_poly_basis::bench_roundtrip, tests::novel_poly_basis::bench_encode);
criterion_group!(name = acc_status_quo; config = adjusted_criterion(); targets =  tests::status_quo::bench_roundtrip, tests::status_quo::bench_encode);

criterion_main!(acc_novel_poly_basis, acc_status_quo);
