use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};
use rs_ec_perf::*;

/// Create a new testset for a particular RS encoding.
macro_rules! instanciate_upper_bound_test {
	($name:literal, $mp:ident) => {
		pub mod $mp {
			/// number of shards we want
			/// equal to number of validators
			/// unrelated to payload size
			const VALIDATOR_COUNT: usize = 2000;
			/// max payload size is 10_000_000
			/// this allows for quicker iterations with smaller
			/// payload sizes.
			const PAYLOAD_SIZE_CUTOFF: usize = 10_000_000;

			use super::super::$mp::{encode, reconstruct};
			use super::super::{BYTES, SMALL_RNG_SEED};
			use crate::drop_random_max;
			use crate::WrappedShard;
			use criterion::{black_box, Criterion};
			use rand::{rngs::SmallRng, SeedableRng};

			#[test]
			fn criterion_roundtrip_integrity() {
				roundtrip(encode, reconstruct, black_box(&BYTES[..PAYLOAD_SIZE_CUTOFF]), VALIDATOR_COUNT);
			}

			pub fn bench_encode(crit: &mut Criterion) {
				crit.bench_function(concat!($name, " encode upper bound"), |b| {
					b.iter(|| {
						let _ = encode(black_box(&BYTES[..PAYLOAD_SIZE_CUTOFF]), VALIDATOR_COUNT);
					})
				});
			}

			pub fn bench_reconstruct(crit: &mut Criterion) {
				crit.bench_function(concat!($name, " decode upper bound"), |b| {
					let encoded = encode(&BYTES[..PAYLOAD_SIZE_CUTOFF], VALIDATOR_COUNT).unwrap();
					let shards = encoded.clone().into_iter().map(Some).collect::<Vec<Option<WrappedShard>>>();

					let mut rng = SmallRng::from_seed(SMALL_RNG_SEED);

					b.iter(|| {
						let mut shards2: Vec<Option<WrappedShard>> = shards.clone();
						drop_random_max(&mut shards2[..], VALIDATOR_COUNT, VALIDATOR_COUNT / 3, &mut rng);
						let _ = reconstruct(black_box(shards2), VALIDATOR_COUNT);
					})
				});
			}
		}
	};
}

pub mod tests {
	instanciate_upper_bound_test!("novel poly basis", novel_poly_basis);

	#[cfg(features = "status-quo")]
	instanciate_upper_bound_test!("status quo", status_quo);
}

macro_rules! instanciate_parameterized_test {
	($name:literal, $mp:ident) => {
		pub mod $mp {
			use super::super::$mp::{encode, reconstruct};
			use super::super::{BYTES, SMALL_RNG_SEED};
			use crate::drop_random_max;
			use crate::WrappedShard;
			use criterion::{black_box, BenchmarkId, Criterion, Throughput};
			use rand::{rngs::SmallRng, SeedableRng};

			pub const fn log2(mut x: usize) -> usize {
				let mut o: usize = 0;
				while x > 1 {
					x >>= 1;
					o += 1;
				}
				o
			}

			pub fn bench_encode(crit: &mut Criterion) {
				for validator_count in (log2(4)..log2(500)).into_iter().map(|i| 1 << i) {
					let mut group =
						crit.benchmark_group(format!(concat!($name, " encode validator_count={}"), validator_count));
					for payload_size in (1_000_usize..log2(10_000_000_usize)).into_iter().map(|i| 1 << i) {
						group.throughput(Throughput::Bytes(payload_size as u64));
						group.bench_with_input(
							BenchmarkId::from_parameter(format!("payload_size={}", payload_size)),
							&payload_size,
							|b, &payload_size| {
								b.iter(|| {
									let _ = encode(black_box(&BYTES[..payload_size]), black_box(validator_count));
								})
							},
						);
					}
					group.finish();
				}
			}

			pub fn bench_reconstruct(crit: &mut Criterion) {
				let mut rng = SmallRng::from_seed(SMALL_RNG_SEED);

				for validator_count in (log2(4)..log2(500)).into_iter().map(|i| 1 << i) {
					let mut group = crit
						.benchmark_group(format!(concat!($name, "reconstruct validator_count={}"), validator_count));
					for payload_size in (1_000_usize..log2(10_000_000)).into_iter().map(|i| 1 << i) {
						group.throughput(Throughput::Bytes(payload_size as u64));
						group.bench_with_input(
							BenchmarkId::from_parameter(format!(" payload_size={}", payload_size)),
							&payload_size,
							|b, &payload_size| {
								let encoded = encode(&BYTES[..payload_size], validator_count).unwrap();
								let shards = encoded.clone().into_iter().map(Some).collect::<Vec<_>>();

								b.iter(|| {
									let mut shards2: Vec<Option<WrappedShard>> = shards.clone();
									drop_random_max(&mut shards2[..], validator_count, validator_count / 3, &mut rng);
									let _ = reconstruct(black_box(shards2), black_box(validator_count));
								})
							},
						);
					}
					group.finish();
				}
			}
		}
	};
}

pub mod parameterized {
	instanciate_parameterized_test!("novel poly basis", novel_poly_basis);

	#[cfg(features = "status-quo")]
	instanciate_parameterized_test!("status quo", status_quo);
}

fn parameterized_criterion() -> Criterion {
	let crit = Criterion::default()
		.sample_size(10)
		.warm_up_time(Duration::from_millis(200))
		.measurement_time(Duration::from_secs(10));
	crit
}

criterion_group!(
	name = plot_paramterized;
	config = parameterized_criterion();
	targets =
	parameterized::novel_poly_basis::bench_encode,
	parameterized::novel_poly_basis::bench_reconstruct,
	// parameterized::status_quo::bench_encode,
	// parameterized::status_quo::bench_reconstruct,
);

fn adjusted_criterion() -> Criterion {
	let crit = Criterion::default()
		.sample_size(10)
		.warm_up_time(Duration::from_secs(1))
		.measurement_time(Duration::from_secs(70));
	crit
}

criterion_group!(
name = upper_bounds;
config = adjusted_criterion();
targets =
	tests::novel_poly_basis::bench_encode,
	tests::novel_poly_basis::bench_reconstruct,
	// too slow, takes 30 minutes for 10 test runs
	// tests::status_quo::bench_encode,
	// tests::status_quo::bench_reconstruct,
);

criterion_main!(plot_paramterized, upper_bounds);
