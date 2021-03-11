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

			use crate::drop_random_max;
			use crate::$mp::{encode, reconstruct};
			use crate::WrappedShard;
			use crate::{BYTES, SMALL_RNG_SEED};
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
	instanciate_upper_bound_test!("novel poly", novel_poly_basis);

	#[cfg(feature = "status-quo")]
	instanciate_upper_bound_test!("status quo", status_quo);
}

pub mod parameterized {
	use std::ops::Range;

	use crate::drop_random_max;
	use crate::WrappedShard;
	use crate::{BYTES, SMALL_RNG_SEED};
	use criterion::{black_box, BenchmarkId, Criterion};
	use rand::{rngs::SmallRng, SeedableRng};

	const STEPS_VALIDATORS: usize = 3;
	const STEPS_PAYLOAD: usize = 7;

	pub fn steped(range: Range<usize>, steps: usize) -> impl Iterator<Item = usize> {
		assert!(steps > 1);
		let step = range.len() / (steps - 1);
		range.into_iter().step_by(step)
	}

	pub fn bench_encode_2d(crit: &mut Criterion) {
		for validator_count in steped(4..100, STEPS_VALIDATORS) {
			let mut group = crit.benchmark_group(format!("parameterized encode validator_count={}", validator_count));

			for payload_size in steped(1_000..10_000_000, STEPS_PAYLOAD) {
				encode_add_to_group(&mut group, payload_size, validator_count, payload_size);
			}
			group.finish();
		}
	}
	pub fn bench_encode_fixed_1mb_payload(crit: &mut Criterion) {
		let payload_size: usize = 1_000_000;

		let mut group = crit.benchmark_group("parameterized encode fixed payload");
		for validator_count in steped(4..1000, STEPS_VALIDATORS * 4) {
			encode_add_to_group(&mut group, validator_count, validator_count, payload_size);
		}
		group.finish();
	}

	pub fn bench_reconstruct_2d(crit: &mut Criterion) {
		let mut rng = SmallRng::from_seed(SMALL_RNG_SEED);

		for validator_count in steped(4..100, STEPS_VALIDATORS) {
			let mut group =
				crit.benchmark_group(format!("parameterized reconstruct validator_count={}", validator_count));
			for payload_size in steped(1_000..10_000_000, STEPS_PAYLOAD) {
				reconstruct_add_to_group(&mut group, payload_size, validator_count, payload_size, &mut rng);
			}
			group.finish();
		}
	}

	pub fn bench_reconstruct_fixed_1mb_payload(crit: &mut Criterion) {
		let mut rng = SmallRng::from_seed(SMALL_RNG_SEED);
		let payload_size: usize = 1_000_000;
		let mut group = crit.benchmark_group("parameterized reconstruct fixed payload");

		for validator_count in steped(4..1000, STEPS_VALIDATORS * 4) {
			reconstruct_add_to_group(&mut group, validator_count, validator_count, payload_size, &mut rng);
		}
		group.finish();
	}

	fn encode_add_to_group<M: criterion::measurement::Measurement>(
		group: &mut criterion::BenchmarkGroup<M>,
		param: impl ToString,
		validator_count: usize,
		payload_size: usize,
	) {
		{
			use crate::novel_poly_basis::encode;

			group.bench_with_input(
				BenchmarkId::new("novel-poly-encode", param.to_string()),
				&payload_size,
				|b, &payload_size| {
					{
						b.iter(|| {
							let _ = encode(black_box(&BYTES[..payload_size]), black_box(validator_count));
						})
					}
				},
			);
		}
		#[cfg(feature = "status-quo")]
		{
			use crate::status_quo::encode;

			group.bench_with_input(
				BenchmarkId::new("status-quo-encode", param.to_string()),
				&payload_size,
				|b, &payload_size| {
					b.iter(|| {
						let _ = encode(black_box(&BYTES[..payload_size]), black_box(validator_count));
					})
				},
			);
		}
	}

	fn reconstruct_add_to_group<M: criterion::measurement::Measurement>(
		group: &mut criterion::BenchmarkGroup<M>,
		param: impl ToString,
		validator_count: usize,
		payload_size: usize,
		rng: &mut SmallRng,
	) {
		{
			use crate::novel_poly_basis::{encode, reconstruct};

			group.bench_with_input(
				BenchmarkId::new("novel-poly-reconstruct", param.to_string()),
				&payload_size,
				|b, &payload_size| {
					let encoded = encode(&BYTES[..payload_size], validator_count).unwrap();
					let shards = encoded.clone().into_iter().map(Some).collect::<Vec<_>>();

					b.iter(|| {
						let mut shards2: Vec<Option<WrappedShard>> = shards.clone();
						drop_random_max(&mut shards2[..], validator_count, validator_count / 3, rng);
						let _ = reconstruct(black_box(shards2), black_box(validator_count));
					})
				},
			);
		}

		#[cfg(feature = "status-quo")]
		{
			use crate::status_quo::{encode, reconstruct};

			group.bench_with_input(
				BenchmarkId::new("status-quo-reconstruct", param.to_string()),
				&payload_size,
				|b, &payload_size| {
					let encoded = encode(&BYTES[..payload_size], validator_count).unwrap();
					let shards = encoded.clone().into_iter().map(Some).collect::<Vec<_>>();

					b.iter(|| {
						let mut shards2: Vec<Option<WrappedShard>> = shards.clone();
						drop_random_max(&mut shards2[..], validator_count, validator_count / 3, rng);
						let _ = reconstruct(black_box(shards2), black_box(validator_count));
					})
				},
			);
		}
	}
}

fn parameterized_criterion() -> Criterion {
	let crit = Criterion::default().sample_size(10).warm_up_time(Duration::from_millis(100));
	crit
}

criterion_group!(
name = plot_paramterized;
config = parameterized_criterion();
targets =
	parameterized::bench_encode_2d,
	parameterized::bench_reconstruct_2d,
	parameterized::bench_encode_fixed_1mb_payload,
	parameterized::bench_reconstruct_fixed_1mb_payload,
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
