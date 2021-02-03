use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rs_ec_perf::*;

fn bench_status_quo_roundtrip(crit: &mut Criterion) {
	crit.bench_function("status_quo roudtrip", |b| {
		b.iter(|| {
			roundtrip(status_quo::encode, status_quo::reconstruct, black_box(BYTES));
		})
	});
}

fn bench_status_quo_encode(crit: &mut Criterion) {
	crit.bench_function("status_quo encode", |b| {
		b.iter(|| {
			let _ = status_quo::encode(black_box(BYTES));
		})
	});
}

fn adjusted_criterion() -> Criterion {
	let crit = Criterion::default()
		.sample_size(10)
		.warm_up_time(Duration::from_secs(1))
		.measurement_time(Duration::from_secs(60));
	crit
}

criterion_group!(name = benches; config = adjusted_criterion(); targets =  bench_status_quo_roundtrip, bench_status_quo_encode);

criterion_main!(benches);
