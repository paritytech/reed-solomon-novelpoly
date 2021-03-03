use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};
use rs_ec_perf::*;

/// Create a new testset for a particular RS encoding.
macro_rules! instanciate_test {
	($name:literal, $mp:ident) => {
		pub mod $mp {
			const N_VALS: usize = 2000;

			use super::super::$mp::{encode, reconstruct};
			use super::super::{roundtrip, BYTES};
			use criterion::{black_box, Criterion};

			pub fn bench_roundtrip(crit: &mut Criterion) {
				crit.bench_function(concat!($name, " roudtrip"), |b| {
					b.iter(|| {
						roundtrip(encode, reconstruct, black_box(&BYTES[..]), N_VALS);
					})
				});
			}

			pub fn bench_encode(crit: &mut Criterion) {
				crit.bench_function(concat!($name, " encode"), |b| {
					b.iter(|| {
						let _ = encode(black_box(&BYTES[..]), N_VALS);
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
