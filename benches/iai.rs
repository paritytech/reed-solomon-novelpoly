use iai::black_box;
use rs_ec_perf::*;

fn bench_status_quo_roundtrip() {
	roundtrip(status_quo::encode, status_quo::reconstruct, black_box(BYTES), 2000);
}

fn bench_status_quo_encode() {
	let _ = status_quo::encode(black_box(BYTES), 2000);
}

iai::main!(bench_status_quo_roundtrip, bench_status_quo_encode);
