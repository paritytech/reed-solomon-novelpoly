use iai::black_box;
use rs_ec_perf::*;


// dummy data, everybody should have some good old `cargo-spellcheck` locally installed
// roughly 11 MB
const BYTES: &[u8] = include_bytes!("/home/bernhard/.cargo/bin/cargo-spellcheck");

fn bench_status_quo_roundtrip() {
    roundtrip(status_quo_encode, status_quo_reconstruct, black_box(&BYTES[0..512]));
}

fn bench_status_quo() {
    let _ = status_quo_encode(black_box(&BYTES[0..512]));
}

iai::main!(bench_status_quo_roundtrip, bench_status_quo);
