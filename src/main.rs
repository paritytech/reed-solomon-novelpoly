use rs_ec_perf::*;

fn main() {
	roundtrip(status_quo::encode, status_quo::reconstruct, BYTES);
}
