use rs_ec_perf::*;

fn main() {
	roundtrip(novel_poly_basis::encode, novel_poly_basis::reconstruct, &BYTES[..32]);
	roundtrip(status_quo::encode, status_quo::reconstruct, &BYTES[..32]);
}
