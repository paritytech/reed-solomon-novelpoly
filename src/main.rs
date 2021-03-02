use rs_ec_perf::*;

fn main() {
	// roundtrip(novel_poly_basis_cxx::encode, novel_poly_basis_cxx::reconstruct, &BYTES[..DATA_SHARDS * 2]);
	roundtrip(novel_poly_basis::encode, novel_poly_basis::reconstruct, &BYTES[..DATA_SHARDS * 2]);
	roundtrip(status_quo::encode, status_quo::reconstruct, &BYTES[..DATA_SHARDS * 2]);
}
