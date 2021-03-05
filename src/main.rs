use rs_ec_perf::*;

fn main() -> Result<()> {
	#[cfg(feature = "cmp-with-cxx")]
	{
		// roundtrip(novel_poly_basis_cxx::encode, novel_poly_basis_cxx::reconstruct, &BYTES[..DATA_SHARDS * 2], DATA_SHARDS)?;
	}

	roundtrip(novel_poly_basis::encode, novel_poly_basis::reconstruct, &BYTES[..], N_VALIDATORS)?;

	#[cfg(feature = "status-quo")]
	{
		roundtrip(status_quo::encode, status_quo::reconstruct, &BYTES[..], N_VALIDATORS)?;
	}

	Ok(())
}
