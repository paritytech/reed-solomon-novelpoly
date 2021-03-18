use color_eyre::Result;
use reed_solomon_performance::WrappedShard;

fn main() -> Result<()> {
	color_eyre::install()?;

	#[allow(unused)]
	use reed_solomon_performance::{BYTES, N_SHARDS, TEST_DATA_CHUNK_SIZE};

	{
		use reed_solomon_performance::novelpoly;
		reed_solomon_tester::roundtrip(
			novelpoly::encode::<WrappedShard>,
			novelpoly::reconstruct::<WrappedShard>,
			&BYTES[..TEST_DATA_CHUNK_SIZE],
			N_SHARDS,
		)?;
	}

	#[cfg(feature = "novelpoly-with-alt-cxx-impl")]
	{
		use reed_solomon_performance::novelpoly;
		reed_solomon_tester::roundtrip(
			novelpoly::cxx::encode,
			novelpoly::cxx::reconstruct,
			&BYTES[..TEST_DATA_CHUNK_SIZE],
			N_SHARDS,
		)?;
	}

	#[cfg(feature = "naive")]
	{
		use reed_solomon_performance::naive;
		reed_solomon_tester::roundtrip(
			naive::encode::<WrappedShard>,
			naive::reconstruct::<WrappedShard>,
			&BYTES[..TEST_DATA_CHUNK_SIZE],
			N_SHARDS,
		)?;
	}

	Ok(())
}
