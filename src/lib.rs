pub use novelpoly::{Shard, WrappedShard};
pub use reed_solomon_novelpoly as novelpoly;

#[cfg(feature = "naive")]
pub mod naive;

pub use reed_solomon_tester::{BYTES, N_SHARDS, TEST_DATA_CHUNK_SIZE};

#[cfg(test)]
mod test {

	use super::*;

	#[cfg(feature = "novelpoly")]
	#[test]
	fn novelpoly_roundtrip() -> std::result::Result<(), novelpoly::Error> {
		reed_solomon_tester::roundtrip(
			novelpoly::encode::<WrappedShard>,
			novelpoly::reconstruct::<WrappedShard>,
			&BYTES[..TEST_DATA_CHUNK_SIZE],
			N_SHARDS,
		)
	}

	#[cfg(feature = "novelpoly-with-alt-cxx-impl")]
	fn novelpoly_cxx_roundtrip() -> std::result::Result<(), novelpoly::Error> {
		reed_solomon_tester::roundtrip(
			novelpoly::cxx::encode,
			novelpoly::cxx::reconstruct,
			&BYTES[..TEST_DATA_CHUNK_SIZE],
			N_SHARDS,
		)?;
	}

	#[cfg(feature = "naive")]
	#[test]
	fn naive_roundtrip() -> std::result::Result<(), naive::Error> {
		reed_solomon_tester::roundtrip(
			naive::encode::<WrappedShard>,
			naive::reconstruct::<WrappedShard>,
			&BYTES[..TEST_DATA_CHUNK_SIZE],
			N_SHARDS,
		)
	}
}
