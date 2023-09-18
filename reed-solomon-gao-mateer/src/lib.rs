#![forbid(unused_crate_dependencies)]

pub mod errors;
pub use errors::*;

pub mod util;
pub use util::*;


pub mod gao_mateer;

pub mod shard;
pub use self::shard::Shard;

pub mod wrapped_shard;
pub use self::wrapped_shard::WrappedShard;

#[cfg(test)]
mod test {
	use super::*;
	use reed_solomon_tester::{roundtrip, BYTES, N_SHARDS};

	#[cfg(feature = "naive")]
	#[test]
	fn status_quo_roundtrip() -> Result<()> {
		roundtrip(status_quo::encode::<WrappedShard>, status_quo::reconstruct::<WrappedShard>, &BYTES[..1337], N_SHARDS)
	}

	#[test]
	fn novel_poly_basis_roundtrip() -> Result<()> {
		roundtrip(
			novel_poly_basis::encode::<WrappedShard>,
			novel_poly_basis::reconstruct::<WrappedShard>,
			&BYTES[..1337],
			N_SHARDS,
		)
	}
}
