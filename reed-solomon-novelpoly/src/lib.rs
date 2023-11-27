#![forbid(unused_crate_dependencies)]

pub mod errors;
pub use errors::*;

pub mod util;
pub use util::*;

pub mod field;
#[cfg(feature = "f256")]
pub use self::field::f256;
pub use self::field::f2e16;

#[cfg(target_feature = "avx")]
pub use self::field::faster8;

mod novel_poly_basis;
pub use self::novel_poly_basis::*;

pub mod shard;
pub use self::shard::Shard;

pub mod wrapped_shard;
pub use self::wrapped_shard::WrappedShard;

#[cfg(feature = "with-alt-cxx-impl")]
pub mod cxx;

#[cfg(test)]
mod test {
	use crate::f2e16::encode_sub;

	use super::*;
	use reed_solomon_tester::{roundtrip, BYTES, N_SHARDS};

	#[test]
	fn novel_poly_basis_roundtrip() -> Result<()> {
		roundtrip(
			novel_poly_basis::encode::<WrappedShard>,
			novel_poly_basis::reconstruct::<WrappedShard>,
			&BYTES[..1337],
			N_SHARDS,
		)
	}

	/// Showcase the systematic nature of the algorithm.
	#[test]
	fn systematic_for_sure() {
		let bytes = [1_u8, 2, 3, 4];
		let bytes = bytes.as_slice();

		let shards = encode_sub(bytes, 8, 4).unwrap();
		let x = WrappedShard::from_iter(shards.iter().map(|x| x.0.to_be_bytes()));
		assert_eq!(&x.into_inner()[..bytes.len()], bytes);
	}
}
