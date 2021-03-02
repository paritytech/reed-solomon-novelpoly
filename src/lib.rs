mod wrapped_shard;
use rand::seq::index::{IndexVec, IndexVecIntoIter};
pub use wrapped_shard::*;

pub mod status_quo;

pub mod novel_poly_basis;
pub mod novel_poly_basis_cxx;

// we want one message per validator, so this is the total number of shards that we should own
// after
pub const N_VALIDATORS: usize = 32; //256;
pub const DATA_SHARDS: usize = 4; // N_VALIDATORS / 3;
pub const PARITY_SHARDS: usize = N_VALIDATORS - DATA_SHARDS;

pub const BYTES: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/rand_data.bin"));

pub fn roundtrip<E, R>(encode: E, reconstruct: R, payload: &[u8])
where
E: for<'r> Fn(&'r [u8]) -> Vec<WrappedShard>,
R: Fn(Vec<Option<WrappedShard>>) -> Option<Vec<u8>>,
{

	// Drop 3 shards
	let mut rng = rand::thread_rng();
	let drop_random = move |shards: &mut [Option<WrappedShard>], n: usize, k: usize| -> IndexVec {
		let iv = rand::seq::index::sample(&mut rng, n, (n-k));
		iv
	};
	roundtrip_w_drop_closure::<E,R,_>(encode, reconstruct, payload, drop_random)
}

pub fn roundtrip_w_drop_closure<E, R, F>(encode: E, reconstruct: R, payload: &[u8], mut drop_rand: F )
where
	E: for<'r> Fn(&'r [u8]) -> Vec<WrappedShard>,
	R: Fn(Vec<Option<WrappedShard>>) -> Option<Vec<u8>>,
	F: for<'z> FnMut(&'z mut [Option<WrappedShard>], usize, usize) -> IndexVec,
{

	assert!(DATA_SHARDS >= payload.len() / 2);
	// Construct the shards
	let encoded = encode(payload);

	// Make a copy and transform it into option shards arrangement
	// for feeding into reconstruct_shards
	let mut shards = encoded.clone().into_iter().map(Some).collect::<Vec<_>>();

	let iv = drop_rand(shards.as_mut_slice(), N_VALIDATORS, DATA_SHARDS);
	iv.into_iter().for_each(|idx| { shards[idx] = None; });

	let result = reconstruct(shards).expect("reconstruction must work");

	// the result might have trailing zeros or non matching trailing data
	assert_eq!(&payload[..], &result[0..payload.len()]);
}

#[cfg(test)]
mod test {
	use super::*;

	#[test]
	fn status_quo_roundtrip() {
		roundtrip(status_quo::encode, status_quo::reconstruct, &BYTES[0..(DATA_SHARDS<<1)])
	}

	#[test]
	fn novel_poly_basis_roundtrip() {
		roundtrip(novel_poly_basis::encode, novel_poly_basis::reconstruct, &BYTES[0..(DATA_SHARDS<<1)])
	}
}
