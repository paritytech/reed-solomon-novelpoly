mod wrapped_shard;
pub use wrapped_shard::*;

pub mod status_quo;

pub mod novel_poly_basis;

// we want one message per validator, so this is the total number of shards that we should own
// after
const N_VALIDATORS: usize = 16; //256;
const DATA_SHARDS: usize = 4; // N_VALIDATORS / 3;
const PARITY_SHARDS: usize = N_VALIDATORS - DATA_SHARDS;

pub const BYTES: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/rand_data.bin"));

pub fn roundtrip<E, R>(encode: E, reconstruct: R, payload: &[u8])
where
	E: Fn(&[u8]) -> Vec<WrappedShard>,
	R: Fn(Vec<Option<WrappedShard>>) -> Option<Vec<u8>>,
{
	// Construct the shards
	let encoded = encode(payload);

	// Make a copy and transform it into option shards arrangement
	// for feeding into reconstruct_shards
	let mut shards = encoded.clone().into_iter().map(Some).collect::<Vec<_>>();

	// Drop 3 shards
	let mut rng = rand::thread_rng();

	// randomly lose `2/3 - eps` of the messages
	let iv = rand::seq::index::sample(&mut rng, N_VALIDATORS, (N_VALIDATORS << 1) / 3);
	iv.into_iter().for_each(|idx| {
		shards[idx] = None;
	});

	let result = reconstruct(shards).expect("reconstruction must work");

	// the result might have trailing zeros
	assert_eq!(&payload[..], &result[0..payload.len()]);
}

#[cfg(test)]
mod test {
	use super::*;

	#[test]
	fn status_quo_roundtrip() {
		roundtrip(status_quo::encode, status_quo::reconstruct, &BYTES[0..32])
	}

	#[test]
	fn novel_poly_basis_roundtrip() {
		roundtrip(novel_poly_basis::encode, novel_poly_basis::reconstruct, &BYTES[0..32])
	}
}
