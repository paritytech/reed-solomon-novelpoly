use std::result;

use crate::novelpoly::Shard;
use reed_solomon_erasure as reed_solomon_naive;
use reed_solomon_naive::galois_16::ReedSolomon;
pub use reed_solomon_naive::Error;

pub fn to_shards<S: Shard>(payload: &[u8], rs: &ReedSolomon) -> Vec<S> {
	let base_len = payload.len();

	// how many bytes we actually need.
	let needed_shard_len = (base_len + rs.data_shard_count() - 1) / rs.data_shard_count();

	// round up, ing GF(2^16) there are only 2 byte values, so each shard must a multiple of 2
	let needed_shard_len = needed_shard_len + (needed_shard_len & 0x01);

	let shard_len = needed_shard_len;

	let mut shards = vec![S::from(vec![0u8; shard_len]); rs.total_shard_count()];
	for (data_chunk, blank_shard) in payload.chunks(shard_len).zip(&mut shards) {
		// fill the empty shards with the corresponding piece of the payload,
		// zero-padded to fit in the shards.
		let len = std::cmp::min(shard_len, data_chunk.len());
		let blank_shard: &mut [u8] = blank_shard.as_mut();
		blank_shard[..len].copy_from_slice(&data_chunk[..len]);
	}

	shards
}

pub fn rs(validator_count: usize) -> ReedSolomon {
	ReedSolomon::new(validator_count, validator_count - validator_count / 3)
		.expect("this struct is not created with invalid shard number; qed")
}

pub fn encode<S: Shard>(data: &[u8], validator_count: usize) -> result::Result<Vec<S>, Error> {
	let encoder = rs(validator_count);
	let mut shards = to_shards::<S>(data, &encoder);
	encoder.encode(&mut shards).unwrap();
	Ok(shards)
}

pub fn reconstruct<S: Shard>(
	mut received_shards: Vec<Option<S>>,
	validator_count: usize,
) -> result::Result<Vec<u8>, Error> {
	let r = rs(validator_count);

	// Try to reconstruct missing shards
	r.reconstruct_data(&mut received_shards).expect("Sufficient shards must be received. qed");

	let result = received_shards.into_iter().filter_map(|x| x).take(r.data_shard_count()).fold(
		Vec::with_capacity(12 << 20),
		|mut acc, reconstructed_shard| {
			acc.extend_from_slice(AsRef::<[u8]>::as_ref(&reconstructed_shard));
			acc
		},
	);

	Ok(result)
}
