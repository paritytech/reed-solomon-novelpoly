
// FIXME use galois_16, but requires AsRef<[u8]>
use reed_solomon_erasure::galois_8::ReedSolomon;


// we want one message per validator
const N_VALIDATORS: usize = 13;
const DATA_SHARDS: usize = N_VALIDATORS / 3;
const PARITY_SHARDS: usize = N_VALIDATORS - DATA_SHARDS;


fn to_shards(payload: &[u8]) -> Vec<Vec<u8>> {

	let base_len = payload.len();

	// how many bytes we actually need.
	let needed_shard_len = base_len / DATA_SHARDS +
		(base_len % DATA_SHARDS != 0) as usize;

	// round up to next even number
	// (no actual space overhead since we are working in GF(2^16)).
	// XXX I find this statement rather questionable, while we work in GF(2^16)
	// XXX that does not necessarily apply to the chunk to be transfered
	// let needed_shard_len = needed_shard_len + (needed_shard_len & 0x01);


	let shard_len = needed_shard_len;

	let mut shards = vec![vec![0u8; shard_len]; N_VALIDATORS];
	for (data_chunk, blank_shard) in payload.chunks(shard_len).zip(&mut shards) {
		// fill the empty shards with the corresponding piece of the payload,
		// zero-padded to fit in the shards.
		let len = std::cmp::min(shard_len, data_chunk.len());
		blank_shard[..len].copy_from_slice(&data_chunk[..len]);
	}

	shards
}

fn rs() -> ReedSolomon {
	ReedSolomon::new(DATA_SHARDS, PARITY_SHARDS)
			.expect("this struct is not created with invalid shard number; qed")
}


pub fn status_quo_encode(data: &[u8]) -> Vec<Vec<u8>> {
	let encoder = rs();
	let mut shards = to_shards(data);
	encoder.encode(&mut shards).unwrap();
	shards
}

pub fn status_quo_reconstruct(mut received_shards: Vec<Option<Vec<u8>>>) -> Option<Vec<u8>> {

	let r = rs();

    // Try to reconstruct missing shards
    r.reconstruct_data(&mut received_shards).expect("Sufficient shards must be received. qed");

    // Convert back to normal shard arrangement
	// let l = received_shards.len();

    // let result_data_shards= received_shards
	// 	.into_iter()
	// 	.filter_map(|x| x)
	// 	.collect::<Vec<Vec<u8>>>();

	let result = received_shards
		.into_iter()
		.filter_map(|x| x)
		.take(DATA_SHARDS)
		.fold(Vec::with_capacity(12<<20), |mut acc, x| {
			acc.extend_from_slice(&x);
			acc
		});

	Some(result)
}


pub fn roundtrip<E,R>(encode: E, reconstruct: R, payload: &[u8])
where
	E: Fn(&[u8]) -> Vec<Vec<u8>>,
	R: Fn(Vec<Option<Vec<u8>>>) -> Option<Vec<u8>>,
{
    // Construct the shards
    let encoded = encode(payload);

    // Make a copy and transform it into option shards arrangement
    // for feeding into reconstruct_shards
    let mut shards = encoded.clone().into_iter().map(Some).collect::<Vec<_>>();

	// Drop 3 shards
	let n = shards.len();
    shards[0] = None;
	// shards[7] = None;
    shards[n.saturating_sub(2)] = None;


	dbg!(shards.len());
	dbg!(shards[1].as_ref().unwrap().len());

	let result = reconstruct(shards).expect("must qork");
    assert_eq!(payload, result);
}


#[cfg(test)]
mod test {
	use super::*;

	const BYTES: &[u8] = include_bytes!("/home/bernhard/.cargo/bin/cargo-spellcheck");

	#[test]
	fn status_quo_roundtrip() {
		roundtrip( status_quo_encode, status_quo_reconstruct, &BYTES[0..100])
	}
}
