#[macro_use(shards)]
extern crate reed_solomon_erasure;

use reed_solomon_erasure::galois_8::ReedSolomon;
// or use the following for Galois 2^16 backend
// use reed_solomon_erasure::galois_16::ReedSolomon;


// we want one message per validator
const N_VALIDATORS: usize = 100;
const DATA_SHARDS: usize = N_VALIDATORS / 3;
const PARITY_SHARDS: usize = N_VALIDATORS - DATA_SHARDS;

// dummy data, everybody should have some good old `cargo-spellcheck` locally installed
// roughly 11 MB
const BYTES: &[u8] = include_bytes!("/home/bernhard/.cargo/bin/cargo-spellcheck");


fn to_shards(payload: &[u8]) -> Vec<Vec<u8>> {

	let base_len = payload.len();

	// how many bytes we actually need.
	let needed_shard_len = base_len / DATA_SHARDS
		+ (base_len % DATA_SHARDS != 0) as usize;

	// round up to next even number
	// (no actual space overhead since we are working in GF(2^16)).
	// XXX I find this statement rather questionable, while we work in GF(2^16)
	// XXX that does not necessarily apply to the chunk to be transfered
	let needed_shard_len = needed_shard_len + (needed_shard_len & 0x01);


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


fn encode(data: &[u8]) -> Vec<Vec<u8>> {
	let mut encoder = rs();
	let mut shards = to_shards(data);
	encoder.encode(&mut shards).unwrap();
	shards
}

fn reconstruct(mut received_shards: Vec<Option<Vec<u8>>>) -> Vec<Vec<u8>> {

	let r = rs();

    // Try to reconstruct missing shards
    r.reconstruct(&mut received_shards).expect("Sufficient shards must be received. qed");

    // Convert back to normal shard arrangement
    let result= received_shards.into_iter().filter_map(|x| x).collect::<Vec<_>>();
	result
}


// struct GF16 {

// }

// impl ff::Field for GF16 {

// }

// fn reconstruct_fft_fast(mut received_shards: Vec<Option<Vec<u8>>>) -> Vec<Vec<u8>> {



// 	ff::Field
// 	fffft::
// }

fn main () {


    // Construct the parity shards
    let encoded = encode(BYTES);

    // Make a copy and transform it into option shards arrangement
    // for feeding into reconstruct_shards
    let mut shards = encoded.clone().into_iter().map(Some).collect::<Vec<_>>();

    // We can remove up to 2 shards, which may be data or parity shards
    shards[0] = None;
    shards[4] = None;


	let result = reconstruct(shards.clone());
    assert!(rs().verify(&result).unwrap());
    assert_eq!(encoded, result);
}
