use rand::prelude::*;
use rand::seq::index::IndexVec;
use std::error;
use std::result;
use std::iter;

pub static SMALL_RNG_SEED: [u8; 32] = [
	0, 6, 0xFA, 0, 0x37, 3, 19, 89, 32, 032, 0x37, 0x77, 77, 0b11, 112, 52, 12, 40, 82, 34, 0, 0, 0, 1, 4, 4, 1, 4, 99,
	127, 121, 107,
];

/// Demo test data, generated via `build.rs`.
pub const BYTES: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/rand_data.bin"));

/// Shared target number of shards for simple, quirk turnaround tests:
pub const N_SHARDS: usize = 123;
/// Shared target number of payload size for simple, quirk turnaround tests:
pub const TEST_DATA_CHUNK_SIZE: usize = 1337;

/// Assert the byte ranges derived from the index vec are recovered properly
pub fn assert_recovery(payload: &[u8], reconstructed_payload: &[u8], dropped_indices: IndexVec) {
	assert!(reconstructed_payload.len() >= payload.len());

	dropped_indices.into_iter().for_each(|dropped_idx| {
		let byteoffset = dropped_idx * 2;
		let range = byteoffset..(byteoffset + 2);
		// dropped indices are over `n`, but our data indices are just of length `k * 2`
		if payload.len() >= range.end {
			assert_eq!(
				&payload[range.clone()],
				&reconstructed_payload[range.clone()],
				"Data at bytes {:?} must match:",
				range
			);
		}
	});
}

/// Drop half the shards at the beginning, and half of them at the end.
pub fn deterministic_drop_shards<T: Sized, G: rand::SeedableRng + rand::Rng>(
	codewords: &mut [Option<T>],
	n: usize,
	k: usize,
	_rng: &mut G,
) -> IndexVec {
	let l = codewords.len();
	let mut v = Vec::with_capacity(n - k);
	// k is a power of 2
	let half = (n - k) >> 1;
	for i in 0..half {
		codewords[i] = None;
		v.push(i);
	}
	// if the codewords is shorter than n
	// the remaining ones were
	// already dropped implicitly
	for i in n - half..n {
		if i < l {
			codewords[i] = None;
			v.push(i);
		}
	}
	IndexVec::from(v)
}

pub fn deterministic_drop_shards_clone<T: Sized + Clone>(
	codewords: &[T],
	n: usize,
	k: usize,
) -> (Vec<Option<T>>, IndexVec) {
	let mut rng = SmallRng::from_seed(SMALL_RNG_SEED);
	let mut codewords = codewords.into_iter().map(|x| Some(x.clone())).collect::<Vec<Option<T>>>();
	let idx = deterministic_drop_shards::<T, SmallRng>(&mut codewords, n, k, &mut rng);
	assert!(idx.len() <= n - k);
	(codewords, idx)
}

pub fn drop_random_max<T: Sized + Clone>(
	shards: &mut [Option<T>],
	n: usize,
	k: usize,
	rng: &mut impl rand::Rng,
) -> IndexVec {
	let l = shards.len();
	let already_dropped = n.saturating_sub(l);
	let iv = rand::seq::index::sample(rng, l, n - k - already_dropped);
	assert_eq!(iv.len(), n - k);
	iv.clone().into_iter().for_each(|idx| {
		shards[idx] = None;
	});
	let kept_count = shards.iter().map(Option::is_some).count();
	assert!(kept_count >= k);
	iv
}

pub fn roundtrip<'s, Enc, Recon, S, E>(
	encode: Enc,
	reconstruct: Recon,
	payload: &'s [u8],
	target_shard_count: usize,
) -> result::Result<(), E>
where
	Enc: Fn(&'s [u8], usize) -> result::Result<Vec<S>, E>,
	Recon: Fn(Vec<Option<S>>, usize) -> result::Result<Vec<u8>, E>,
	E: error::Error + Send + Sync + 'static,
	S: Clone + AsRef<[u8]> + AsMut<[[u8; 2]]> + AsRef<[[u8; 2]]> + iter::FromIterator<[u8; 2]> + From<Vec<u8>>,
{
	let v = roundtrip_w_drop_closure::<'s, Enc, Recon, _, SmallRng, S, E>(
		encode,
		reconstruct,
		payload,
		target_shard_count,
		drop_random_max,
	)?;
	Ok(v)
}

pub fn roundtrip_w_drop_closure<'s, Enc, Recon, DropFun, RandGen, S, E>(
	encode: Enc,
	reconstruct: Recon,
	payload: &'s [u8],
	target_shard_count: usize,
	mut drop_rand: DropFun,
) -> result::Result<(), E>
where
	E: error::Error + Send + Sync + 'static,
	S: Clone + AsRef<[u8]> + AsMut<[[u8; 2]]> + AsRef<[[u8; 2]]> + iter::FromIterator<[u8; 2]> + From<Vec<u8>>,
	Enc: Fn(&'s [u8], usize) -> result::Result<Vec<S>, E>,
	Recon: Fn(Vec<Option<S>>, usize) -> result::Result<Vec<u8>, E>,
	DropFun: for<'z> FnMut(&'z mut [Option<S>], usize, usize, &mut RandGen) -> IndexVec,
	RandGen: rand::Rng + rand::SeedableRng<Seed = [u8; 32]>,
{
	let mut rng = <RandGen as rand::SeedableRng>::from_seed(SMALL_RNG_SEED);

	// Construct the shards
	let shards = encode(payload, target_shard_count)?;

	// Make a copy and transform it into option shards arrangement
	// for feeding into reconstruct_shards
	let mut received_shards = shards.into_iter().map(Some).collect::<Vec<Option<S>>>();

	let dropped_indices =
		drop_rand(received_shards.as_mut_slice(), target_shard_count, target_shard_count / 3, &mut rng);

	let recovered_payload = reconstruct(received_shards, target_shard_count)?;

	assert_recovery(&payload[..], &recovered_payload[..], dropped_indices);
	Ok(())
}
