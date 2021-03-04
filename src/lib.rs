pub static SMALL_RNG_SEED: [u8; 32] = [
	0, 6, 0xFA, 0, 0x37, 3, 19, 89, 32, 032, 0x37, 0x77, 77, 0b11, 112, 52, 12, 40, 82, 34, 0, 0, 0, 1, 4, 4, 1, 4, 99,
	127, 121, 107,
];

mod wrapped_shard;

use rand::prelude::*;
use rand::seq::index::IndexVec;

pub use wrapped_shard::*;

pub mod status_quo;

pub mod novel_poly_basis;
#[cfg(feature = "cmp-with-cxx")]
pub mod novel_poly_basis_cxx;

pub const N_VALIDATORS: usize = 128;

pub const BYTES: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/rand_data.bin"));

/// Assert the byte ranges derived from the index vec are recovered properly
pub fn assert_recovery(payload: &[u8], reconstructed_payload: &[u8], dropped_indices: IndexVec) {
	assert!(reconstructed_payload.len() >= payload.len());

	dropped_indices.into_iter().for_each(|dropped_idx| {
		let byteoffset = dropped_idx * 2;
		let range = byteoffset..(byteoffset+2);
		// dropped indices are over `n`, but our data indices are just of length `k * 2`
		if payload.len() >= range.end {
			assert_eq!(&payload[range.clone()], &reconstructed_payload[range.clone()],
				"Data at bytes {:?} must match:", range);
		}
	});
}

pub fn drop_random_max(shards: &mut [Option<WrappedShard>], n: usize, k: usize, rng: &mut impl rand::Rng) -> IndexVec {
	let l = shards.len();
	let already_dropped = n.saturating_sub(l);
	let iv = rand::seq::index::sample(rng, l, n - k - already_dropped);
	assert_eq!(iv.len(), n-k);
	iv.clone().into_iter().for_each(|idx| {
		shards[idx] = None;
	});
	let kept_count = shards.iter().map(Option::is_some).count();
	assert!(kept_count >= k);
	iv
}

pub fn roundtrip<E, R>(encode: E, reconstruct: R, payload: &[u8], validator_count: usize)
where
	E: for<'r> Fn(&'r [u8], usize) -> Vec<WrappedShard>,
	R: Fn(Vec<Option<WrappedShard>>, usize) -> Option<Vec<u8>>,
{
	roundtrip_w_drop_closure::<E, R, _, SmallRng>(encode, reconstruct, payload, validator_count, drop_random_max)
}

pub fn roundtrip_w_drop_closure<E, R, F, G>(
	encode: E,
	reconstruct: R,
	payload: &[u8],
	validator_count: usize,
	mut drop_rand: F,
) where
	E: for<'r> Fn(&'r [u8], usize) -> Vec<WrappedShard>,
	R: Fn(Vec<Option<WrappedShard>>, usize) -> Option<Vec<u8>>,
	F: for<'z> FnMut(&'z mut [Option<WrappedShard>], usize, usize, &mut G) -> IndexVec,
	G: rand::Rng + rand::SeedableRng<Seed = [u8; 32]>,
{
	let mut rng = <G as rand::SeedableRng>::from_seed(SMALL_RNG_SEED);

	// Construct the shards
	let shards = encode(payload, validator_count);

	// Make a copy and transform it into option shards arrangement
	// for feeding into reconstruct_shards
	let mut received_shards = shards
		.into_iter()
		.map(Some)
		.collect::<Vec<Option<WrappedShard>>>();

	let dropped_indices = drop_rand(received_shards.as_mut_slice(), validator_count, validator_count / 3 , &mut rng);

	let recovered_payload = reconstruct(received_shards, validator_count).expect("reconstruction must work");

	assert_recovery(&payload[..], &recovered_payload[..], dropped_indices);
}

#[cfg(test)]
mod test {
	use super::*;

	#[test]
	fn status_quo_roundtrip() {
		roundtrip(status_quo::encode, status_quo::reconstruct, &BYTES[..1337], N_VALIDATORS)
	}

	#[test]
	fn novel_poly_basis_roundtrip() {
		roundtrip(novel_poly_basis::encode, novel_poly_basis::reconstruct, &BYTES[..1337], N_VALIDATORS)
	}
}
