pub static SMALL_RNG_SEED: [u8; 32] = [
	0, 6, 0xFA, 0, 0x37, 3, 19, 89, 32, 032, 0x37, 0x77, 77, 0b11, 112, 52, 12, 40, 82, 34, 0, 0, 0, 1, 4, 4, 1, 4, 99,
	127, 121, 107,
];

mod wrapped_shard;

use rand::prelude::*;
pub use wrapped_shard::*;

pub mod status_quo;

pub mod novel_poly_basis;
#[cfg(feature = "cmp-with-cxx")]
pub mod novel_poly_basis_cxx;

pub const N_VALIDATORS: usize = 32;

pub const BYTES: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/rand_data.bin"));

pub fn drop_random_max(shards: &mut [Option<WrappedShard>], n: usize, k: usize, rng: &mut impl rand::Rng) {
	let iv = rand::seq::index::sample(rng, n, n - k);
	iv.into_iter().for_each(|idx| {
		shards[idx] = None;
	});
}

#[inline(always)]
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
	F: for<'z> FnMut(&'z mut [Option<WrappedShard>], usize, usize, &mut G),
	G: rand::Rng + rand::SeedableRng<Seed = [u8; 32]>,
{
	let mut rng = <G as rand::SeedableRng>::from_seed(SMALL_RNG_SEED);

	// Construct the shards
	let encoded = encode(payload, validator_count);

	// Make a copy and transform it into option shards arrangement
	// for feeding into reconstruct_shards
	let mut shards = encoded.clone().into_iter().map(Some).collect::<Vec<_>>();

	drop_rand(shards.as_mut_slice(), validator_count, validator_count / 3, &mut rng);

	let result = reconstruct(shards, validator_count).expect("reconstruction must work");

	assert!(payload.len() <= result.len());

	// the result might have trailing zeros or non matching trailing data
	itertools::assert_equal(payload[..].iter(), result[0..payload.len()].iter());
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
