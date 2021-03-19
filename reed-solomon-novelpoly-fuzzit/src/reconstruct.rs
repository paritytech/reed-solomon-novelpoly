use honggfuzz::fuzz;


use novelpoly::{Shard, WrappedShard};

use arbitrary::*;

use rand::prelude::*;


#[derive(Debug, Clone)]
struct ReconstructionFeed {
	validator_count: usize,
	received: Vec<Option<WrappedShard>>,
}

impl<'a> Arbitrary<'a> for ReconstructionFeed {
	fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {

		let validator_count = u.int_in_range(0_usize..=2200)?;
		let shard_drop_count = u.int_in_range(0_usize..=validator_count)?;

		let n_chunks = validator_count - shard_drop_count;
        let bytes_per_shard = if n_chunks > 0 {
			u.len() / n_chunks
		} else {
			0
		};

		let mut rng = rand_chacha::ChaCha8Rng::from_seed([0u8;32]);
		let iv = rand::seq::index::sample(&mut rng, validator_count, validator_count - n_chunks).into_vec();

		let mut received = (0..validator_count).into_iter()
			.map(|idx| {
				if iv.contains(&idx) {
					None
				} else {
					Some(WrappedShard::new(u.bytes(bytes_per_shard).ok()?.to_vec()))
				}
			})
			.take(validator_count.saturating_sub(1))
			.collect::<Vec<Option<_>>>();

		if u.is_empty() || u.len() > bytes_per_shard / 2 {
			let data = u.bytes(u.len())?.to_vec();
			received.push(Some(WrappedShard::new(data)));
		}

		Ok(Self {
			validator_count,
			received,
		})
	}
}


fn main() {
	// You have full control over the loop but
	// you're supposed to call `fuzz` ad vitam aeternam
	loop {
		// The fuzz macro gives an arbitrary object (see `arbitrary crate`)
		// to a closure-like block of code.
		// For performance reasons, it is recommended that you use the native type
		// `&[u8]` when possible.
		// Here, this slice will contain a "random" quantity of "random" data.
		fuzz!(|feed: ReconstructionFeed| {
            let _ = novelpoly::reconstruct::<WrappedShard>(feed.received, feed.validator_count);
		});
	}
}
