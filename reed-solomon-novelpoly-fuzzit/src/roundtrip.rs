use honggfuzz::fuzz;

use novelpoly::WrappedShard;

use arbitrary::*;


#[derive(Debug, Clone, Copy)]
struct RoundtripFeed<'a> {
	validator_count: usize,
	data: &'a [u8],
}



impl<'a> Arbitrary<'a> for RoundtripFeed<'a> {
	fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
		let validator_count = u.int_in_range(0..=2200)?;
		Ok(Self {
			validator_count,
			data: u.bytes(u.len())?,
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
		fuzz!(|feed: RoundtripFeed| {
            let _ = rstester::roundtrip(
                novelpoly::encode::<WrappedShard>,
                novelpoly::reconstruct::<WrappedShard>,
                feed.data,
                feed.validator_count,
            );
		});
	}
}
