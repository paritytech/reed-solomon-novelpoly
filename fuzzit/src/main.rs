use honggfuzz::fuzz;

use rs::*;

use arbitrary::*;

#[derive(Debug, Clone, Copy)]
struct ValidatorCount(usize);

impl std::ops::Deref for ValidatorCount {
	type Target = usize;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}

impl<'a> Arbitrary<'a> for ValidatorCount {
	fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
		let data = u16::arbitrary(u)?;
		if data > 2200 {
			Err(arbitrary::Error::IncorrectFormat)
		} else {
			Ok(Self(data as usize))
		}
	}
}

#[derive(Debug, Clone, Copy, Arbitrary)]
struct Feed<'a> {
	validator_count: ValidatorCount,
	data: &'a [u8],
}

fn main() {
	novel_poly_basis::setup();

	// You have full control over the loop but
	// you're supposed to call `fuzz` ad vitam aeternam
	loop {
		// The fuzz macro gives an arbitrary object (see `arbitrary crate`)
		// to a closure-like block of code.
		// For performance reasons, it is recommended that you use the native type
		// `&[u8]` when possible.
		// Here, this slice will contain a "random" quantity of "random" data.
		fuzz!(|feed: Feed| {
			let _ =
				roundtrip(novel_poly_basis::encode, novel_poly_basis::reconstruct, feed.data, *feed.validator_count);
		});
	}
}
