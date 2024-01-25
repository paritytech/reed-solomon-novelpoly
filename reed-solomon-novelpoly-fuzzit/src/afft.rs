#![allow(dead_code)]
#![allow(unused_imports)]
use hongg::fuzz;

use novelpoly::f2e16::*;

use arbitrary::*;

use rand::prelude::*;

#[derive(Debug, Clone)]
struct AfftParams {
	k: usize,
	f: usize,
	shift: usize,
}

impl<'a> Arbitrary<'a> for AfftParams {
	fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
		let k = 1 << u.int_in_range(4..=12)?;
		let f = u.int_in_range(2..=4)?;
		let shift = k * u.choose_index(f)?;

		Ok(Self { k, f, shift })
	}
}

fn main() {
	#[cfg(all(target_feature = "avx", feature = "avx"))]
	run();

	#[cfg(not(target_feature = "avx"))]
	panic!("Nothing to do for non avx enabled targets")
}

#[cfg(all(target_feature = "avx", feature = "avx"))]
fn run() {
	// You have full control over the loop but
	// you're supposed to call `fuzz` ad vitam aeternam

	loop {
		// The fuzz macro gives an arbitrary object (see `arbitrary crate`)
		// to a closure-like block of code.
		// For performance reasons, it is recommended that you use the native type
		// `&[u8]` when possible.
		// Here, this slice will contain a "random" quantity of "random" data.
		fuzz!(|params: AfftParams| {
			let AfftParams { k, shift, f: _ } = params;

			let mut data_plain = rstester::gen_plain::<SmallRng>(k);
			let mut data_faster8 = rstester::gen_faster8_from_plain(&data_plain);

			// k := size , so we need to have shift as steps of k
			unsafe { &AFFT }.afft(&mut data_plain, k, shift);
			unsafe { &AFFT }.afft_faster8(&mut data_faster8[..], k, shift);

			rstester::assert_plain_eq_faster8(data_plain, data_faster8);
		});
	}
}
