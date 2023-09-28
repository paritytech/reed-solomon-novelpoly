#![allow(dead_code)]
#![allow(unused_imports)]

use hongg::fuzz;

use novelpoly::f2e16::*;

use arbitrary::*;

#[derive(Debug, Clone)]
struct FieldMpyParams {
	additive: Additive,
	mpy: Multiplier,
	idx_to_test: usize,
}

#[cfg(target_feature = "avx")]
impl<'a> Arbitrary<'a> for FieldMpyParams {
	fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
		let additive = Additive(u.int_in_range(0..=u16::MAX)?);
		let mpy = Multiplier(u.int_in_range(0..=u16::MAX)?);
		let idx_to_test = u.choose_index(Additive8x::LANE)?;

		Ok(Self { additive, mpy, idx_to_test })
	}
}

fn main() {
	#[cfg(target_feature = "avx")]
	run();

	#[cfg(not(target_feature = "avx"))]
	panic!("Nothing to do for non avx enabled targets")
}

#[cfg(target_feature = "avx")]
fn run() {
	// You have full control over the loop but
	// you're supposed to call `fuzz` ad vitam aeternam
	loop {
		// The fuzz macro gives an arbitrary object (see `arbitrary crate`)
		// to a closure-like block of code.
		// For performance reasons, it is recommended that you use the native type
		// `&[u8]` when possible.
		// Here, this slice will contain a "random" quantity of "random" data.
		fuzz!(|params: FieldMpyParams| {
			let FieldMpyParams { idx_to_test, additive, mpy } = params;
			let values = [additive; 8];
			let values8x = Additive8x::from(values);
			let res_faster8 = values8x.mul(mpy);
			let res_plain = values[idx_to_test].mul(mpy);
			assert_eq!(res_plain, Additive8x::unpack(&res_faster8)[idx_to_test]);
		});
	}
}
