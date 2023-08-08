decl_field!("f256", u8, u16, 8, gen = 0x1D, cantor = [1, 214, 152, 146, 86, 200, 88, 230]);

include!("inc_log_mul.rs");

#[cfg(table_bootstrap_complete)]
include!("inc_afft.rs");

#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Additive8x(());

impl Additive8x {
	pub const LANE: usize = 8;

	pub fn zero() -> Self {
		Self(())
	}

	pub fn mul(&self, other: Multiplier) -> Self {
		unreachable!()
	}

	pub fn load(_data: &[u8]) -> Self {
		unreachable!()
	}

	pub fn copy_to(dest: &mut [Additive]) {
		unreachable!()
	}
}
