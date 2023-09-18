decl_field!("f256", u8, u16, 8, gen = 0x1D, cantor = [1, 214, 152, 146, 86, 200, 88, 230]);

include!("inc_log_mul.rs");

#[cfg(table_bootstrap_complete)]
include!("inc_afft.rs");

impl std::fmt::Display for Additive {
	fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		write!(f, "{:02x}", self.0)
	}
}

impl std::fmt::Debug for Additive {
	fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		write!(f, "{:02x}", self.0)
	}
}

/// A bare bones gap fill dummy implementation.
#[repr(transparent)]
#[derive(Clone, Copy)]
pub struct Additive8x([Additive; Self::LANE]);

impl From<[Additive; Self::LANE]> for Additive8x {
	fn from(x: [Additive; Self::LANE]) -> Self {
		Self(x)
	}
}

impl Additive8x {
	pub const LANE: usize = 8;

	pub fn zero() -> Self {
		Self([Additive::zero(); 8])
	}

	#[cfg(table_bootstrap_complete)]
	pub fn mul(&self, other: Multiplier) -> Self {
		let mut copy = self.0.clone();
		copy.iter_mut().for_each(move |x| {
			*x = x.clone().mul(other);
		});
		Self(copy)
	}

	pub fn load(data: &[Additive]) -> Self {
		let mut dest = [Additive::zero(); Additive8x::LANE];
		dest.copy_from_slice(data);
		Self(dest)
	}

	pub fn copy_to_slice(&self, dest: &mut [Additive]) {
		dest.copy_from_slice(self.0.as_slice());
	}

	pub fn unpack(&self) -> [Additive; Self::LANE] {
		self.0.clone()
	}
}

impl std::fmt::Debug for Additive8x {
	fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		let unpacked = &self.0;
		write!(
			f,
			"{}, {}, {}, {}, {}, {}, {}, {}",
			unpacked[0], unpacked[1], unpacked[2], unpacked[3], unpacked[4], unpacked[5], unpacked[6], unpacked[7]
		)
	}
}

use std::ops::{BitXor, BitXorAssign};
use std::cmp::{PartialEq, Eq};

impl PartialEq<Self> for Additive8x {
	/// This includes a fallacy, when migrating from [Additive;8*y] to [Additive8x;y]
	#[inline(always)]
	fn eq(&self, other: &Self) -> bool {
		self.0.iter().zip(other.0.iter()).find(|(a, b)| a != b).is_none()
	}
}

impl Eq for Additive8x {}

impl BitXor for Additive8x {
	type Output = Self;

	#[inline(always)]
	fn bitxor(self, rhs: Self) -> Self::Output {
		let mut x = Self::zero();
		for i in 0..Self::LANE {
			x.0[i] = self.0[i] ^ rhs.0[i];
		}
		x
	}
}

impl BitXorAssign for Additive8x {
	fn bitxor_assign(&mut self, rhs: Self) {
		for i in 0..Self::LANE {
			self.0[i] ^= rhs.0[i];
		}
	}
}
