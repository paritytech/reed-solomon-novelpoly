#[cfg(table_bootstrap_complete)]
use super::*;

decl_field!(
	"f2e16",
	u16,
	u32,
	16,
	gen = 0x2D,
	cantor =
		[1_u16, 44234, 15374, 5694, 50562, 60718, 37196, 16402, 27800, 4312, 27250, 47360, 64952, 64308, 65336, 39198]
);

include!("inc_log_mul.rs");

#[cfg(table_bootstrap_complete)]
include!("inc_afft.rs");

#[cfg(table_bootstrap_complete)]
include!("inc_encode.rs");

#[cfg(table_bootstrap_complete)]
include!("inc_reconstruct.rs");

use core::ops::{BitXor, BitXorAssign};
use core::arch::x86_64::*;

#[allow(non_camel_case_types)]
pub type u16x8 = __m128i;

#[allow(non_camel_case_types)]
pub type u32x8 = __m256i;

#[repr(transparent)]
#[derive(Debug, Clone, Copy)]
pub struct Additive8x(pub u16x8);

impl PartialEq<Self> for Additive8x {
	fn eq(&self, other: &Self) -> bool {
		unsafe {
			// set 0xFFFF for equality, 0x000 for inequality
			let eq = _mm_cmpeq_epi16(self.0, other.0);
			let mask = splat_u16x8(ONEMASK);
			// 0xFFFF for INEQUALITY, 0x0000 for EQUALITY
			let inverted = _mm_xor_si128(eq, mask);
			// if all bits are 0s, things are equal
			_mm_testz_si128(inverted, mask) == 1
		}
	}
}

impl Eq for Additive8x {}

impl BitXor for Additive8x {
	type Output = Self;

	fn bitxor(self, rhs: Self) -> Self::Output {
		Self(unsafe { _mm_xor_si128(self.0, rhs.0) })
	}
}

impl BitXorAssign for Additive8x {
	fn bitxor_assign(&mut self, rhs: Self) {
		unsafe {
			self.0 = _mm_xor_si128(self.0, rhs.0);
		}
	}
}

fn splat_u16x8(v: u16) -> u16x8 {
	unsafe { _mm_set1_epi16(v as i16) }
}

fn splat_u32x8(v: u32) -> u32x8 {
	unsafe { _mm256_set1_epi32(v as i32) }
}

unsafe fn unpack_u32x8(v: u32x8) -> [u32; 8] {
	#[repr(C, align(64))]
	struct Aligned;

	#[repr(C, align(64))]
	struct Stack {
		_alignment: [Aligned; 0],
		data: [u32; 8],
	}

	let mut dest = Stack { _alignment: [], data: [0u32; 8] };

	_mm256_store_si256(core::intrinsics::transmute((&mut dest.data) as *mut [u32; 8]), v);

	dest.data
}

unsafe fn unpack_u16x8(v: u16x8) -> [u16; 8] {
	#[repr(C, align(32))]
	struct Aligned;

	#[repr(C, align(32))]
	struct Stack {
		_alignment: [Aligned; 0],
		data: [u16; 8],
	}

	let mut dest = Stack { _alignment: [], data: [0u16; 8] };

	_mm_store_si128(core::intrinsics::transmute((&mut dest.data) as *mut [u16; 8]), v);

	dest.data
}

unsafe fn load_u16x8(src: &[u16], idx: [u16; 8]) -> u16x8 {
	_mm_set_epi16(
		src[idx[0] as usize] as i16,
		src[idx[1] as usize] as i16,
		src[idx[2] as usize] as i16,
		src[idx[3] as usize] as i16,
		src[idx[4] as usize] as i16,
		src[idx[5] as usize] as i16,
		src[idx[6] as usize] as i16,
		src[idx[7] as usize] as i16,
	)
}

impl Additive8x {
	const LANE: usize = 8;

	pub fn zero() -> Self {
		Self(unsafe { _mm_setzero_si128() })
	}

	#[cfg(table_bootstrap_complete)]
	pub fn mul(&self, other: Multiplier) -> Self {
		// let log = (LOG_TABLE[self.0 as usize] as Wide) + other.0 as Wide;
		// let offset = (log & ONEMASK as Wide) + (log >> FIELD_BITS);
		// Additive(EXP_TABLE[offset as usize])

		unsafe {
			// spread the multiplier
			let other = splat_u32x8(other.0 as Wide);
			unpack_u32x8(other);

			// get the indices
			let idx = unpack_u16x8(self.0);
			// load from the table based on the indices
			let logtable = load_u16x8(LOG_TABLE.as_slice(), idx);
			unpack_u16x8(logtable);

			let logtable = _mm256_cvtepu16_epi32(logtable);

			let log = _mm256_add_epi16(logtable, other);
			unpack_u32x8(log);

			// (log & ONEMASK) + (log >> shift)
			let onemasks = splat_u32x8(ONEMASK as Wide);
			let offset_sum_left = _mm256_and_si256(log, onemasks);
			let offset_sum_right = _mm256_srli_epi32(log, FIELD_BITS as i32);
			let offset = _mm256_add_epi32(offset_sum_left, offset_sum_right);

			// lut;
			// u32 to u16 with 0x0000FFFF mask per element
			let offset = _mm256_packus_epi32(offset, offset);
			let offset = _mm256_castsi256_si128(offset);
			let offset = unpack_u16x8(offset);
			let res = load_u16x8(EXP_TABLE.as_slice(), offset);
			Self(res)
		}
	}

	pub fn load(data: &[Additive]) -> Self {
		assert_eq!(data.len(), 8);
		let v = unsafe { _mm_loadu_si128(core::intrinsics::transmute(data.as_ptr())) };
		Self(v)
	}

	pub fn copy_to_slice(&self, dest: &mut [Additive]) {
		assert!(dest.len() >= 8);
		todo!()
	}

	pub fn copy_from_slice(&mut self, src: &[Additive]) {
		*self = Self::load(&src[..8]);
	}
}

impl From<[Additive; 8]> for Additive8x {
	fn from(a: [Additive; 8]) -> Self {
		Self::load(&a[..])
	}
}
impl Into<[Additive; 8]> for Additive8x {
	fn into(self) -> [Additive; 8] {
		unsafe { std::intrinsics::transmute(unpack_u16x8(self.0)) }
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn identical() {
		let m = Multiplier(12048);
		let v = [Additive(2345); 8];

		let b = v[0].clone();
		let b = b.mul(m);

		let a = Additive8x::from(v);
		let a = a.mul(m);
		let a: [Additive; 8] = a.into();
		assert_eq!(a[0], b);
	}
}
