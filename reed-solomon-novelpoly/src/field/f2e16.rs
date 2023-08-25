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

impl std::fmt::Display for Additive {
	fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		write!(f, "{}", self.0)
	}
}

impl std::fmt::Debug for Additive {
	fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		write!(f, "{}", self.0)
	}
}

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

#[inline(always)]
fn splat_u16x8(v: u16) -> u16x8 {
	unsafe { _mm_set1_epi16(v as i16) }
}

#[inline(always)]
fn splat_u32x8(v: u32) -> u32x8 {
	unsafe { _mm256_set1_epi32(v as i32) }
}

#[inline(always)]
fn unpack_u32x8(v: u32x8) -> [u32; 8] {
	#[repr(C, align(32))]
	struct Aligned;

	#[repr(C, align(32))]
	struct Stack {
		_alignment: [Aligned; 0],
		data: [u32; 8],
	}

	let mut dest = Stack { _alignment: [], data: [0u32; 8] };

	unsafe {
		_mm256_store_si256(core::intrinsics::transmute((&mut dest.data) as *mut [u32; 8]), v);
	}

	dest.data
}

#[inline(always)]
fn unpack_u16x8(v: u16x8) -> [u16; 8] {
	#[repr(C, align(16))]
	struct Aligned;

	#[repr(C, align(16))]
	struct Stack {
		_alignment: [Aligned; 0],
		data: [u16; 8],
	}

	let mut dest = Stack { _alignment: [], data: [0u16; 8] };

	unsafe {
		_mm_store_si128(core::intrinsics::transmute((&mut dest.data) as *mut [u16; 8]), v);
	}
	dest.data
}

#[inline(always)]
fn load_u16x8(src: &[u16], idx: [u16; 8]) -> u16x8 {
	unsafe {
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
}

impl Additive8x {
	pub const LANE: usize = 8;

	#[inline(always)]
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

	#[inline(always)]
	pub fn load(data: &[Additive]) -> Self {
		assert!(data.len() >= Self::LANE);
		let v = unsafe { _mm_loadu_si128(core::intrinsics::transmute(data.as_ptr())) };
		Self(v)
	}

	// shift the data b y `shift` u16 positions, then load the
	#[inline(always)]
	pub fn override_partial_continuous(&mut self, shift: usize, data: &[Additive]) {
		assert!(data.len() <= Additive8x::LANE);
		let mut value = [Additive::zero(); Additive8x::LANE];
		let mut mask = [Additive::zero(); Additive8x::LANE];
		let max_shift = data.len() + shift;
		let end = if max_shift > Additive8x::LANE {
			eprintln!("Dropping elements: shift={shift:?}");
			Additive8x::LANE.saturating_sub(shift)
		} else {
			data.len()
		};
		for i in 0..end {
			value[i + shift] = data[i];
			mask[i + shift] = Additive(0xFFFF);
		}
		let value = Self::load(&value[..]).0;
		let mask = Self::load(&mask[..]).0;
		let v = unsafe {
			let masked_value = _mm_and_si128(mask, value);
			let masked_self = _mm_andnot_si128(mask, self.0);
			_mm_or_si128(masked_self, masked_value)
		};
		self.0 = v;
	}

	#[inline(always)]
	pub fn copy_to_slice(&self, dest: &mut [Additive]) {
		assert!(dest.len() >= Self::LANE);
		unsafe {
			_mm_storeu_si128(core::intrinsics::transmute(dest.as_ptr()), self.0);
		}
	}

	#[inline(always)]
	pub fn copy_from_slice(&mut self, src: &[Additive]) {
		assert!(src.len() >= Self::LANE);
		*self = Self::load(&src[..Self::LANE]);
	}

	#[inline(always)]
	pub fn unpack(&self) -> [Additive; Self::LANE] {
		unsafe { std::mem::transmute::<_, [Additive; Self::LANE]>(unpack_u16x8(self.0)) }
	}

	#[inline(always)]
	pub fn splat(value: Additive) -> Self {
		Self(splat_u16x8(value.0))
	}
}

impl From<[Additive; Additive8x::LANE]> for Additive8x {
	fn from(a: [Additive; Additive8x::LANE]) -> Self {
		Self::load(&a[..])
	}
}
impl Into<[Additive; Additive8x::LANE]> for Additive8x {
	fn into(self) -> [Additive; Additive8x::LANE] {
		unsafe { std::intrinsics::transmute(unpack_u16x8(self.0)) }
	}
}

pub fn convert_to_faster8(data: &[Additive], dest: &mut [Additive8x]) {
	assert_eq!(data.len() % Additive8x::LANE, 0);
	assert!(data.len() <= dest.len() * Additive8x::LANE);
	(0..(data.len() + Additive8x::LANE - 1) / Additive8x::LANE)
		.into_iter()
		.for_each(|i| dest[i] = Additive8x::load(&data[(i * Additive8x::LANE)..][..Additive8x::LANE]));
}

pub fn convert_from_faster8(data: &[Additive8x], dest: &mut [Additive]) {
	data.iter()
		.enumerate()
		.for_each(|(idx, x8)| x8.copy_to_slice(&mut dest[idx * Additive8x::LANE..][..Additive8x::LANE]));
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn identical_copy_in_copy_out() {
		assert!(cfg!(target_feature = "avx2"), "Tests are meaningless without avx2 target feature");

		// setup some data
		let v = [Additive(2345); 8];

		// load it
		let mut a = Additive8x::from(v);

		// construct some dummy scratch space
		let mut scratch = vec![Additive(0xFFFF_u16); 8];

		// copy date to scratch space
		a.copy_to_slice(&mut scratch);

		assert_eq!(scratch, v);

		// copy in _different_ values , although

		let v = [Additive(1177); 8];
		a.copy_from_slice(&v);

		let array = Into::<[Additive; Additive8x::LANE]>::into(a);
		for i in 0..Additive8x::LANE {
			assert_eq!(array[i], Additive(1177));
		}
	}

	#[test]
	fn identical_slice_elements() {
		assert!(cfg!(target_feature = "avx2"), "Tests are meaningless without avx2 target feature");

		let m = Multiplier(12048);
		let v = [Additive(2345); 8];

		let b = v[0].clone();
		let b = b.mul(m);

		let a = Additive8x::from(v);
		let a = a.mul(m);
		let a: [Additive; 8] = a.into();
		assert_eq!(a[0], b);
		assert_eq!(a[1], b);
		assert_eq!(a[2], b);
		assert_eq!(a[3], b);
		assert_eq!(a[4], b);
		assert_eq!(a[5], b);
		assert_eq!(a[6], b);
		assert_eq!(a[7], b);
	}

	#[test]
	fn identical_mul() {
		assert!(cfg!(target_feature = "avx2"), "Tests are meaningless without avx2 target feature");

		let m = Multiplier(7);
		let faster8 = Additive8x::from([Additive(2345); 8]);
		let single = Additive(2345);

		let res_faster8 = faster8.mul(m);
		let res_faster8 = Into::<[Additive; Additive8x::LANE]>::into(res_faster8);
		let res = single.mul(m);

		for i in 0..Additive8x::LANE {
			assert_eq!(res_faster8[i], res);
		}
	}

	#[test]
	fn identical_splat_u16x8() {
		assert!(cfg!(target_feature = "avx2"), "Tests are meaningless without avx2 target feature");
		let value = 0xABCD;
		let eight_values = splat_u16x8(value);
		let eight_values = unpack_u16x8(eight_values);
		for i in 0..(eight_values.len()) {
			assert_eq!(eight_values[i], value);
		}
	}

	#[test]
	fn identical_splat_u32x8() {
		assert!(cfg!(target_feature = "avx2"), "Tests are meaningless without avx2 target feature");
		let value = 0x09ABCDEF;
		let eight_values = splat_u32x8(value);
		let eight_values = unpack_u32x8(eight_values);
		for i in 0..(eight_values.len()) {
			assert_eq!(eight_values[i], value);
		}
	}

	#[test]
	fn partial_works() {
		let mut base = Additive8x::splat(Additive(42));
		let updated = base.override_partial_continuous(3, &[Additive(7), Additive(6)]);
		let unpacked = base.unpack();
		assert_eq!(unpacked[0], Additive(42));
		assert_eq!(unpacked[1], Additive(42));
		assert_eq!(unpacked[2], Additive(42));
		assert_eq!(unpacked[3], Additive(7));
		assert_eq!(unpacked[4], Additive(6));
		assert_eq!(unpacked[5], Additive(42));
		assert_eq!(unpacked[6], Additive(42));
		assert_eq!(unpacked[7], Additive(42));
	}
}
