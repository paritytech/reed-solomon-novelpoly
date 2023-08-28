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
		write!(f, "{:04x}", self.0)
	}
}

impl std::fmt::Debug for Additive {
	fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		write!(f, "{:04x}", self.0)
	}
}

use core::ops::{BitXor, BitXorAssign};
use core::arch::x86_64::*;

#[allow(non_camel_case_types)]
pub type u16x8 = __m128i;

#[allow(non_camel_case_types)]
pub type u32x8 = __m256i;

#[repr(transparent)]
#[derive(Clone, Copy)]
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

impl std::fmt::Debug for Additive8x {
	fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		let unpacked = self.unpack();
		write!(
			f,
			"{}, {}, {}, {}, {}, {}, {}, {}",
			unpacked[0], unpacked[1], unpacked[2], unpacked[3], unpacked[4], unpacked[5], unpacked[6], unpacked[7]
		)
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
pub(crate) fn splat_u16x8(v: u16) -> u16x8 {
	unsafe { _mm_set1_epi16(v as i16) }
}

#[inline(always)]
pub(crate) fn splat_u32x8(v: u32) -> u32x8 {
	unsafe { _mm256_set1_epi32(v as i32) }
}

#[allow(dead_code)]
#[inline(always)]
pub(crate) fn unpack_u32x8(v: u32x8) -> [u32; 8] {
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
pub(crate) fn unpack_u16x8(v: u16x8) -> [u16; 8] {
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
pub(crate) fn load_u16x8(src: &[u16], idx: [u16; 8]) -> u16x8 {
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

#[inline(always)]
pub(crate) fn clipping_cast(data: u32x8) -> u16x8 {
	unsafe {
		let mask = splat_u32x8(0x0000_FFFF);
		// need
		let data = _mm256_and_si256(data, mask);
		// data_lo zero_lo data_hi zero_hi
		let packed = _mm256_packus_epi32(data, data);

		// use `a` as above in the lower and higher 128 bits
		const PERMUTE: i32 = 0b11_01_10_00;
		let data = _mm256_permute4x64_epi64(packed, PERMUTE);
		let data = _mm256_castsi256_si128(data);
		data
	}
}

#[inline(always)]
pub(crate) fn expand_cast(data: u16x8) -> u32x8 {
	unsafe { _mm256_cvtepu16_epi32(data) }
}

impl Additive8x {
	pub const LANE: usize = 8;

	#[inline(always)]
	pub fn zero() -> Self {
		Self(unsafe { _mm_setzero_si128() })
	}

	#[cfg(table_bootstrap_complete)]
	pub fn mul(&self, other: Multiplier) -> Self {
		if *self == Self::zero() {
			return Self::zero();
		}
		// let log = (LOG_TABLE[self.0 as usize] as Wide) + other.0 as Wide;
		// let offset = (log & ONEMASK as Wide) + (log >> FIELD_BITS);
		// Additive(EXP_TABLE[offset as usize])

		unsafe {
			// spread the multiplier
			let other = splat_u32x8(other.0 as Wide);
			// unpack_u32x8(other);

			// get the indices
			let idx = unpack_u16x8(self.0);
			// load from the table based on the indices
			let logtable = load_u16x8(LOG_TABLE.as_slice(), idx);
			// dbg!(unpack_u16x8(logtable));

			let logtable = _mm256_cvtepu16_epi32(logtable);

			let log = _mm256_add_epi32(logtable, other);
			// dbg!(unpack_u32x8(log));

			// (log & ONEMASK) + (log >> shift)
			let onemasks = splat_u32x8(ONEMASK as Wide);
			let offset_sum_left = _mm256_and_si256(log, onemasks);
			// dbg!(unpack_u32x8(offset_sum_left));
			let offset_sum_right = _mm256_srli_epi32(log, FIELD_BITS as i32);
			// dbg!(unpack_u32x8(offset_sum_right));
			let offset = _mm256_add_epi32(offset_sum_left, offset_sum_right);
			// dbg!(unpack_u32x8(offset));

			let offset = clipping_cast(offset);
			let offset = unpack_u16x8(offset);
			// dbg!(offset);
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

	fn test_mul(additive: Additive, mpy: Multiplier) {
		assert!(cfg!(target_feature = "avx2"), "Tests are meaningless without avx2 target feature");

		println!("0x{additive:?} .mul( {mpy:?} ) ...",);
		let single = additive;
		let faster8 = Additive8x::from([single; 8]);

		let res_faster8 = faster8.mul(mpy);
		let res_faster8 = Into::<[Additive; Additive8x::LANE]>::into(res_faster8);
		let res = single.mul(mpy);

		for i in 0..Additive8x::LANE {
			assert_eq!(res_faster8[i], res, " @ [{i}]");
		}
	}

	#[test]
	fn identical_mul_regressions() {
		assert!(cfg!(target_feature = "avx2"), "Tests are meaningless without avx2 target feature");
		test_mul(Additive(0x0003), Multiplier(20182));

		test_mul(Additive(0xFA1C), Multiplier(63493));

		test_mul(Additive(0), Multiplier(1));

		test_mul(Additive(1), Multiplier(0));

		test_mul(Additive(0x16e7), Multiplier(18124));

		test_mul(Additive(0x3d3d), Multiplier(15677));
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
		base.override_partial_continuous(3, &[Additive(7), Additive(6)]);
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

	#[test]
	fn expand_and_clip_cast() {
		assert!(cfg!(target_feature = "avx2"), "Tests are meaningless without avx2 target feature");

		let v = dbg!(splat_u32x8(0xFF11_AABB));
		let x = dbg!(clipping_cast(v));
		assert_eq!(unpack_u16x8(splat_u16x8(0xAA_BB)), unpack_u16x8(x));

		let y = expand_cast(x);
		assert_eq!(unpack_u32x8(splat_u32x8(0x0000_AA_BB)), unpack_u32x8(y));
	}
}
