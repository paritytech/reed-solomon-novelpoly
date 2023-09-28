use derive_more::{Add, AddAssign, BitXor, BitXorAssign, Sub, SubAssign};

/// Additive via XOR form of f2e16
#[repr(transparent)]
#[derive(Clone, Copy, BitXor, BitXorAssign, PartialEq, Eq)] // PartialOrd,Ord
pub struct Additive(pub Elt);

impl AsRef<Elt> for Additive {
	fn as_ref(&self) -> &Elt {
		&self.0
	}
}

impl Additive {
    #[inline(always)]
	pub fn to_wide(self) -> Wide {
		self.0 as Wide
	}
    #[inline(always)]
	pub fn from_wide(x: Wide) -> Additive {
		Additive(x as Elt)
	}

	pub const ZERO: Additive = Additive(0);
	
	pub fn zero() -> Self {
		Self(0)
	}
}

#[cfg(table_bootstrap_complete)]
impl Additive {
	/// Return multiplier prepared form
    #[inline(always)]
	pub fn to_multiplier(self) -> Multiplier {
		Multiplier(LOG_TABLE[self.0 as usize])
	}

	/// Return a*EXP_TABLE[b] over GF(2^r)
    #[inline(always)]
	pub fn mul(self, other: Multiplier) -> Additive {
		if self == Self::ZERO {
			return Self::ZERO;
		}
		let log = (LOG_TABLE[self.0 as usize] as Wide) + other.0 as Wide;
		let offset = (log & ONEMASK as Wide) + (log >> FIELD_BITS);
		Additive(EXP_TABLE[offset as usize])
	}

	/// Multiply field elements by a single multiplier, using SIMD if available
    #[inline(always)]
	pub fn mul_assign_slice(selfy: &mut [Self], other: Multiplier) {
		for s in selfy {
			*s = s.mul(other);
		}
	}
}


/// Multiplicaiton friendly LOG form of f2e16
#[derive(Clone, Debug, Copy, Add, AddAssign, Sub, SubAssign, PartialEq, Eq)] // Default, PartialOrd,Ord
pub struct Multiplier(pub Elt);

impl Multiplier {
    #[inline(always)]
	pub fn to_wide(self) -> Wide {
		self.0 as Wide
	}
    #[inline(always)]
	pub fn from_wide(x: Wide) -> Multiplier {
		Multiplier(x as Elt)
	}
}

impl std::fmt::Display for Multiplier {
	fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		write!(f, "_{:04x}", self.0)
	}
}


/// Fast Walsh–Hadamard transform over modulo `ONEMASK`
#[inline(always)]
pub fn walsh(data: &mut [Multiplier], size: usize) {
	#[cfg(all(target_feature = "avx", table_bootstrap_complete))]
	walsh_faster8(data, size);
	#[cfg(not(all(target_feature = "avx", table_bootstrap_complete)))]
	walsh_plain(data, size);
}

pub fn walsh_plain(data: &mut [Multiplier], size: usize) {
	let mask = ONEMASK as Wide;
	let mut depart_no = 1_usize;
	while depart_no < size {
		let mut j = 0;
		while j < size {
			for i in j..(j+depart_no) {
				// We deal with data in log form here, but field form looks like:
				//			 data[i] := data[i] / data[i+depart_no]
				// data[i+depart_no] := data[i] * data[i+depart_no]

				let tmp2: Wide = data[i].to_wide() + mask - data[i + depart_no].to_wide();
				let tmp1: Wide = data[i].to_wide() + data[i + depart_no].to_wide();
				data[i] = Multiplier(((tmp1 & mask) + (tmp1 >> FIELD_BITS)) as Elt);
				data[i + depart_no] = Multiplier(((tmp2 & mask) + (tmp2 >> FIELD_BITS)) as Elt);
			}
			j += depart_no << 1;
		}
		depart_no <<= 1;
	}
}


#[cfg(table_bootstrap_complete)]
#[cfg(target_feature = "avx")]
pub fn walsh_faster8(data: &mut [Multiplier], size: usize) {
	const LANE: usize = 8;
	
	let mut depart_no = 1_usize;
	
	while depart_no < LANE {
		let mut j = 0;
		while j < size {
			for i in j..(j+depart_no) {
				// We deal with data in log form here, but field form looks like:
				//			 data[i] := data[i] / data[i+depart_no]
				// data[i+depart_no] := data[i] * data[i+depart_no]
				let mask = ONEMASK as Wide;
				let tmp2: Wide = data[i].to_wide() + mask - data[i + depart_no].to_wide();
				let tmp1: Wide = data[i].to_wide() + data[i + depart_no].to_wide();
				data[i] = Multiplier(((tmp1 & mask) + (tmp1 >> FIELD_BITS)) as Elt);
				data[i + depart_no] = Multiplier(((tmp2 & mask) + (tmp2 >> FIELD_BITS)) as Elt);
			}
			j += depart_no << 1;
		}
		depart_no <<= 1;
	}
	
	use core::arch::x86_64::*;
	use crate::f2e16::{u16x8, splat_u32x8, clipping_cast, expand_cast};

	// let mask = ONEMASK as Wide;
	let mask = splat_u32x8(ONEMASK as _);
	
	while depart_no < size {
		let mut j = 0;
		while j < size {

			for i in (j..(j+depart_no)).step_by(LANE) {
				// We deal with data in log form here, but field form looks like:
				//			 data[i] := data[i] / data[i+depart_no]
				// data[i+depart_no] := data[i] * data[i+depart_no]
				unsafe {
					let data0 = _mm_loadu_si128(core::intrinsics::transmute::<_, *const u16x8>(data[i..][..LANE].as_ptr()));
					// dbg!(unpack_u16x8(data0));

					let data1 = _mm_loadu_si128(core::intrinsics::transmute::<_, *const u16x8>(data[(i+depart_no)..][..LANE].as_ptr()));
					// dbg!(unpack_u16x8(data1));
					
					let data0 = expand_cast(data0);
					// dbg!(unpack_u32x8(data0));
					
					let data1 = expand_cast(data1);
					// dbg!(unpack_u32x8(data1));

					// let tmp2: Wide = data[i].to_wide() + mask - data[i + depart_no].to_wide();
					let tmp2 = _mm256_sub_epi32(_mm256_add_epi32(data0, mask), data1);
					// dbg!(unpack_u32x8(tmp2));
					
					// let tmp1: Wide = data[i].to_wide() + data[i + depart_no].to_wide();
					let tmp1 = _mm256_add_epi32(data0, data1);
					// dbg!(unpack_u32x8(tmp1));
					
					{
						// data[i] = Multiplier(((tmp1 & mask) + (tmp1 >> FIELD_BITS)) as Elt);
						let i_a = _mm256_and_si256(tmp1, mask);
						let i_b = _mm256_srli_epi32(tmp1, FIELD_BITS as _);
						let res0 = _mm256_add_epi32(i_a,i_b);
						// dbg!(unpack_u32x8(res0));
						let res0 = clipping_cast(res0);
						// dbg!(unpack_u16x8(res0));
					
						// store to &mut data[i..][..8]
						_mm_storeu_si128(core::intrinsics::transmute(data[i..][..LANE].as_ptr()), res0);
					}
					
					{
						// data[i + depart_no] = Multiplier(((tmp2 & mask) + (tmp2 >> FIELD_BITS)) as Elt);
						let i_c = _mm256_and_si256(tmp2, mask);
						let i_d = _mm256_srli_epi32(tmp2, FIELD_BITS as _);
						let res1 = _mm256_add_epi32(i_c,i_d);
						// dbg!(unpack_u32x8(res1));
						let res1 = clipping_cast(res1);
						// dbg!(unpack_u16x8(res1));

						// store to &mut data[(i+depart_no)..][..8]
						_mm_storeu_si128(core::intrinsics::transmute(data[(i+depart_no)..][..LANE].as_ptr()), res1);
					}
				}
			}
			j += depart_no << 1;
		}
		depart_no <<= 1;
	}
}

#[allow(unused)]
fn bitpoly_mul16(a: Wide, b: Wide) -> Wide {
    let mut r: Wide =0;
    for i in 0..FIELD_BITS {
        if (b>>i) & 1 != 0 {
			r ^= a<<i;
		}
    }
    r
}

#[allow(unused)]
fn gf_mul_bitpoly_reduced(a: Elt, b: Elt) -> Elt {
    use core::convert::TryInto;
    let len = FIELD_BITS;
    let mut r: Wide = bitpoly_mul16(a as Wide,b as Wide);
    let red : Wide = (1 << FIELD_BITS) + (GENERATOR as Wide);
    for i in (len..=(len*2-1)).rev() {
        if r & (1<<i) != 0 {
			r ^= red<<(i-len);
		}
    }
    r.try_into().unwrap()
}

#[test]
fn cantor_basis() {
    for w in BASE.windows(2) {
        let b = w[1];
        let square = gf_mul_bitpoly_reduced(b,b);
        let a = w[0];
        // let eq = if a == (square ^ b) { "==" } else { "!=" };
        // println!("{:#b} {} {:#b}\n", a, eq, square ^ b);
        assert_eq!(a, square ^ b);
    }
}


#[cfg(table_bootstrap_complete)]
#[cfg(target_feature = "avx")]
#[test]
fn walsh_output_plain_eq_faster8() {
	use reed_solomon_tester::*;
	use rand::prelude::*;
	
	const SZ: usize = FIELD_SIZE/2;
	
	let mut rng = rand::rngs::SmallRng::from_seed(SMALL_RNG_SEED);
	let mut data = [Multiplier(0_u16); SZ];
	{
		let data = unsafe {
			core::mem::transmute::<_, &mut [u16;SZ]>(&mut data)
		};
		rng.fill(&mut data[..]);
	}
	let mut d0 = data;
	let mut d1 = data;

	walsh_plain(&mut d0, SZ);
	eprintln!("\n\n\n\n <<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>\n\n\n\n");
	walsh_faster8(&mut d1, SZ);

	assert_eq!(d0, d1);
}
