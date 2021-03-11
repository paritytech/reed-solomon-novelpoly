use derive_more::{BitXor,BitXorAssign,Add,AddAssign,Sub,SubAssign};

#[cfg(not(table_bootstrap_complete))]
pub(crate) const LOG_TABLE: [u16; FIELD_SIZE] = [0; FIELD_SIZE];
#[cfg(not(table_bootstrap_complete))]
pub(crate) const LOG_WALSH: [u16; FIELD_SIZE] = [0; FIELD_SIZE];
#[cfg(not(table_bootstrap_complete))]
pub(crate) const EXP_TABLE: [u16; FIELD_SIZE] = [0; FIELD_SIZE];

// must be placed in a separate file, such that the preproc never tries to eval OUT_DIR
// in env which does not exist in the build.rs case
#[cfg(table_bootstrap_complete)]
include!(concat!(env!("OUT_DIR"), "/table_f2e16.rs"));


pub type Elt = u16;
pub type Wide = u32;

pub const FIELD_BITS: usize = 16;
pub const FIELD_SIZE: usize = 1_usize << FIELD_BITS;

#[derive(Clone,Copy,Debug,Default,BitXor,BitXorAssign,PartialEq,Eq)] // PartialOrd,Ord
pub struct Additive(pub Elt);
impl Additive {
    pub fn to_wide(self) -> Wide { self.0 as Wide }
    pub fn from_wide(x: Wide) -> Additive { Additive(x as Elt) }

    pub const ZERO: Additive = Additive(0u16);
    // pub const ONE: Additive = Additive(???);

    /// Return multiplier prepared form
    pub fn to_multiplier(self) -> Multiplier {
        Multiplier( LOG_TABLE[self.0 as usize] )
    }

    /// Return a*EXP_TABLE[b] over GF(2^r)
    pub fn mul(self, other: Multiplier) -> Additive {
    	if self == Self::ZERO { return Self::ZERO; }
    	let log = (LOG_TABLE[self.0 as usize] as u32) + other.0 as u32;
    	let offset = (log & ONEMASK as u32) + (log >> FIELD_BITS);
    	Additive(EXP_TABLE[offset as usize])
    }

    /// Multiply field elements by a single multiplier, using SIMD if available
    pub fn mul_assign_slice(selfy: &mut [Self], other: Multiplier) {
        // TODO: SIMD
        for s in selfy { *s = s.mul(other); }
    }
}

#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,PartialEq,Eq)] // Default, PartialOrd,Ord
pub struct Multiplier(pub u16);
impl Multiplier {
    pub fn to_wide(self) -> u32 { self.0 as u32 }
    pub fn from_wide(x: u32) -> Multiplier { Multiplier(x as u16) }
}


/// Fast Walsh–Hadamard transform over modulo ONEMASK
pub fn walsh(data: &mut [Multiplier], size: usize) {
	let mut depart_no = 1_usize;
	while depart_no < size {
		let mut j = 0;
		let depart_no_next = depart_no << 1;
		while j < size {
			for i in j..(depart_no + j) {
				// We deal with data in log form here, but field form looks like:
				//			 data[i] := data[i] / data[i+depart_no]
				// data[i+depart_no] := data[i] * data[i+depart_no]
				let mask = ONEMASK as Wide;
				let tmp2: Wide = data[i].to_wide() + mask - data[i + depart_no].to_wide();
				let tmp1: Wide = data[i].to_wide() + data[i + depart_no].to_wide();
				data[i] = Multiplier((
                    (tmp1 & mask) + (tmp1 >> FIELD_BITS)
                ) as Elt);
				data[i + depart_no] = Multiplier((
                    (tmp2 & mask) + (tmp2 >> FIELD_BITS)
                ) as Elt);
			}
			j += depart_no_next;
		}
		depart_no = depart_no_next;
	}
}


/* Needs Cleanup  */

pub type GFSymbol = Elt;
pub const ONEMASK: GFSymbol = (FIELD_SIZE - 1) as GFSymbol;

/// Quotient ideal generator given by tail of irreducible polynomial
pub const GENERATOR: GFSymbol = 0x2D; // x^16 + x^5 + x^3 + x^2 + 1

// Cantor basis
pub const BASE: [GFSymbol; FIELD_BITS] =
	[1_u16, 44234, 15374, 5694, 50562, 60718, 37196, 16402, 27800, 4312, 27250, 47360, 64952, 64308, 65336, 39198];
