use derive_more::{BitXor,BitXorAssign,Add,AddAssign,Sub,SubAssign};


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


include!("f2e16.tables.rs");


/// Compute tables determined solely by the field, which never depend
/// upon the FFT domain or erasure coding paramaters.
///
/// We compute `LOG_TABLE` and `EXP_TABLE` here of course.  We compute
/// the Walsh transform table `LOG_WALSH` here too because we never figured
/// out how to shrink `LOG_WALSH` below the size of the full field (TODO).
/// We thus assume it depends only upon the field for now.
#[cfg(not(table_build_done))]
fn write_field_tables<W: std::io::Write>(mut w: W) -> std::io::Result<()> {
    let mut log_table: [GFSymbol; FIELD_SIZE] = [0_u16; FIELD_SIZE];
    let mut exp_table: [GFSymbol; FIELD_SIZE] = [0_u16; FIELD_SIZE];

	let mas: Elt = (1 << FIELD_BITS - 1) - 1;
	let mut state: usize = 1;
	for i in 0_usize..(ONEMASK as usize) {
		exp_table[state] = i as Elt;
		if (state >> FIELD_BITS - 1) != 0 {
			state &= mas as usize;
			state = state << 1_usize ^ GENERATOR as usize;
		} else {
			state <<= 1;
		}
	}
	exp_table[0] = ONEMASK;

	log_table[0] = 0;
	for i in 0..FIELD_BITS {
		for j in 0..(1 << i) {
			log_table[j + (1 << i)] = log_table[j] ^ BASE[i];
		}
	}
	for i in 0..FIELD_SIZE {
		log_table[i] = exp_table[log_table[i] as usize];
	}

	for i in 0..FIELD_SIZE {
		exp_table[log_table[i] as usize] = i as GFSymbol;
	}
	exp_table[ONEMASK as usize] = exp_table[0];

    write_const(&mut w,"LOG_TABLE",&log_table,Some("[u16; FIELD_SIZE]")) ?;
    write_const(&mut w,"EXP_TABLE",&exp_table,Some("[u16; FIELD_SIZE]")) ?;

	// mem_cpy(&mut log_walsh[..], &log_table[..]);
    let log_walsh = log_table.clone();
    let mut log_walsh = unsafe { core::mem::transmute::<_,[Multiplier; FIELD_SIZE]>(log_walsh) };
	log_walsh[0] = Multiplier(0);
	walsh(&mut log_walsh[..], FIELD_SIZE);

    write_const(w,"LOG_WALSH",&log_walsh,Some("[Multiplier; FIELD_SIZE]"))
}


