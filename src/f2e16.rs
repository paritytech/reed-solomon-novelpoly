use derive_more::{BitXor,BitXorAssign};


pub type Elt = u16;
pub type Wide = u32;

pub const FIELD_BITS: usize = 16;
pub const FIELD_SIZE: usize = 1_usize << FIELD_BITS;

#[derive(Clone,Copy,Debug,BitXor,BitXorAssign,PartialEq,Eq)] // PartialOrd,Ord
pub struct Additive(pub Elt);
impl Additive {
    pub fn to_wide(self) -> Wide { self.0 as Wide }
    pub fn from_wide(x: Wide) -> Additive { Additive(x as Elt) }

    pub const ZERO: Additive = Additive(0u16);
    // pub const ONE: Additive = Additive(???);

    /// Return multiplier prepared form
    pub fn to_multiplier(self) -> Multiplier {
        Multiplier( unsafe { LOG_TABLE[self.0 as usize] } )
    }

    /// Return a*EXP_TABLE[b] over GF(2^r)
    pub fn mul(self, other: Multiplier) -> Additive {
    	if self == Self::ZERO { return Self::ZERO; }
    	unsafe {
    		let log = (LOG_TABLE[self.0 as usize] as u32) + other.0 as u32;
    		let offset = (log & ONEMASK as u32) + (log >> FIELD_BITS);
			Additive(EXP_TABLE[offset as usize])
    	}
    }

    /// Multiply field elements by a single multiplier, using SIMD if available
    pub fn mul_assign_slice(selfy: &mut [Self], other: Multiplier) {
        // TODO: SIMD
        for s in selfy { *s = s.mul(other); }
    }
}

#[derive(Clone,Copy,Debug,BitXor,BitXorAssign,PartialEq,Eq)] // PartialOrd,Ord
pub struct Multiplier(pub u16);
impl Multiplier {
    pub fn to_wide(self) -> u32 { self.0 as u32 }
    pub fn from_wide(x: u32) -> Multiplier { Multiplier(x as u16) }
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

// Compute LOG_TABLE and EXP_TABLE
#[cfg(not(table_build_done))]

fn write_tables<W: std::io::Write>(mut w: W) -> Result<(),std::io::Error> {
    let mut log_table: [GFSymbol; FIELD_SIZE] = [0_u16; FIELD_SIZE];
    let mut exp_table: [GFSymbol; FIELD_SIZE] = [0_u16; FIELD_SIZE];

	let mas: GFSymbol = (1 << FIELD_BITS - 1) - 1;
	let mut state: usize = 1;
	for i in 0_usize..(ONEMASK as usize) {
		exp_table[state] = i as GFSymbol;
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

    write_const(&mut w,"LOG_TABLE",&log_table,None) ?;
    write_const(&mut w,"EXP_TABLE",&exp_table,None)
}

