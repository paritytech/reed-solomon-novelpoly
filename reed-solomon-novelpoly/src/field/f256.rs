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

#[cfg(table_bootstrap_complete)]
#[cfg(all(target_feature = "avx", feature = "avx"))]
pub use faster8::f256::*;
