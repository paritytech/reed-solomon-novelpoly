#![cfg(all(target_feature = "avx", feature = "avx"))]

#[cfg(feature = "f256")]
pub mod f256;
#[cfg(feature = "f256")]
pub use self::f256::*;

pub mod f2e16;
pub use self::f2e16::*;
