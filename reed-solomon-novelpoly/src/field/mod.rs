use crate::errors::*;
use crate::util::*;

#[macro_use]
mod gen;

#[cfg(feature = "f256")]
pub mod f256;

pub mod f2e16;

#[cfg(all(target_feature = "avx", feature = "avx"))]
pub mod faster8;
