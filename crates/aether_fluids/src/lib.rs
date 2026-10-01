#![cfg_attr(all(feature = "no_std", not(feature = "std")), no_std)]

pub mod tanks;
pub use tanks::*;
