#![feature(f16, f128)]
#![cfg_attr(target_os = "none", no_std)]

pub mod trig;
pub use trig::*;
pub mod bench;
pub mod policy;
pub use policy::policy::*;
pub use policy::*;

#[cfg(feature = "simd")]
pub mod simd;

#[cfg(all(feature = "simd", not(target_os = "none")))]
pub mod bench_paths;

#[cfg(all(feature = "simd", not(target_os = "none")))]
pub mod public_simd;

#[cfg(feature = "ocp")]
pub mod ocp;
#[cfg(feature = "ocp")]
pub use ocp::*;
