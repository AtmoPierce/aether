#![cfg_attr(feature = "f16", feature(f16))]
#![cfg_attr(feature = "f128", feature(f128))]
#![cfg_attr(all(feature = "no_std", not(feature = "std")), no_std)]

pub mod attitude;
pub mod coordinate;
pub mod math;
pub mod models;
pub mod numerical_methods;
pub mod real;
pub mod reference_frame;
pub mod utils;

#[cfg(target_arch = "aarch64")]
pub use math::arch::arm::neon as arch_neon;
