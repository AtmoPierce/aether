#![cfg_attr(all(feature = "no_std", not(feature = "std")), no_std)]

pub use aether_core::{attitude, coordinate, math, numerical_methods, reference_frame, utils};
pub use aether_core::real;
pub mod models;
pub mod celestial;
pub mod terrestrial;
pub mod lunar;