#![cfg_attr(all(feature = "no_std", not(feature = "std")), no_std)]

pub use aether_core::real;
pub use aether_core::{attitude, coordinate, math, numerical_methods, reference_frame, utils};
pub mod celestial;
pub mod lunar;
pub mod models;
pub mod terrestrial;
