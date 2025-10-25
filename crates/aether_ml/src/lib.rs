#![cfg_attr(not(feature = "std"), no_std)]
//! Machine learning models ported from CS640
pub mod regression;
pub use regression::*;
pub mod classification;
pub use classification::*;
