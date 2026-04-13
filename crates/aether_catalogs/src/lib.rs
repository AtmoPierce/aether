#![cfg_attr(all(feature = "no_std", not(feature = "std")), no_std)]

#[cfg(feature = "std")]
pub mod star_catalogs;
#[cfg(feature = "std")]
pub use star_catalogs::*;

#[cfg(feature = "world-magnetic-model")]
pub mod wmm;
#[cfg(feature = "world-magnetic-model")]
pub use wmm::*;