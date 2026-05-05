#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(test)]
extern crate std;

mod fft;
#[cfg(feature = "std")]
mod spectral;

pub use aether_core::math::{Complex, ComplexField};
pub use fft::{
	dft_into,
	fft_in_place,
	ifft_in_place,
	transform_in_place_with_scratch,
	FftAlgorithm,
	FftDirection,
	FftError,
};
#[cfg(feature = "std")]
pub use spectral::{real_periodogram, PsdEstimate, SpectralError};