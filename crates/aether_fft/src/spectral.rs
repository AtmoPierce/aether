use std::vec::Vec;

use aether_core::{math::Complex, real::Real};

use crate::{fft::fft_in_place_slice, FftError};

#[derive(Debug, Clone, PartialEq)]
pub struct PsdEstimate<T: Real> {
    pub frequencies_hz: Vec<T>,
    pub power_density: Vec<T>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpectralError {
    NonPositiveSampleRate,
    Fft(FftError),
}

impl core::fmt::Display for SpectralError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonPositiveSampleRate => write!(f, "sample rate must be positive"),
            Self::Fft(err) => write!(f, "{err}"),
        }
    }
}

impl std::error::Error for SpectralError {}

impl From<FftError> for SpectralError {
    fn from(value: FftError) -> Self {
        Self::Fft(value)
    }
}

pub fn real_periodogram<T: Real>(samples: &[T], sample_rate_hz: T) -> Result<PsdEstimate<T>, SpectralError> {
    if sample_rate_hz <= T::ZERO {
        return Err(SpectralError::NonPositiveSampleRate);
    }

    let len = samples.len();
    let mut spectrum = samples
        .iter()
        .copied()
        .map(Complex::from_real)
        .collect::<Vec<_>>();
    fft_in_place_slice(&mut spectrum)?;

    let n = T::from_usize(len);
    let bin_count = len / 2 + 1;
    let mut frequencies_hz = Vec::with_capacity(bin_count);
    let mut power_density = Vec::with_capacity(bin_count);

    let scale = sample_rate_hz * n;
    for (index, value) in spectrum.iter().take(bin_count).copied().enumerate() {
        let mut power = value.norm_sqr() / scale;
        let is_dc = index == 0;
        let is_nyquist = len % 2 == 0 && index == len / 2;
        if !is_dc && !is_nyquist {
            power = power * T::from_f64(2.0);
        }

        frequencies_hz.push(T::from_usize(index) * sample_rate_hz / n);
        power_density.push(power);
    }

    Ok(PsdEstimate {
        frequencies_hz,
        power_density,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn periodogram_of_zero_signal_is_zero() {
        let samples = vec![0.0_f64; 8];
        let psd = real_periodogram(&samples, 8.0).unwrap();
        assert_eq!(psd.frequencies_hz.len(), 5);
        assert!(psd.power_density.iter().all(|value| value.abs() <= 1e-12));
    }

    #[test]
    fn rejects_non_positive_sample_rate() {
        let samples = vec![0.0_f64; 8];
        let result = real_periodogram(&samples, 0.0);
        assert_eq!(result, Err(SpectralError::NonPositiveSampleRate));
    }
}