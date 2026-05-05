use core::fmt;

use aether_core::{math::ComplexField, real::real::{Real, RealCast}};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FftAlgorithm {
    Auto,
    Radix2,
    Dft,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FftDirection {
    Forward,
    Inverse,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FftError {
    EmptyInput,
    NonPowerOfTwoLength,
    LengthMismatch,
    ScratchLengthMismatch,
}

impl fmt::Display for FftError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "fft input must not be empty"),
            Self::NonPowerOfTwoLength => write!(f, "fft input length must be a power of two"),
            Self::LengthMismatch => write!(f, "input and output lengths must match"),
            Self::ScratchLengthMismatch => write!(f, "scratch length must match data length"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for FftError {}

pub fn fft_in_place<C: ComplexField, const N: usize>(data: &mut [C; N]) -> Result<(), FftError> {
    transform_in_place_slice(data.as_mut_slice(), FftDirection::Forward)
}

pub fn ifft_in_place<C: ComplexField, const N: usize>(data: &mut [C; N]) -> Result<(), FftError> {
    transform_in_place_slice(data.as_mut_slice(), FftDirection::Inverse)
}

#[cfg(feature = "std")]
pub(crate) fn fft_in_place_slice<C: ComplexField>(data: &mut [C]) -> Result<(), FftError> {
    transform_in_place_slice(data, FftDirection::Forward)
}

pub fn transform_in_place_with_scratch<C: ComplexField, const N: usize>(
    data: &mut [C; N],
    scratch: &mut [C; N],
    direction: FftDirection,
    algorithm: FftAlgorithm,
) -> Result<(), FftError> {
    transform_in_place_with_scratch_slice(
        data.as_mut_slice(),
        scratch.as_mut_slice(),
        direction,
        algorithm,
    )
}

pub(crate) fn transform_in_place_with_scratch_slice<C: ComplexField>(
    data: &mut [C],
    scratch: &mut [C],
    direction: FftDirection,
    algorithm: FftAlgorithm,
) -> Result<(), FftError> {
    let len = data.len();
    if len == 0 {
        return Err(FftError::EmptyInput);
    }
    if scratch.len() != len {
        return Err(FftError::ScratchLengthMismatch);
    }

    match algorithm {
        FftAlgorithm::Auto => {
            if len.is_power_of_two() {
                radix2_transform_in_place(data, direction)
            } else {
                dft_into_slice(data, scratch, direction)?;
                data.copy_from_slice(scratch);
                Ok(())
            }
        }
        FftAlgorithm::Radix2 => radix2_transform_in_place(data, direction),
        FftAlgorithm::Dft => {
            dft_into_slice(data, scratch, direction)?;
            data.copy_from_slice(scratch);
            Ok(())
        }
    }
}

pub fn dft_into<C: ComplexField, const N: usize>(
    input: &[C; N],
    output: &mut [C; N],
    direction: FftDirection,
) -> Result<(), FftError> {
    dft_into_slice(input.as_slice(), output.as_mut_slice(), direction)
}

pub(crate) fn dft_into_slice<C: ComplexField>(
    input: &[C],
    output: &mut [C],
    direction: FftDirection,
) -> Result<(), FftError> {
    let len = input.len();
    if len == 0 {
        return Err(FftError::EmptyInput);
    }
    if output.len() != len {
        return Err(FftError::LengthMismatch);
    }

    let sign = match direction {
        FftDirection::Forward => -C::RealPart::ONE,
        FftDirection::Inverse => C::RealPart::ONE,
    };
    let len_real = C::RealPart::from_usize(len);

    let mut k = 0usize;
    while k < len {
        let mut sum = C::zero();
        let mut n = 0usize;
        while n < len {
            let angle = sign
                * C::RealPart::from_f64(2.0)
                * C::RealPart::PI
                * C::RealPart::from_usize(k * n)
                / len_real;
            let twiddle = C::from_polar(C::RealPart::ONE, angle);
            sum += input[n] * twiddle;
            n += 1;
        }

        if direction == FftDirection::Inverse {
            sum /= len_real;
        }
        output[k] = sum;
        k += 1;
    }

    Ok(())
}

fn transform_in_place_slice<C: ComplexField>(
    data: &mut [C],
    direction: FftDirection,
) -> Result<(), FftError> {
    radix2_transform_in_place(data, direction)
}

fn radix2_transform_in_place<C: ComplexField>(
    data: &mut [C],
    direction: FftDirection,
) -> Result<(), FftError> {
    let len = data.len();
    if len == 0 {
        return Err(FftError::EmptyInput);
    }
    if !len.is_power_of_two() {
        return Err(FftError::NonPowerOfTwoLength);
    }

    bit_reverse_permute(data);

    let sign = match direction {
        FftDirection::Forward => -C::RealPart::ONE,
        FftDirection::Inverse => C::RealPart::ONE,
    };

    let mut block_len = 2usize;
    while block_len <= len {
        let half_block = block_len / 2;
        let theta = sign * C::RealPart::from_f64(2.0) * C::RealPart::PI / C::RealPart::from_usize(block_len);
        let w_len = C::from_polar(C::RealPart::ONE, theta);

        let mut block_start = 0usize;
        while block_start < len {
            let mut w = C::from_real(C::RealPart::ONE);
            let mut offset = 0usize;
            while offset < half_block {
                let even_index = block_start + offset;
                let odd_index = even_index + half_block;

                let even = data[even_index];
                let odd = data[odd_index] * w;

                data[even_index] = even + odd;
                data[odd_index] = even - odd;

                w *= w_len;
                offset += 1;
            }

            block_start += block_len;
        }

        block_len *= 2;
    }

    if direction == FftDirection::Inverse {
        let scale = C::RealPart::from_usize(len);
        for value in data.iter_mut() {
            *value /= scale;
        }
    }

    Ok(())
}

fn bit_reverse_permute<C: ComplexField>(data: &mut [C]) {
    let len = data.len();
    let shift = usize::BITS - len.trailing_zeros();

    let mut index = 0usize;
    while index < len {
        let reversed = index.reverse_bits() >> shift;
        if reversed > index {
            data.swap(index, reversed);
        }
        index += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use aether_core::math::Complex;

    fn assert_close(actual: Complex<f64>, expected: Complex<f64>, tolerance: f64) {
        assert!(
            (actual.real() - expected.real()).abs() <= tolerance,
            "real mismatch: {:?} vs {:?}",
            actual,
            expected
        );
        assert!(
            (actual.imag() - expected.imag()).abs() <= tolerance,
            "imag mismatch: {:?} vs {:?}",
            actual,
            expected
        );
    }

    fn coherent_real_tone<const N: usize>(bin: usize, amplitude: f64) -> [Complex<f64>; N] {
        core::array::from_fn(|index| {
                let angle = 2.0 * core::f64::consts::PI * bin as f64 * index as f64 / N as f64;
                Complex::from_real(amplitude * angle.cos())
            })
    }

    fn magnitude_spectrum<const N: usize>(values: &[Complex<f64>; N]) -> [f64; N] {
        core::array::from_fn(|index| values[index].norm())
    }

    #[test]
    fn fft_places_a_coherent_tone_in_the_expected_bins() {
        let len = 16usize;
        let tone_bin = 3usize;
        let mut data = coherent_real_tone::<16>(tone_bin, 1.0);

        fft_in_place(&mut data).unwrap();

        let magnitudes = magnitude_spectrum(&data);
        let expected_magnitude = len as f64 / 2.0;

        for (index, magnitude) in magnitudes.into_iter().enumerate() {
            if index == tone_bin || index == len - tone_bin {
                assert!(
                    (magnitude - expected_magnitude).abs() <= 1e-9,
                    "expected dominant tone at bin {index}, got {magnitude}"
                );
            } else {
                assert!(magnitude <= 1e-9, "unexpected leakage at bin {index}: {magnitude}");
            }
        }
    }

    #[test]
    fn fft_then_ifft_round_trips_a_signal_like_input() {
        let original: [Complex<f64>; 16] = core::array::from_fn(|index| {
            coherent_real_tone::<16>(3, 1.0)[index] + Complex::from_real(0.25 * index as f64 / 16.0)
        });

        let mut transformed = original;
        fft_in_place(&mut transformed).unwrap();
        ifft_in_place(&mut transformed).unwrap();

        for (actual, expected) in transformed.into_iter().zip(original) {
            assert_close(actual, expected, 1e-9);
        }
    }

    #[test]
    fn rejects_non_power_of_two_length() {
        let mut data = [Complex::<f64>::zero(); 3];
        let result = fft_in_place(&mut data);
        assert_eq!(result, Err(FftError::NonPowerOfTwoLength));
    }

    #[test]
    fn dft_matches_fft_for_power_of_two_tone() {
        let input = coherent_real_tone::<8>(2, 1.0);

        let mut fft_output = input;
        fft_in_place(&mut fft_output).unwrap();

        let mut dft_output = [Complex::zero(); 8];
        dft_into(&input, &mut dft_output, FftDirection::Forward).unwrap();

        for (actual, expected) in fft_output.into_iter().zip(dft_output) {
            assert_close(actual, expected, 1e-9);
        }
    }

    #[test]
    fn auto_falls_back_to_dft_for_non_power_of_two_tone() {
        let input = coherent_real_tone::<6>(2, 1.0);
        let mut data = input;
        let mut scratch = [Complex::zero(); 6];

        transform_in_place_with_scratch(
            &mut data,
            &mut scratch,
            FftDirection::Forward,
            FftAlgorithm::Auto,
        )
        .unwrap();

        let mut expected = [Complex::zero(); 6];
        dft_into(&input, &mut expected, FftDirection::Forward).unwrap();

        for (actual, expected) in data.into_iter().zip(expected) {
            assert_close(actual, expected, 1e-9);
        }
    }

    #[test]
    fn dft_inverse_round_trips_a_short_signal() {
        let base = coherent_real_tone::<5>(1, 1.0);
        let input = core::array::from_fn(|index| base[index] + Complex::from_real(0.1 * index as f64));

        let mut forward = [Complex::zero(); 5];
        dft_into(&input, &mut forward, FftDirection::Forward).unwrap();

        let mut inverse = [Complex::zero(); 5];
        dft_into(&forward, &mut inverse, FftDirection::Inverse).unwrap();

        for (actual, expected) in inverse.into_iter().zip(input) {
            assert_close(actual, expected, 1e-9);
        }
    }

    #[test]
    fn public_const_api_works_for_fixed_window_sizes() {
        let input = coherent_real_tone::<5>(2, 1.0);
        let mut data = input;
        let mut scratch = [Complex::zero(); 5];

        transform_in_place_with_scratch(&mut data, &mut scratch, FftDirection::Forward, FftAlgorithm::Auto)
            .unwrap();

        let mut expected = [Complex::zero(); 5];
        dft_into(&input, &mut expected, FftDirection::Forward).unwrap();

        for (actual, expected) in data.into_iter().zip(expected) {
            assert_close(actual, expected, 1e-9);
        }
    }

    #[test]
    fn const_fft_api_places_a_coherent_tone_in_the_expected_bins() {
        let mut data = [
            Complex::from_real(1.0),
            Complex::from_real(0.0),
            Complex::from_real(-1.0),
            Complex::from_real(0.0),
            Complex::from_real(1.0),
            Complex::from_real(0.0),
            Complex::from_real(-1.0),
            Complex::from_real(0.0),
        ];

        fft_in_place(&mut data).unwrap();

        let magnitudes = magnitude_spectrum(&data);
        for (index, magnitude) in magnitudes.into_iter().enumerate() {
            if index == 2 || index == 6 {
                assert!((magnitude - 4.0).abs() <= 1e-9);
            } else {
                assert!(magnitude <= 1e-9, "unexpected leakage at bin {index}: {magnitude}");
            }
        }
    }

    #[test]
    fn const_auto_api_uses_stack_scratch_for_non_power_of_two() {
        let input = [
            Complex::from_real(1.0),
            Complex::from_real(-0.5),
            Complex::from_real(-0.5),
            Complex::from_real(1.0),
            Complex::from_real(-0.5),
        ];
        let mut data = input;
        let mut scratch = [Complex::zero(); 5];

        transform_in_place_with_scratch(
            &mut data,
            &mut scratch,
            FftDirection::Forward,
            FftAlgorithm::Auto,
        )
        .unwrap();

        let mut expected = [Complex::zero(); 5];
        dft_into(&input, &mut expected, FftDirection::Forward).unwrap();

        for (actual, expected) in data.into_iter().zip(expected) {
            assert_close(actual, expected, 1e-9);
        }
    }
}