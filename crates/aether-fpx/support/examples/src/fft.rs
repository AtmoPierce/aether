use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use aether_core::math::{Complex, ComplexField};
use aether_core::real::Real;
use aether_fft::{fft_in_place, ifft_in_place};
use aether_viz::{plot_series_with_config, PlotConfig, PlotStyle, XYSeries};
use aether_fpx::ocp::fp4::Fp4E2M1;
use aether_fpx::ocp::fp6::Fp6E3M2;
use aether_fpx::ocp::fp8::{Fp8E4M3, Fp8E5M2};

const N: usize = 4;
const PRECISION_COUNT: usize = 8;
const INPUT: [f64; N] = [0.5, -0.5, 1.0, -1.0];
const REFERENCE_SPECTRUM: [[f64; 2]; N] = [[0.0, 0.0], [-0.5, -0.5], [3.0, 0.0], [-0.5, 0.5]];
const OUTPUT_DIR: &str = "outputs/fft";
const CSV_PATH: &str = "outputs/fft/fft.csv";
const SPECTRUM_SVG_PATH: &str = "outputs/fft/spectrum.svg";

const LABELS: [&str; PRECISION_COUNT] = [
    "FP4 E2M1", "FP6 E3M2", "FP8 E4M3", "FP8 E5M2", "f16", "f32", "f64", "f128",
];

const CSV_NAMES: [&str; PRECISION_COUNT] = [
    "fpx_fp4_e2m1",
    "fpx_fp6_e3m2",
    "fpx_fp8_e4m3",
    "fpx_fp8_e5m2",
    "fpx_f16",
    "fpx_f32",
    "fpx_f64",
    "fpx_f128",
];

#[derive(Clone, Copy, Debug)]
struct FftResult {
    spectrum: [[f64; 2]; N],
    magnitudes: [f64; N],
    relative_spectrum_error: f64,
    round_trip_rmse: f64,
    parseval_relative_error: f64,
}

fn magnitude(value: [f64; 2]) -> f64 {
    f64::hypot(value[0], value[1])
}

fn evaluate<T: Real + Copy>() -> FftResult {
    let original: [Complex<T>; N] =
        std::array::from_fn(|index| Complex::from_real(T::from_f64(INPUT[index])));
    let mut transformed = original;
    fft_in_place(&mut transformed).expect("four-point input is a supported radix-2 length");

    let spectrum = std::array::from_fn(|index| {
        [
            transformed[index].real().to_f64(),
            transformed[index].imag().to_f64(),
        ]
    });
    let magnitudes = spectrum.map(magnitude);
    let reference_energy = REFERENCE_SPECTRUM
        .iter()
        .map(|value| value[0] * value[0] + value[1] * value[1])
        .sum::<f64>();
    let error_energy = spectrum
        .iter()
        .zip(REFERENCE_SPECTRUM)
        .map(|(actual, expected)| {
            let real_error = actual[0] - expected[0];
            let imag_error = actual[1] - expected[1];
            real_error * real_error + imag_error * imag_error
        })
        .sum::<f64>();

    let time_energy = INPUT.iter().map(|value| value * value).sum::<f64>();
    let spectral_energy = spectrum
        .iter()
        .map(|value| value[0] * value[0] + value[1] * value[1])
        .sum::<f64>()
        / N as f64;

    let mut reconstructed = transformed;
    ifft_in_place(&mut reconstructed).expect("four-point input is a supported radix-2 length");
    let round_trip_rmse = f64::sqrt(
        reconstructed
            .iter()
            .zip(INPUT)
            .map(|(actual, expected)| {
                let real_error = actual.real().to_f64() - expected;
                let imag_error = actual.imag().to_f64();
                real_error * real_error + imag_error * imag_error
            })
            .sum::<f64>()
            / N as f64,
    );

    FftResult {
        spectrum,
        magnitudes,
        relative_spectrum_error: f64::sqrt(error_energy / reference_energy),
        round_trip_rmse,
        parseval_relative_error: (spectral_energy - time_energy).abs() / time_energy,
    }
}

fn evaluate_all() -> [FftResult; PRECISION_COUNT] {
    [
        evaluate::<Fp4E2M1>(),
        evaluate::<Fp6E3M2>(),
        evaluate::<Fp8E4M3>(),
        evaluate::<Fp8E5M2>(),
        evaluate::<f16>(),
        evaluate::<f32>(),
        evaluate::<f64>(),
        evaluate::<f128>(),
    ]
}

fn write_csv(results: &[FftResult; PRECISION_COUNT], path: &Path) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    writeln!(
        output,
        "format,bin,real,imag,magnitude,relative_spectrum_error,round_trip_rmse,parseval_relative_error"
    )?;
    for (precision, result) in results.iter().enumerate() {
        for bin in 0..N {
            writeln!(
                output,
                "{},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                CSV_NAMES[precision],
                bin,
                result.spectrum[bin][0],
                result.spectrum[bin][1],
                result.magnitudes[bin],
                result.relative_spectrum_error,
                result.round_trip_rmse,
                result.parseval_relative_error,
            )?;
        }
    }
    output.flush()
}

fn draw_spectrum_plot(
    results: &[FftResult; PRECISION_COUNT],
) -> Result<(), Box<dyn std::error::Error>> {
    let bins = (0..N).map(|bin| bin as f64).collect::<Vec<_>>();
    let reference_magnitudes = REFERENCE_SPECTRUM.map(magnitude);
    let mut series = Vec::with_capacity(PRECISION_COUNT + 1);
    series.push(
        XYSeries::new(&bins, &reference_magnitudes)
            .with_label("analytic DFT reference")
            .with_style(PlotStyle::DashedLineWithPalette(8)),
    );
    for (precision, result) in results.iter().enumerate() {
        series.push(
            XYSeries::new(&bins, &result.magnitudes)
                .with_label(LABELS[precision])
                .with_style(PlotStyle::LineWithPalette(precision)),
        );
    }
    plot_series_with_config(
        &series,
        &PlotConfig::publication("Four-point FFT magnitude spectrum")
            .with_font_scale(1.6)
            .with_x_label("frequency bin")
            .with_y_label("magnitude")
            .with_save_path(SPECTRUM_SVG_PATH),
    )?;
    Ok(())
}

pub(crate) fn run() -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(OUTPUT_DIR)?;
    let results = evaluate_all();
    write_csv(&results, Path::new(CSV_PATH))?;
    draw_spectrum_plot(&results)?;

    println!("wrote {CSV_PATH}");
    println!("wrote {SPECTRUM_SVG_PATH}");
    println!("input = {:?}", INPUT);
    for (label, result) in LABELS.iter().zip(results) {
        println!(
            "{label:>9}: relative spectrum error = {:.9e}, round-trip RMSE = {:.9e}, Parseval error = {:.9e}",
            result.relative_spectrum_error,
            result.round_trip_rmse,
            result.parseval_relative_error,
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_precision_produces_finite_fft_metrics() {
        for (label, result) in LABELS.iter().zip(evaluate_all()) {
            assert!(
                result.spectrum.into_iter().flatten().all(f64::is_finite)
                    && result.magnitudes.into_iter().all(f64::is_finite)
                    && result.relative_spectrum_error.is_finite()
                    && result.round_trip_rmse.is_finite()
                    && result.parseval_relative_error.is_finite(),
                "{label} produced a non-finite FFT result",
            );
        }
    }
}
