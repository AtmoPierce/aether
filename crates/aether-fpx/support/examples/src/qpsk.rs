use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use aether_core::real::Real;
use aether_viz::{plot_series_with_config, PlotConfig, PlotStyle, XYSeries};
use aether_fpx::ocp::fp4::Fp4E2M1;
use aether_fpx::ocp::fp6::Fp6E3M2;
use aether_fpx::ocp::fp8::{Fp8E4M3, Fp8E5M2};

const PRECISION_COUNT: usize = 8;
const SYMBOL_COUNT: usize = 32;
const CHANNEL_PHASE_RAD: f64 = 0.5;
const NOISE_AMPLITUDE: f64 = 0.04;

// Gray-coded QPSK symbols. Each two-bit value maps to a constellation quadrant.
const PACKET: [u8; SYMBOL_COUNT] = [
    0, 1, 3, 2, 0, 3, 1, 2, 2, 1, 0, 3, 1, 1, 2, 0, 3, 2, 3, 0, 0, 2, 1, 3, 2, 0, 1, 2, 3, 1, 0, 3,
];

const OUTPUT_DIR: &str = "outputs/qpsk";
const CSV_PATH: &str = "outputs/qpsk/qpsk.csv";
const CONSTELLATION_SVG_PATH: &str = "outputs/qpsk/constellation.svg";

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
struct QpskResult {
    recovered: [[f64; 2]; SYMBOL_COUNT],
    decoded: [u8; SYMBOL_COUNT],
    symbol_errors: usize,
    bit_errors: usize,
    evm: f64,
    minimum_decision_margin: f64,
}

fn symbol_coordinates(symbol: u8) -> [f64; 2] {
    match symbol {
        0b00 => [0.5, 0.5],
        0b01 => [-0.5, 0.5],
        0b11 => [-0.5, -0.5],
        0b10 => [0.5, -0.5],
        _ => unreachable!("a QPSK symbol contains two bits"),
    }
}

fn channel_sample(index: usize, symbol: u8) -> [f64; 2] {
    let [in_phase, quadrature] = symbol_coordinates(symbol);
    let cosine = f64::cos(CHANNEL_PHASE_RAD);
    let sine = f64::sin(CHANNEL_PHASE_RAD);
    let noise_i = NOISE_AMPLITUDE * f64::sin((index as f64 + 1.0) * 1.7);
    let noise_q = NOISE_AMPLITUDE * f64::cos((index as f64 + 1.0) * 1.3);
    [
        in_phase * cosine - quadrature * sine + noise_i,
        in_phase * sine + quadrature * cosine + noise_q,
    ]
}

fn decode(in_phase: f64, quadrature: f64) -> u8 {
    match (in_phase >= 0.0, quadrature >= 0.0) {
        (true, true) => 0b00,
        (false, true) => 0b01,
        (false, false) => 0b11,
        (true, false) => 0b10,
    }
}

fn evaluate<T: Real + Copy>() -> QpskResult {
    let phase = T::from_f64(CHANNEL_PHASE_RAD);
    let cosine = phase.cos();
    let sine = phase.sin();
    let mut recovered = [[0.0; 2]; SYMBOL_COUNT];
    let mut decoded = [0; SYMBOL_COUNT];

    for (index, symbol) in PACKET.into_iter().enumerate() {
        let sample = channel_sample(index, symbol);
        let received_i = T::from_f64(sample[0]);
        let received_q = T::from_f64(sample[1]);

        // Rotate the received sample by the known negative carrier phase.
        let recovered_i = received_i * cosine + received_q * sine;
        let recovered_q = -received_i * sine + received_q * cosine;
        recovered[index] = [recovered_i.to_f64(), recovered_q.to_f64()];
        decoded[index] = decode(recovered[index][0], recovered[index][1]);
    }

    let symbol_errors = decoded
        .iter()
        .zip(PACKET)
        .filter(|(actual, expected)| **actual != *expected)
        .count();
    let bit_errors = decoded
        .iter()
        .zip(PACKET)
        .map(|(actual, expected)| (actual ^ expected).count_ones() as usize)
        .sum();
    let error_energy = recovered
        .iter()
        .zip(PACKET)
        .map(|(actual, symbol)| {
            let expected = symbol_coordinates(symbol);
            let error_i = actual[0] - expected[0];
            let error_q = actual[1] - expected[1];
            error_i * error_i + error_q * error_q
        })
        .sum::<f64>();
    let reference_energy = PACKET
        .into_iter()
        .map(|symbol| {
            let expected = symbol_coordinates(symbol);
            expected[0] * expected[0] + expected[1] * expected[1]
        })
        .sum::<f64>();
    let minimum_decision_margin = recovered
        .iter()
        .map(|point| point[0].abs().min(point[1].abs()))
        .fold(f64::INFINITY, f64::min);

    QpskResult {
        recovered,
        decoded,
        symbol_errors,
        bit_errors,
        evm: f64::sqrt(error_energy / reference_energy),
        minimum_decision_margin,
    }
}

fn evaluate_all() -> [QpskResult; PRECISION_COUNT] {
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

fn write_csv(results: &[QpskResult; PRECISION_COUNT], path: &Path) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    writeln!(
        output,
        "format,index,transmitted_symbol,received_i,received_q,recovered_i,recovered_q,decoded_symbol,symbol_errors,bit_errors,evm,minimum_decision_margin"
    )?;
    for (precision, result) in results.iter().enumerate() {
        for (index, symbol) in PACKET.into_iter().enumerate() {
            let received = channel_sample(index, symbol);
            writeln!(
                output,
                "{},{},{},{:.17e},{:.17e},{:.17e},{:.17e},{},{},{},{:.17e},{:.17e}",
                CSV_NAMES[precision],
                index,
                symbol,
                received[0],
                received[1],
                result.recovered[index][0],
                result.recovered[index][1],
                result.decoded[index],
                result.symbol_errors,
                result.bit_errors,
                result.evm,
                result.minimum_decision_margin,
            )?;
        }
    }
    output.flush()
}

fn draw_constellation(
    results: &[QpskResult; PRECISION_COUNT],
) -> Result<(), Box<dyn std::error::Error>> {
    let ideal_i = [-0.5, -0.5, 0.5, 0.5];
    let ideal_q = [-0.5, 0.5, -0.5, 0.5];
    let fp4_i = results[0]
        .recovered
        .iter()
        .map(|point| point[0])
        .collect::<Vec<_>>();
    let fp4_q = results[0]
        .recovered
        .iter()
        .map(|point| point[1])
        .collect::<Vec<_>>();
    let series = [
        XYSeries::new(&ideal_i, &ideal_q)
            .with_label("ideal symbols")
            .with_style(PlotStyle::Points),
        XYSeries::new(&fp4_i, &fp4_q)
            .with_label("FP4 recovered packet")
            .with_style(PlotStyle::Points),
    ];
    plot_series_with_config(
        &series,
        &PlotConfig::publication("QPSK carrier recovery")
            .with_font_scale(1.4)
            .with_x_label("in-phase component")
            .with_y_label("quadrature component")
            .with_save_path(CONSTELLATION_SVG_PATH),
    )?;
    Ok(())
}

pub(crate) fn run() -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(OUTPUT_DIR)?;
    let results = evaluate_all();
    write_csv(&results, Path::new(CSV_PATH))?;
    draw_constellation(&results)?;

    println!("wrote {CSV_PATH}");
    println!("wrote {CONSTELLATION_SVG_PATH}");
    println!(
        "symbols = {SYMBOL_COUNT}; phase offset = {CHANNEL_PHASE_RAD} rad; noise amplitude = {NOISE_AMPLITUDE}"
    );
    for (label, result) in LABELS.iter().zip(results) {
        println!(
            "{label:>9}: symbol errors = {:>2}, bit errors = {:>2}, EVM = {:.9e}, minimum margin = {:.9e}",
            result.symbol_errors, result.bit_errors, result.evm, result.minimum_decision_margin,
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_precision_decodes_the_packet_without_error() {
        for (label, result) in LABELS.iter().zip(evaluate_all()) {
            assert_eq!(result.symbol_errors, 0, "{label} symbol errors");
            assert_eq!(result.bit_errors, 0, "{label} bit errors");
            assert!(result.evm.is_finite(), "{label} non-finite EVM");
            assert!(
                result.minimum_decision_margin > 0.0,
                "{label} reached a decision boundary",
            );
        }
    }
}
