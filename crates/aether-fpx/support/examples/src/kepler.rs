use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use aether_core::real::Real;
use aether_viz::{plot_series_with_config, PlotConfig, PlotStyle, XYSeries};
use aether_fpx::ocp::fp4::Fp4E2M1;
use aether_fpx::ocp::fp6::Fp6E3M2;
use aether_fpx::ocp::fp8::{Fp8E4M3, Fp8E5M2};

const PRECISION_COUNT: usize = 8;
const SAMPLE_COUNT: usize = 6;
const ITERATIONS: usize = 8;
const ECCENTRICITY: f64 = 0.5;
const MEAN_ANOMALIES: [f64; SAMPLE_COUNT] = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0];
const TRACE_SAMPLE: usize = 3;
const OUTPUT_DIR: &str = "outputs/kepler";
const CSV_PATH: &str = "outputs/kepler/kepler.csv";
const ERROR_SVG_PATH: &str = "outputs/kepler/error.svg";
const CONVERGENCE_SVG_PATH: &str = "outputs/kepler/convergence.svg";

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
struct KeplerResult {
    eccentric_anomaly: [f64; SAMPLE_COUNT],
    solution_error: [f64; SAMPLE_COUNT],
    residual: [f64; SAMPLE_COUNT],
    trace_error: [f64; ITERATIONS + 1],
}

fn reference_solution(mean_anomaly: f64) -> f64 {
    let mut eccentric_anomaly = mean_anomaly;
    for _ in 0..32 {
        let residual =
            eccentric_anomaly - ECCENTRICITY * f64::sin(eccentric_anomaly) - mean_anomaly;
        let derivative = 1.0 - ECCENTRICITY * f64::cos(eccentric_anomaly);
        eccentric_anomaly -= residual / derivative;
    }
    eccentric_anomaly
}

fn solve<T: Real + Copy>(mean_anomaly: f64) -> [f64; ITERATIONS + 1] {
    let mean_anomaly = T::from_f64(mean_anomaly);
    let eccentricity = T::from_f64(ECCENTRICITY);
    let mut eccentric_anomaly = mean_anomaly;
    let mut trace = [0.0; ITERATIONS + 1];
    trace[0] = eccentric_anomaly.to_f64();

    for iteration in 0..ITERATIONS {
        let residual = eccentric_anomaly - eccentricity * eccentric_anomaly.sin() - mean_anomaly;
        let derivative = T::ONE - eccentricity * eccentric_anomaly.cos();
        eccentric_anomaly = eccentric_anomaly - residual / derivative;
        trace[iteration + 1] = eccentric_anomaly.to_f64();
    }
    trace
}

fn evaluate<T: Real + Copy>() -> KeplerResult {
    let mut eccentric_anomaly = [0.0; SAMPLE_COUNT];
    let mut solution_error = [0.0; SAMPLE_COUNT];
    let mut residual = [0.0; SAMPLE_COUNT];
    let mut trace_error = [0.0; ITERATIONS + 1];

    for (sample, mean_anomaly) in MEAN_ANOMALIES.into_iter().enumerate() {
        let trace = solve::<T>(mean_anomaly);
        let estimate = trace[ITERATIONS];
        let reference = reference_solution(mean_anomaly);
        eccentric_anomaly[sample] = estimate;
        solution_error[sample] = (estimate - reference).abs();
        residual[sample] = (estimate - ECCENTRICITY * f64::sin(estimate) - mean_anomaly).abs();

        if sample == TRACE_SAMPLE {
            trace_error = trace.map(|value| (value - reference).abs());
        }
    }

    KeplerResult {
        eccentric_anomaly,
        solution_error,
        residual,
        trace_error,
    }
}

fn evaluate_all() -> [KeplerResult; PRECISION_COUNT] {
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

fn write_csv(results: &[KeplerResult; PRECISION_COUNT], path: &Path) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    write!(
        output,
        "format,mean_anomaly,eccentricity,eccentric_anomaly,solution_error,residual"
    )?;
    for iteration in 0..=ITERATIONS {
        write!(output, ",trace_error_iteration_{iteration}")?;
    }
    writeln!(output)?;

    for (precision, result) in results.iter().enumerate() {
        for sample in 0..SAMPLE_COUNT {
            write!(
                output,
                "{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                CSV_NAMES[precision],
                MEAN_ANOMALIES[sample],
                ECCENTRICITY,
                result.eccentric_anomaly[sample],
                result.solution_error[sample],
                result.residual[sample],
            )?;
            for error in result.trace_error {
                write!(output, ",{error:.17e}")?;
            }
            writeln!(output)?;
        }
    }
    output.flush()
}

fn draw_error_plot(
    results: &[KeplerResult; PRECISION_COUNT],
) -> Result<(), Box<dyn std::error::Error>> {
    let errors: [Vec<f64>; PRECISION_COUNT] =
        std::array::from_fn(|precision| results[precision].solution_error.into_iter().collect());
    let series = std::array::from_fn::<_, PRECISION_COUNT, _>(|precision| {
        XYSeries::new(&MEAN_ANOMALIES, &errors[precision])
            .with_label(LABELS[precision])
            .with_style(PlotStyle::LineWithPalette(precision))
    });
    plot_series_with_config(
        &series,
        &PlotConfig::publication("Kepler solution error")
            .with_font_scale(1.4)
            .with_x_label("mean anomaly M [rad]")
            .with_y_label("absolute E error [rad]")
            .with_save_path(ERROR_SVG_PATH),
    )?;
    Ok(())
}

fn draw_convergence_plot(
    results: &[KeplerResult; PRECISION_COUNT],
) -> Result<(), Box<dyn std::error::Error>> {
    let iterations = (0..=ITERATIONS)
        .map(|iteration| iteration as f64)
        .collect::<Vec<_>>();
    let errors: [Vec<f64>; PRECISION_COUNT] =
        std::array::from_fn(|precision| results[precision].trace_error.into_iter().collect());
    let series = std::array::from_fn::<_, PRECISION_COUNT, _>(|precision| {
        XYSeries::new(&iterations, &errors[precision])
            .with_label(LABELS[precision])
            .with_style(PlotStyle::LineWithPalette(precision))
    });
    plot_series_with_config(
        &series,
        &PlotConfig::publication("Newton convergence: e = 0.5, M = 1.5")
            .with_font_scale(1.4)
            .with_x_label("Newton iteration")
            .with_y_label("absolute E error [rad]")
            .with_save_path(CONVERGENCE_SVG_PATH),
    )?;
    Ok(())
}

pub(crate) fn run() -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(OUTPUT_DIR)?;
    let results = evaluate_all();
    write_csv(&results, Path::new(CSV_PATH))?;
    draw_error_plot(&results)?;
    draw_convergence_plot(&results)?;

    println!("wrote {CSV_PATH}");
    println!("wrote {ERROR_SVG_PATH}");
    println!("wrote {CONVERGENCE_SVG_PATH}");
    println!("e = {ECCENTRICITY}; Newton iterations = {ITERATIONS}");
    for (label, result) in LABELS.iter().zip(results) {
        let max_error = result.solution_error.into_iter().fold(0.0, f64::max);
        let max_residual = result.residual.into_iter().fold(0.0, f64::max);
        println!(
            "{label:>9}: max solution error = {max_error:.9e}, max residual = {max_residual:.9e}, E(M=1.5) = {:.9e}",
            result.eccentric_anomaly[TRACE_SAMPLE],
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_precision_produces_finite_kepler_solutions() {
        for (label, result) in LABELS.iter().zip(evaluate_all()) {
            assert!(
                result.eccentric_anomaly.into_iter().all(f64::is_finite)
                    && result.solution_error.into_iter().all(f64::is_finite)
                    && result.residual.into_iter().all(f64::is_finite)
                    && result.trace_error.into_iter().all(f64::is_finite),
                "{label} produced a non-finite Kepler result",
            );
        }
    }
}
