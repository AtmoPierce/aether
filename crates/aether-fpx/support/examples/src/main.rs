#![feature(f16, f128)]

mod fft;
mod kepler;
mod linear_system;
mod perceptron;
mod qpsk;

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use aether_core::attitude::Quaternion;
use aether_core::math::Vector;
use aether_core::real::Real;
use aether_core::reference_frame::Unknown;
use aether_viz::{plot_series_with_config, AxisScale, PlotConfig, PlotStyle, XYSeries};
use aether_fpx::ocp::fp4::Fp4E2M1;
use aether_fpx::ocp::fp6::Fp6E3M2;
use aether_fpx::ocp::fp8::{Fp8E4M3, Fp8E5M2};

// A 1/2-second step is the finest shared cadence representable by every type:
// 0.5 is fp4 e2m1's smallest positive value.
const SAMPLE_COUNT: usize = 41;
const DURATION_S: f64 = 20.0;
const PRECISION_COUNT: usize = 8;
const PLOT_ERROR_FLOOR: f64 = 1.0e-18;

const PRECISION_LABELS: [&str; PRECISION_COUNT] = [
    "FPX fp4 e2m1",
    "FPX fp6 e3m2",
    "FPX fp8 e4m3",
    "FPX fp8 e5m2",
    "FPX f16",
    "FPX f32",
    "FPX f64",
    "FPX f128",
];
const PLOT_PRECISION_LABELS: [&str; PRECISION_COUNT] = [
    "FP4 E2M1", "FP6 E3M2", "FP8 E4M3", "FP8 E5M2", "f16", "f32", "f64", "f128",
];
const PRECISION_CSV_NAMES: [&str; PRECISION_COUNT] = [
    "fpx_fp4_e2m1",
    "fpx_fp6_e3m2",
    "fpx_fp8_e4m3",
    "fpx_fp8_e5m2",
    "fpx_f16",
    "fpx_f32",
    "fpx_f64",
    "fpx_f128",
];

const OUTPUT_DIR: &str = "outputs/rotational_dynamics";
const CSV_PATH: &str = "outputs/rotational_dynamics/rotational_dynamics.csv";
const ATTITUDE_SVG_PATH: &str = "outputs/rotational_dynamics/rotational_dynamics.svg";
const ATTITUDE_LINEAR_SVG_PATH: &str = "outputs/rotational_dynamics/rotational_dynamics_linear.svg";
const RATE_SVG_PATH: &str = "outputs/rotational_dynamics/angular_rate.svg";
const NORM_SVG_PATH: &str = "outputs/rotational_dynamics/quaternion_norm_error.svg";

type Attitude<T> = Quaternion<T, Unknown, Unknown>;

#[derive(Clone, Copy, Debug)]
struct Sample {
    time_s: f64,
    omega_mag_rad_s: f64,
    attitude_error_deg: [f64; PRECISION_COUNT],
    norm_error: [f64; PRECISION_COUNT],
}

struct FpxState<T: Real> {
    attitude: Attitude<T>,
    dt: T,
}

impl<T: Real> FpxState<T> {
    fn new(dt: f64) -> Self {
        Self {
            attitude: Attitude::identity(),
            dt: T::from_f64(dt),
        }
    }

    fn advance(&mut self, eval_time_s: f64) {
        let omega = angular_rate(T::from_f64(eval_time_s));
        let delta = delta_quaternion(omega, self.dt);
        self.attitude = self.attitude.mul_hamilton_same(&delta).normalized();
    }
}

struct FpxStates {
    fp4_e2m1: FpxState<Fp4E2M1>,
    fp6_e3m2: FpxState<Fp6E3M2>,
    fp8_e4m3: FpxState<Fp8E4M3>,
    fp8_e5m2: FpxState<Fp8E5M2>,
    f16: FpxState<f16>,
    f32: FpxState<f32>,
    f64: FpxState<f64>,
    f128: FpxState<f128>,
}

impl FpxStates {
    fn new(dt: f64) -> Self {
        Self {
            fp4_e2m1: FpxState::new(dt),
            fp6_e3m2: FpxState::new(dt),
            fp8_e4m3: FpxState::new(dt),
            fp8_e5m2: FpxState::new(dt),
            f16: FpxState::new(dt),
            f32: FpxState::new(dt),
            f64: FpxState::new(dt),
            f128: FpxState::new(dt),
        }
    }

    fn advance(&mut self, eval_time_s: f64) {
        self.fp4_e2m1.advance(eval_time_s);
        self.fp6_e3m2.advance(eval_time_s);
        self.fp8_e4m3.advance(eval_time_s);
        self.fp8_e5m2.advance(eval_time_s);
        self.f16.advance(eval_time_s);
        self.f32.advance(eval_time_s);
        self.f64.advance(eval_time_s);
        self.f128.advance(eval_time_s);
    }

    fn attitude_errors(&self, reference: &Attitude<f64>) -> [f64; PRECISION_COUNT] {
        [
            attitude_error_deg(reference, &self.fp4_e2m1.attitude),
            attitude_error_deg(reference, &self.fp6_e3m2.attitude),
            attitude_error_deg(reference, &self.fp8_e4m3.attitude),
            attitude_error_deg(reference, &self.fp8_e5m2.attitude),
            attitude_error_deg(reference, &self.f16.attitude),
            attitude_error_deg(reference, &self.f32.attitude),
            attitude_error_deg(reference, &self.f64.attitude),
            attitude_error_deg(reference, &self.f128.attitude),
        ]
    }

    fn norm_errors(&self) -> [f64; PRECISION_COUNT] {
        [
            norm_error(&self.fp4_e2m1.attitude),
            norm_error(&self.fp6_e3m2.attitude),
            norm_error(&self.fp8_e4m3.attitude),
            norm_error(&self.fp8_e5m2.attitude),
            norm_error(&self.f16.attitude),
            norm_error(&self.f32.attitude),
            norm_error(&self.f64.attitude),
            norm_error(&self.f128.attitude),
        ]
    }
}

// A single generic model exercises Aether's FPX-backed Real implementation for
// fp4 through f128. The changing axis prevents the trajectory from collapsing
// into a one-dimensional rotation.
fn angular_rate<T: Real>(time_s: T) -> Vector<T, 3> {
    Vector::new([
        T::from_f64(0.8) * (T::from_f64(0.7) * time_s).sin() + T::from_f64(0.2),
        T::from_f64(0.5) * (T::from_f64(0.3) * time_s).cos(),
        T::from_f64(0.4) * (T::from_f64(1.1) * time_s).sin() + T::from_f64(0.1),
    ])
}

fn angular_rate_reference(time_s: f64) -> Vector<f64, 3> {
    Vector::new([
        0.8 * f64::sin(0.7 * time_s) + 0.2,
        0.5 * f64::cos(0.3 * time_s),
        0.4 * f64::sin(1.1 * time_s) + 0.1,
    ])
}

fn delta_quaternion<T: Real>(omega: Vector<T, 3>, dt: T) -> Attitude<T> {
    let omega_mag = omega.norm();
    if omega_mag <= T::EPSILON {
        return Attitude::identity();
    }

    let half_theta = omega_mag * dt * T::from_f64(0.5);
    let scale = half_theta.sin() / omega_mag;
    Attitude::new(
        half_theta.cos(),
        omega[0] * scale,
        omega[1] * scale,
        omega[2] * scale,
    )
}

fn delta_quaternion_reference(omega: Vector<f64, 3>, dt: f64) -> Attitude<f64> {
    let omega_mag = omega.norm();
    if omega_mag <= f64::EPSILON {
        return Attitude::identity();
    }

    let half_theta = omega_mag * dt * 0.5;
    let scale = f64::sin(half_theta) / omega_mag;
    Attitude::new(
        f64::cos(half_theta),
        omega[0] * scale,
        omega[1] * scale,
        omega[2] * scale,
    )
}

fn attitude_error_deg<T: Real>(reference: &Attitude<f64>, test: &Attitude<T>) -> f64 {
    let reference_components = [reference.w(), reference.i(), reference.j(), reference.k()];
    let test_components = [
        test.w().to_f64(),
        test.i().to_f64(),
        test.j().to_f64(),
        test.k().to_f64(),
    ];

    // q and -q encode the same attitude. The smaller chord distance is
    // 2*sin(error/4), which stays well conditioned when the attitudes agree.
    let direct_distance_squared = reference_components
        .iter()
        .zip(test_components)
        .map(|(reference, test)| (reference - test) * (reference - test))
        .sum::<f64>();
    let antipodal_distance_squared = reference_components
        .iter()
        .zip(test_components)
        .map(|(reference, test)| (reference + test) * (reference + test))
        .sum::<f64>();
    let chord = f64::sqrt(direct_distance_squared.min(antipodal_distance_squared));
    f64::to_degrees(4.0 * f64::asin((0.5 * chord).clamp(0.0, 1.0)))
}

fn norm_error<T: Real>(attitude: &Attitude<T>) -> f64 {
    (attitude.norm() - T::ONE).abs().to_f64()
}

fn make_sample(
    time_s: f64,
    omega_mag_rad_s: f64,
    reference: &Attitude<f64>,
    states: &FpxStates,
) -> Sample {
    Sample {
        time_s,
        omega_mag_rad_s,
        attitude_error_deg: states.attitude_errors(reference),
        norm_error: states.norm_errors(),
    }
}

fn propagate() -> Vec<Sample> {
    let dt = DURATION_S / (SAMPLE_COUNT - 1) as f64;
    let mut reference = Attitude::<f64>::identity();
    let mut states = FpxStates::new(dt);

    let mut samples = Vec::with_capacity(SAMPLE_COUNT);
    samples.push(make_sample(
        0.0,
        angular_rate_reference(0.0).norm(),
        &reference,
        &states,
    ));

    for index in 1..SAMPLE_COUNT {
        let time_s = index as f64 * dt;
        let eval_time_s = time_s - 0.5 * dt;

        let reference_omega = angular_rate_reference(eval_time_s);
        let reference_delta = delta_quaternion_reference(reference_omega, dt);
        reference = reference.mul_hamilton_same(&reference_delta).normalized();

        states.advance(eval_time_s);

        samples.push(make_sample(
            time_s,
            reference_omega.norm(),
            &reference,
            &states,
        ));
    }

    samples
}

fn write_csv(samples: &[Sample], csv_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(csv_path)?;
    let mut writer = BufWriter::new(file);
    write!(writer, "time_s,omega_mag_rad_s")?;
    for name in PRECISION_CSV_NAMES {
        write!(writer, ",{name}_error_deg")?;
    }
    for name in PRECISION_CSV_NAMES {
        write!(writer, ",{name}_norm_error")?;
    }
    writeln!(writer)?;

    for sample in samples {
        write!(writer, "{:.9},{:.9}", sample.time_s, sample.omega_mag_rad_s)?;
        for value in sample.attitude_error_deg {
            write!(writer, ",{value:.12e}")?;
        }
        for value in sample.norm_error {
            write!(writer, ",{value:.12e}")?;
        }
        writeln!(writer)?;
    }

    writer.flush()?;
    Ok(())
}

fn draw_plots(samples: &[Sample]) -> Result<(), Box<dyn std::error::Error>> {
    let times = samples
        .iter()
        .map(|sample| sample.time_s)
        .collect::<Vec<_>>();
    let rates = samples
        .iter()
        .map(|sample| sample.omega_mag_rad_s)
        .collect::<Vec<_>>();
    let attitude_errors: [Vec<f64>; PRECISION_COUNT] = std::array::from_fn(|precision| {
        samples
            .iter()
            .map(|sample| sample.attitude_error_deg[precision].max(PLOT_ERROR_FLOOR))
            .collect()
    });
    let norm_errors: [Vec<f64>; PRECISION_COUNT] = std::array::from_fn(|precision| {
        samples
            .iter()
            .map(|sample| sample.norm_error[precision].max(PLOT_ERROR_FLOOR))
            .collect()
    });

    plot_series_with_config(
        &[XYSeries::new(&times, &rates)
            .with_label("native f64 reference")
            .with_style(PlotStyle::LineWithPalette(0))],
        &PlotConfig::publication("Time-varying body angular rate")
            .with_font_scale(1.6)
            .with_x_label("time [s]")
            .with_y_label("|omega| [rad/s]")
            .with_save_path(RATE_SVG_PATH),
    )?;

    let attitude_series = std::array::from_fn::<_, PRECISION_COUNT, _>(|precision| {
        XYSeries::new(&times, &attitude_errors[precision])
            .with_label(PLOT_PRECISION_LABELS[precision])
            .with_style(PlotStyle::LineWithPalette(precision))
    });
    plot_series_with_config(
        &attitude_series,
        &PlotConfig::publication("Attitude error vs native f64 (log scale)")
            .with_font_scale(1.6)
            .with_x_label("time [s]")
            .with_y_label("attitude error [deg]")
            .with_save_path(ATTITUDE_SVG_PATH)
            .with_y_scale(AxisScale::Log10),
    )?;
    plot_series_with_config(
        &attitude_series,
        &PlotConfig::publication("Attitude error vs native f64 (linear scale)")
            .with_font_scale(1.6)
            .with_x_label("time [s]")
            .with_y_label("attitude error [deg]")
            .with_save_path(ATTITUDE_LINEAR_SVG_PATH),
    )?;

    let norm_series = std::array::from_fn::<_, PRECISION_COUNT, _>(|precision| {
        XYSeries::new(&times, &norm_errors[precision])
            .with_label(PLOT_PRECISION_LABELS[precision])
            .with_style(PlotStyle::LineWithPalette(precision))
    });
    plot_series_with_config(
        &norm_series,
        &PlotConfig::publication("Quaternion normalization error")
            .with_font_scale(1.6)
            .with_x_label("time [s]")
            .with_y_label("abs(|q| - 1)")
            .with_save_path(NORM_SVG_PATH)
            .with_y_scale(AxisScale::Log10),
    )?;

    Ok(())
}

fn max_for_precision(samples: &[Sample], precision: usize) -> f64 {
    samples
        .iter()
        .map(|sample| sample.attitude_error_deg[precision])
        .fold(0.0, f64::max)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(OUTPUT_DIR)?;

    let samples = propagate();
    write_csv(&samples, Path::new(CSV_PATH))?;
    draw_plots(&samples)?;

    let final_sample = samples.last().expect("trajectory must contain samples");
    println!("wrote {CSV_PATH}");
    println!("wrote {ATTITUDE_SVG_PATH}");
    println!("wrote {ATTITUDE_LINEAR_SVG_PATH}");
    println!("wrote {RATE_SVG_PATH}");
    println!("wrote {NORM_SVG_PATH}");
    println!();
    println!("duration = {DURATION_S:.3} s");
    println!("samples  = {SAMPLE_COUNT}");
    println!("dt       = {:.9} s", DURATION_S / (SAMPLE_COUNT - 1) as f64);
    println!();

    for (precision, label) in PRECISION_LABELS.iter().enumerate() {
        println!(
            "{label:>13}: max error = {:.9e} deg, final error = {:.9e} deg",
            max_for_precision(&samples, precision),
            final_sample.attitude_error_deg[precision],
        );
    }

    println!();
    perceptron::run()?;
    println!();
    fft::run()?;
    println!();
    kepler::run()?;
    println!();
    linear_system::run()?;
    println!();
    qpsk::run()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_fpx_precisions_produce_finite_samples() {
        let samples = propagate();
        assert_eq!(samples.len(), SAMPLE_COUNT);
        for sample in samples {
            assert!(sample.time_s.is_finite());
            assert!(sample.omega_mag_rad_s.is_finite());
            for (precision, label) in PRECISION_LABELS.iter().enumerate() {
                assert!(
                    sample.attitude_error_deg[precision].is_finite(),
                    "{label} attitude error became non-finite at t={}",
                    sample.time_s,
                );
                assert!(
                    sample.norm_error[precision].is_finite(),
                    "{label} norm error became non-finite at t={}",
                    sample.time_s,
                );
            }
        }
    }
}
