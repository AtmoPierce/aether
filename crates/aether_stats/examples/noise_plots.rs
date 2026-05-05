use std::env;
use std::fs;
use std::path::{Path, PathBuf};

use aether_fft::real_periodogram;
use aether_stats::{GaussMarkov1, RandomWalk, WhiteNoise};
use aether_viz::{plot_histograms, plot_psd, plot_series, HistogramSeries, XYSeries};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output_dir = output_dir_from_args();
    fs::create_dir_all(&output_dir)?;

    let sample_dt = 0.05_f64;
    let sample_count = 2_048_usize;
    let times = sample_times(sample_count, sample_dt);
    let sample_rate_hz = 1.0_f64 / sample_dt;

    let white_noise = sample_white_noise(sample_count);
    let random_walk = sample_random_walk(sample_count, sample_dt);
    let gauss_markov = sample_gauss_markov(sample_count, sample_dt);
    let random_walk_increments = successive_differences(&random_walk);
    let white_noise_psd = real_periodogram(&white_noise, sample_rate_hz)?;
    let random_walk_psd = real_periodogram(&random_walk, sample_rate_hz)?;
    let gauss_markov_psd = real_periodogram(&gauss_markov, sample_rate_hz)?;

    write_plot(
        &output_dir.join("white_noise.svg"),
        "White Noise",
        "time [s]",
        "value",
        &[XYSeries::new(&times, &white_noise).with_label("white noise")],
    )?;

    write_plot(
        &output_dir.join("random_walk.svg"),
        "Random Walk",
        "time [s]",
        "bias",
        &[XYSeries::new(&times, &random_walk).with_label("random walk")],
    )?;

    write_plot(
        &output_dir.join("gauss_markov.svg"),
        "First-Order Gauss-Markov",
        "time [s]",
        "bias",
        &[XYSeries::new(&times, &gauss_markov).with_label("gauss-markov")],
    )?;

    write_plot(
        &output_dir.join("noise_comparison.svg"),
        "Noise Process Comparison",
        "time [s]",
        "value",
        &[
            XYSeries::new(&times, &white_noise).with_label("white noise"),
            XYSeries::new(&times, &random_walk).with_label("random walk"),
            XYSeries::new(&times, &gauss_markov).with_label("gauss-markov"),
        ],
    )?;

    write_histogram(
        &output_dir.join("white_noise_histogram.svg"),
        "White Noise Histogram",
        "value",
        "count",
        &[HistogramSeries::new(&white_noise)
            .with_bins(50)
            .with_label("white noise")],
    )?;

    write_histogram(
        &output_dir.join("random_walk_increment_histogram.svg"),
        "Random Walk Increment Histogram",
        "value",
        "count",
        &[HistogramSeries::new(&random_walk_increments)
            .with_bins(50)
            .with_label("random walk increments")],
    )?;

    write_histogram(
        &output_dir.join("gauss_markov_histogram.svg"),
        "Gauss-Markov Histogram",
        "value",
        "count",
        &[HistogramSeries::new(&gauss_markov)
            .with_bins(50)
            .with_label("gauss-markov")],
    )?;

    write_histogram(
        &output_dir.join("noise_histogram_comparison.svg"),
        "Noise Histogram Comparison",
        "value",
        "count",
        &[
            HistogramSeries::new(&white_noise)
                .with_bins(50)
                .with_label("white noise"),
            HistogramSeries::new(&random_walk_increments)
                .with_bins(50)
                .with_label("random walk increments"),
            HistogramSeries::new(&gauss_markov)
                .with_bins(50)
                .with_label("gauss-markov"),
        ],
    )?;

    write_psd(
        &output_dir.join("white_noise_psd.svg"),
        "White Noise PSD",
        &[XYSeries::new(&white_noise_psd.frequencies_hz, &white_noise_psd.power_density)
            .with_label("white noise")],
    )?;

    write_psd(
        &output_dir.join("random_walk_psd.svg"),
        "Random Walk PSD",
        &[XYSeries::new(&random_walk_psd.frequencies_hz, &random_walk_psd.power_density)
            .with_label("random walk")],
    )?;

    write_psd(
        &output_dir.join("gauss_markov_psd.svg"),
        "Gauss-Markov PSD",
        &[XYSeries::new(&gauss_markov_psd.frequencies_hz, &gauss_markov_psd.power_density)
            .with_label("gauss-markov")],
    )?;

    write_psd(
        &output_dir.join("noise_psd_comparison.svg"),
        "Noise PSD Comparison",
        &[
            XYSeries::new(&white_noise_psd.frequencies_hz, &white_noise_psd.power_density)
                .with_label("white noise"),
            XYSeries::new(&random_walk_psd.frequencies_hz, &random_walk_psd.power_density)
                .with_label("random walk"),
            XYSeries::new(&gauss_markov_psd.frequencies_hz, &gauss_markov_psd.power_density)
                .with_label("gauss-markov"),
        ],
    )?;

    println!("wrote plots to {}", output_dir.display());
    Ok(())
}

fn output_dir_from_args() -> PathBuf {
    match env::args().nth(1) {
        Some(path) => PathBuf::from(path),
        None => PathBuf::from("target/aether_stats/noise_plots"),
    }
}

fn sample_times(sample_count: usize, sample_dt: f64) -> Vec<f64> {
    (0..sample_count)
        .map(|index| index as f64 * sample_dt)
        .collect()
}

fn sample_white_noise(sample_count: usize) -> Vec<f64> {
    let mut process = WhiteNoise::new(0.0_f64, 0.15_f64, 1);
    (0..sample_count).map(|_| process.sample()).collect()
}

fn sample_random_walk(sample_count: usize, sample_dt: f64) -> Vec<f64> {
    let mut process = RandomWalk::new(0.0_f64, 0.03_f64, 2);
    (0..sample_count)
        .map(|_| process.advance(sample_dt))
        .collect()
}

fn sample_gauss_markov(sample_count: usize, sample_dt: f64) -> Vec<f64> {
    let mut process = GaussMarkov1::new(0.0_f64, 8.0_f64, 0.25_f64, 3);
    (0..sample_count)
        .map(|_| process.advance(sample_dt))
        .collect()
}

fn successive_differences(values: &[f64]) -> Vec<f64> {
    values
        .windows(2)
        .map(|pair| pair[1] - pair[0])
        .collect()
}

fn write_plot(
    output_path: &Path,
    title: &str,
    x_label: &str,
    y_label: &str,
    series: &[XYSeries<'_>],
) -> Result<(), Box<dyn std::error::Error>> {
    let save_path = output_path.to_str().ok_or("output path is not valid UTF-8")?;
    plot_series(series, title, Some(x_label), Some(y_label), Some(save_path))
}

fn write_histogram(
    output_path: &Path,
    title: &str,
    x_label: &str,
    y_label: &str,
    series: &[HistogramSeries<'_>],
) -> Result<(), Box<dyn std::error::Error>> {
    let save_path = output_path.to_str().ok_or("output path is not valid UTF-8")?;
    plot_histograms(series, title, Some(x_label), Some(y_label), Some(save_path))
}

fn write_psd(
    output_path: &Path,
    title: &str,
    series: &[XYSeries<'_>],
) -> Result<(), Box<dyn std::error::Error>> {
    let save_path = output_path.to_str().ok_or("output path is not valid UTF-8")?;
    plot_psd(series, title, Some(save_path))
}