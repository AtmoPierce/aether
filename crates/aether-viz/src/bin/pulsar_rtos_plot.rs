use std::error::Error;
use std::fs;
use std::io::{Error as IoError, ErrorKind};
use std::path::PathBuf;

use aether_viz::rtos::{
    compare_metric, load_rtos_run, plot_wake_latency_histogram, plot_wake_latency_series,
    DistributionStats, RtosRun,
};

fn usage_error() -> Box<dyn Error> {
    IoError::new(
        ErrorKind::InvalidInput,
        "usage: pulsar_rtos_plot <run-directory> <output-directory> [baseline-run-directory] [--max-p99-regression-percent N]",
    )
    .into()
}

fn report(run: &RtosRun) -> Result<(), Box<dyn Error>> {
    println!(
        "run: {} (origin: {})",
        run.run_id,
        run.metadata.get("data_origin").unwrap_or("unknown")
    );
    for warning in &run.warnings {
        eprintln!("warning [{}]: {warning}", run.run_id);
    }
    match run.wake_latency_summary()? {
        Some(summary) => {
            let (stats, unit) = summary
                .nanoseconds
                .map(|stats| (stats, "ns"))
                .unwrap_or((summary.ticks, "ticks"));
            print_stats(stats, unit);
            if let Some(duration_ns) = summary.capture_duration_ns {
                println!("  capture duration: {duration_ns:.0} ns");
            } else if let Some(duration_ticks) = summary.capture_duration_ticks {
                println!("  capture duration: {duration_ticks} ticks");
            }
        }
        None => println!("  no correlated wake-to-run samples"),
    }
    Ok(())
}

fn print_stats(stats: DistributionStats, unit: &str) {
    println!(
        "  wake-to-run: count={} min={:.3} median={:.3} p99={:.3} max={:.3} {}",
        stats.count, stats.min, stats.median, stats.p99, stats.max, unit
    );
}

fn compare(candidate: &RtosRun, baseline: &RtosRun) -> Result<Option<f64>, Box<dyn Error>> {
    let Some(candidate_summary) = candidate.wake_latency_summary()? else {
        return Ok(None);
    };
    let Some(baseline_summary) = baseline.wake_latency_summary()? else {
        return Ok(None);
    };
    let (candidate_stats, baseline_stats, unit) =
        match (candidate_summary.nanoseconds, baseline_summary.nanoseconds) {
            (Some(candidate), Some(baseline)) => (candidate, baseline, "ns"),
            _ => (candidate_summary.ticks, baseline_summary.ticks, "ticks"),
        };

    let mut p99_change_percent = None;
    for (name, candidate_value, baseline_value) in [
        ("median", candidate_stats.median, baseline_stats.median),
        ("p99", candidate_stats.p99, baseline_stats.p99),
        ("max", candidate_stats.max, baseline_stats.max),
    ] {
        let comparison = compare_metric(candidate_value, baseline_value);
        if name == "p99" {
            p99_change_percent = comparison.relative_change_percent;
        }
        match comparison.relative_change_percent {
            Some(percent) => println!(
                "comparison {name}: {:+.3} {unit} ({percent:+.2}%) candidate vs baseline",
                comparison.absolute_change
            ),
            None => println!(
                "comparison {name}: {:+.3} {unit}; percentage undefined for zero baseline",
                comparison.absolute_change
            ),
        }
    }
    Ok(p99_change_percent)
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut arguments = std::env::args_os().skip(1);
    let run_directory = PathBuf::from(arguments.next().ok_or_else(usage_error)?);
    let output_directory = PathBuf::from(arguments.next().ok_or_else(usage_error)?);
    let mut baseline_directory = None;
    let mut max_p99_regression_percent = None;
    while let Some(argument) = arguments.next() {
        if argument == "--max-p99-regression-percent" {
            let value = arguments.next().ok_or_else(usage_error)?;
            let value = value
                .to_str()
                .ok_or_else(usage_error)?
                .parse::<f64>()
                .map_err(|_| usage_error())?;
            if !value.is_finite() || value < 0.0 {
                return Err(usage_error());
            }
            max_p99_regression_percent = Some(value);
        } else if baseline_directory.is_none() {
            baseline_directory = Some(PathBuf::from(argument));
        } else {
            return Err(usage_error());
        }
    }

    let run = load_rtos_run(&run_directory)?;
    let baseline = baseline_directory
        .as_deref()
        .map(load_rtos_run)
        .transpose()?;
    fs::create_dir_all(&output_directory)?;

    plot_wake_latency_series(&run, output_directory.join("wake_latency.svg"))?;
    let candidate_label = format!("candidate: {}", run.run_id);
    if let Some(baseline) = &baseline {
        let baseline_label = format!("baseline: {}", baseline.run_id);
        plot_wake_latency_histogram(
            &[
                (&run, candidate_label.as_str()),
                (baseline, baseline_label.as_str()),
            ],
            output_directory.join("wake_latency_distribution.svg"),
        )?;
    } else {
        plot_wake_latency_histogram(
            &[(&run, candidate_label.as_str())],
            output_directory.join("wake_latency_distribution.svg"),
        )?;
    }

    report(&run)?;
    if let Some(baseline) = &baseline {
        report(baseline)?;
        let observed_p99_change = compare(&run, baseline)?;
        if let Some(limit) = max_p99_regression_percent {
            let observed = observed_p99_change.ok_or_else(|| {
                IoError::new(
                    ErrorKind::InvalidData,
                    "cannot enforce a percentage threshold with missing samples or a zero baseline",
                )
            })?;
            if observed > limit {
                return Err(IoError::other(format!(
                    "p99 regression {observed:.2}% exceeds configured limit {limit:.2}%"
                ))
                .into());
            }
            println!("threshold passed: p99 change {observed:+.2}% <= {limit:.2}%");
        }
    } else if max_p99_regression_percent.is_some() {
        return Err(IoError::new(
            ErrorKind::InvalidInput,
            "--max-p99-regression-percent requires a baseline run",
        )
        .into());
    }
    Ok(())
}
