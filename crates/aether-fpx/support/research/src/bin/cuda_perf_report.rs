use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
struct Args {
    csv_path: PathBuf,
    out_path: Option<PathBuf>,
    reference_csv_path: Option<PathBuf>,
    baseline_mode: String,
}

#[derive(Clone, Debug)]
struct PerfRow {
    op: String,
    mode: String,
    precision: String,
    run: usize,
    ns_per_launch: f64,
    ns_per_element: f64,
    elements_per_second: f64,
}

#[derive(Clone, Debug)]
struct SummaryRow {
    op: String,
    precision: String,
    mode: String,
    runs: usize,
    mean_ns_per_launch: f64,
    stdev_ns_per_launch: f64,
    min_ns_per_launch: f64,
    max_ns_per_launch: f64,
    mean_ns_per_element: f64,
    stdev_ns_per_element: f64,
    mean_elements_per_second: f64,
    stdev_elements_per_second: f64,
    baseline_mode: String,
    baseline_mean_ns_per_launch: Option<f64>,
    normalized_ns_per_launch: Option<f64>,
    speedup_baseline_over_mode: Option<f64>,
    reference_mean_ns_per_launch: Option<f64>,
    current_over_reference: Option<f64>,
    reference_over_current: Option<f64>,
}

fn fail(message: impl Into<String>) -> ! {
    eprintln!("cuda_perf_report: {}", message.into());
    std::process::exit(1);
}

fn print_usage() {
    println!(
        "usage: cargo cuda-perf-report -- --csv <path> [options]\n\n\
Summarize a raw CUDA performance CSV, normalize results against a baseline mode\n\
from the same CSV, and optionally compare against another CSV captured elsewhere.\n\n\
Options:\n\
    --csv <path>             Raw CUDA performance CSV to summarize\n\
  --out <path>             Summary CSV output path (default: <csv stem>_report.csv)\n\
  --reference-csv <path>   Optional second raw CSV to compare against\n\
  --baseline-mode <mode>   Mode used for normalization (default: nat)\n\
  --help                   Show this help message\n\n\
Examples:\n\
  cargo cuda-perf-report -- --csv target/cuda_perf_123.csv\n\
  cargo cuda-perf-report -- --csv target/cuda_perf_123.csv --reference-csv ~/other_host/cuda_perf_456.csv\n"
    );
}

fn parse_args() -> Args {
    let mut csv_path = None;
    let mut out_path = None;
    let mut reference_csv_path = None;
    let mut baseline_mode = "nat".to_string();

    let cli_args: Vec<String> = env::args().skip(1).collect();
    let mut index = 0usize;
    while index < cli_args.len() {
        let flag = &cli_args[index];
        let require_value = |name: &str, args: &[String], index: &mut usize| -> String {
            if *index + 1 >= args.len() {
                fail(format!("missing value for {name}"));
            }
            *index += 1;
            args[*index].clone()
        };

        match flag.as_str() {
            "--csv" => {
                csv_path = Some(PathBuf::from(require_value("--csv", &cli_args, &mut index)))
            }
            "--out" => {
                out_path = Some(PathBuf::from(require_value("--out", &cli_args, &mut index)))
            }
            "--reference-csv" => {
                reference_csv_path = Some(PathBuf::from(require_value(
                    "--reference-csv",
                    &cli_args,
                    &mut index,
                )))
            }
            "--baseline-mode" => {
                baseline_mode = require_value("--baseline-mode", &cli_args, &mut index)
            }
            "--help" | "-h" => {
                print_usage();
                std::process::exit(0);
            }
            _ => fail(format!("unknown argument: {flag}")),
        }

        index += 1;
    }

    let csv_path = csv_path.unwrap_or_else(|| fail("--csv is required"));
    Args {
        csv_path,
        out_path,
        reference_csv_path,
        baseline_mode,
    }
}

fn default_out_path(csv_path: &Path) -> PathBuf {
    let parent = csv_path.parent().unwrap_or_else(|| Path::new("."));
    let stem = csv_path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("cuda_perf");
    parent.join(format!("{stem}_report.csv"))
}

fn parse_f64(value: &str, field: &str, line_number: usize) -> f64 {
    value
        .parse::<f64>()
        .unwrap_or_else(|_| fail(format!("invalid {field} at line {line_number}: {value}")))
}

fn parse_usize(value: &str, field: &str, line_number: usize) -> usize {
    value
        .parse::<usize>()
        .unwrap_or_else(|_| fail(format!("invalid {field} at line {line_number}: {value}")))
}

fn parse_csv(path: &Path) -> Vec<PerfRow> {
    let text = fs::read_to_string(path)
        .unwrap_or_else(|err| fail(format!("failed to read {}: {err}", path.display())));

    let mut rows = Vec::new();
    for (line_index, line) in text.lines().enumerate() {
        if line_index == 0 {
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() != 11 {
            fail(format!(
                "expected 11 CSV fields at line {} in {}, found {}",
                line_index + 1,
                path.display(),
                fields.len()
            ));
        }

        rows.push(PerfRow {
            op: fields[0].to_string(),
            mode: fields[1].to_string(),
            precision: fields[2].to_string(),
            run: parse_usize(fields[3], "run", line_index + 1),
            ns_per_launch: parse_f64(fields[7], "ns_per_launch", line_index + 1),
            ns_per_element: parse_f64(fields[8], "ns_per_element", line_index + 1),
            elements_per_second: parse_f64(fields[9], "elements_per_second", line_index + 1),
        });
    }

    if rows.is_empty() {
        fail(format!(
            "{} did not contain any benchmark rows",
            path.display()
        ));
    }

    rows
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn stdev(values: &[f64], avg: f64) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }

    let sum_sq = values
        .iter()
        .map(|value| {
            let delta = value - avg;
            delta * delta
        })
        .sum::<f64>();
    (sum_sq / (values.len() as f64 - 1.0)).sqrt()
}

fn summarize(rows: &[PerfRow], baseline_mode: &str) -> Vec<SummaryRow> {
    type Key = (String, String, String);
    type BaseKey = (String, String);

    let mut grouped: BTreeMap<Key, Vec<&PerfRow>> = BTreeMap::new();
    for row in rows {
        grouped
            .entry((row.op.clone(), row.precision.clone(), row.mode.clone()))
            .or_default()
            .push(row);
    }

    let mut baseline_by_group: BTreeMap<BaseKey, f64> = BTreeMap::new();
    for ((op, precision, mode), group_rows) in &grouped {
        if mode == baseline_mode {
            let launches: Vec<f64> = group_rows.iter().map(|row| row.ns_per_launch).collect();
            baseline_by_group.insert((op.clone(), precision.clone()), mean(&launches));
        }
    }

    let mut summary_rows = Vec::new();
    for ((op, precision, mode), group_rows) in grouped {
        let runs: Vec<usize> = group_rows.iter().map(|row| row.run).collect();
        let launches: Vec<f64> = group_rows.iter().map(|row| row.ns_per_launch).collect();
        let elements: Vec<f64> = group_rows.iter().map(|row| row.ns_per_element).collect();
        let throughput: Vec<f64> = group_rows
            .iter()
            .map(|row| row.elements_per_second)
            .collect();

        let mean_launch = mean(&launches);
        let mean_element = mean(&elements);
        let mean_throughput = mean(&throughput);
        let baseline_mean = baseline_by_group
            .get(&(op.clone(), precision.clone()))
            .copied();
        let normalized = baseline_mean.map(|baseline| mean_launch / baseline);
        let speedup = baseline_mean.map(|baseline| baseline / mean_launch);

        summary_rows.push(SummaryRow {
            op,
            precision,
            mode,
            runs: runs.len(),
            mean_ns_per_launch: mean_launch,
            stdev_ns_per_launch: stdev(&launches, mean_launch),
            min_ns_per_launch: launches.iter().copied().fold(f64::INFINITY, f64::min),
            max_ns_per_launch: launches.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            mean_ns_per_element: mean_element,
            stdev_ns_per_element: stdev(&elements, mean_element),
            mean_elements_per_second: mean_throughput,
            stdev_elements_per_second: stdev(&throughput, mean_throughput),
            baseline_mode: baseline_mode.to_string(),
            baseline_mean_ns_per_launch: baseline_mean,
            normalized_ns_per_launch: normalized,
            speedup_baseline_over_mode: speedup,
            reference_mean_ns_per_launch: None,
            current_over_reference: None,
            reference_over_current: None,
        });
    }

    summary_rows
}

fn merge_reference(summary_rows: &mut [SummaryRow], reference_rows: &[PerfRow]) {
    let reference_summary = summarize(reference_rows, "nat");
    let mut reference_by_key: BTreeMap<(String, String, String), f64> = BTreeMap::new();
    for row in reference_summary {
        reference_by_key.insert((row.op, row.precision, row.mode), row.mean_ns_per_launch);
    }

    for row in summary_rows {
        if let Some(reference_mean) = reference_by_key
            .get(&(row.op.clone(), row.precision.clone(), row.mode.clone()))
            .copied()
        {
            row.reference_mean_ns_per_launch = Some(reference_mean);
            row.current_over_reference = Some(row.mean_ns_per_launch / reference_mean);
            row.reference_over_current = Some(reference_mean / row.mean_ns_per_launch);
        }
    }
}

fn fmt_opt(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.17}"))
        .unwrap_or_else(|| "NA".to_string())
}

fn print_summary(summary_rows: &[SummaryRow], reference_csv_path: Option<&Path>) {
    println!("Summary (mean +/- stdev ns/launch):");
    for row in summary_rows {
        println!(
            "op={} mode={} prec={} runs={} ns/launch={:.3} +/- {:.3} ns/element={:.6} speedup({}/{})={}",
            row.op,
            row.mode,
            row.precision,
            row.runs,
            row.mean_ns_per_launch,
            row.stdev_ns_per_launch,
            row.mean_ns_per_element,
            row.baseline_mode,
            row.mode,
            fmt_opt(row.speedup_baseline_over_mode),
        );
    }

    if let Some(reference_path) = reference_csv_path {
        println!(
            "\nReference comparison against {}:",
            reference_path.display()
        );
        for row in summary_rows {
            if row.reference_mean_ns_per_launch.is_some() {
                println!(
                    "op={} mode={} prec={} current/reference={} reference/current={}",
                    row.op,
                    row.mode,
                    row.precision,
                    fmt_opt(row.current_over_reference),
                    fmt_opt(row.reference_over_current),
                );
            }
        }
    }
}

fn write_summary_csv(path: &Path, summary_rows: &[SummaryRow]) {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).unwrap_or_else(|err| {
                fail(format!("failed to create {}: {err}", parent.display()))
            });
        }
    }

    let mut out = String::new();
    out.push_str(
        "op,precision,mode,runs,mean_ns_per_launch,stdev_ns_per_launch,min_ns_per_launch,max_ns_per_launch,mean_ns_per_element,stdev_ns_per_element,mean_elements_per_second,stdev_elements_per_second,baseline_mode,baseline_mean_ns_per_launch,normalized_ns_per_launch,speedup_baseline_over_mode,reference_mean_ns_per_launch,current_over_reference,reference_over_current\n",
    );

    for row in summary_rows {
        out.push_str(&format!(
            "{},{},{},{},{:.17},{:.17},{:.17},{:.17},{:.17},{:.17},{:.17},{:.17},{},{},{},{},{},{},{}\n",
            row.op,
            row.precision,
            row.mode,
            row.runs,
            row.mean_ns_per_launch,
            row.stdev_ns_per_launch,
            row.min_ns_per_launch,
            row.max_ns_per_launch,
            row.mean_ns_per_element,
            row.stdev_ns_per_element,
            row.mean_elements_per_second,
            row.stdev_elements_per_second,
            row.baseline_mode,
            fmt_opt(row.baseline_mean_ns_per_launch),
            fmt_opt(row.normalized_ns_per_launch),
            fmt_opt(row.speedup_baseline_over_mode),
            fmt_opt(row.reference_mean_ns_per_launch),
            fmt_opt(row.current_over_reference),
            fmt_opt(row.reference_over_current),
        ));
    }

    fs::write(path, out)
        .unwrap_or_else(|err| fail(format!("failed to write {}: {err}", path.display())));
}

fn main() {
    let args = parse_args();

    let rows = parse_csv(&args.csv_path);
    let mut summary_rows = summarize(&rows, &args.baseline_mode);
    if let Some(reference_csv_path) = args.reference_csv_path.as_deref() {
        let reference_rows = parse_csv(reference_csv_path);
        merge_reference(&mut summary_rows, &reference_rows);
    }

    let out_path = args
        .out_path
        .unwrap_or_else(|| default_out_path(&args.csv_path));
    write_summary_csv(&out_path, &summary_rows);
    print_summary(&summary_rows, args.reference_csv_path.as_deref());
    println!("\nWrote {}", out_path.display());
}
