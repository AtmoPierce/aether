use std::env;
use std::path::{Path, PathBuf};

use aether_fpx_bench::{export::write_csv, perf::run_perf_stat_counters_only, types::Row};

const OPS: &[&str] = &[
    "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "exp", "expm1", "log",
    "log1p",
];
const MODES: &[&str] = &[
    "fpx_kernel_scalar",
    "fpx_kernel_fma",
    "fpx_kernel_simd",
    "fpx_kernel_simd_estrin",
    "fpx_kernel_raw_scalar",
    "fpx_kernel_raw_fma",
    "fpx_kernel_raw_simd",
    "fpx_kernel_raw_simd_estrin",
];
const PRECISIONS: &[&str] = &["f32", "f64"];
const EVENTS: &[&str] = &["instructions", "cycles", "branches", "branch-misses"];

struct Args {
    probe: PathBuf,
    out: PathBuf,
    runs: usize,
    iters: usize,
    ops: Vec<String>,
    precisions: Vec<String>,
    modes: Vec<String>,
}

fn parse_list(value: &str) -> Vec<String> {
    value
        .split(',')
        .map(str::trim)
        .filter(|item| !item.is_empty())
        .map(str::to_owned)
        .collect()
}

fn parse_args() -> Result<Args, String> {
    let mut out = PathBuf::from("target/kernel_path_perf.csv");
    let mut runs = 10_usize;
    let mut iters = 100_000_usize;
    let mut ops = OPS
        .iter()
        .map(|item| (*item).to_owned())
        .collect::<Vec<_>>();
    let mut precisions = PRECISIONS
        .iter()
        .map(|item| (*item).to_owned())
        .collect::<Vec<_>>();
    let mut modes = MODES
        .iter()
        .map(|item| (*item).to_owned())
        .collect::<Vec<_>>();
    let mut probe = None;
    let args = env::args().skip(1).collect::<Vec<_>>();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--out" => {
                index += 1;
                out = PathBuf::from(args.get(index).ok_or("missing value for --out")?);
            }
            "--runs" => {
                index += 1;
                runs = args
                    .get(index)
                    .ok_or("missing value for --runs")?
                    .parse()
                    .map_err(|_| "invalid --runs")?;
            }
            "--iters" => {
                index += 1;
                iters = args
                    .get(index)
                    .ok_or("missing value for --iters")?
                    .parse()
                    .map_err(|_| "invalid --iters")?;
            }
            "--ops" => {
                index += 1;
                ops = parse_list(args.get(index).ok_or("missing value for --ops")?);
            }
            "--precisions" => {
                index += 1;
                precisions = parse_list(args.get(index).ok_or("missing value for --precisions")?);
            }
            "--modes" => {
                index += 1;
                modes = parse_list(args.get(index).ok_or("missing value for --modes")?);
            }
            "-h" | "--help" => {
                println!(
                    "usage: kernel_path_perf [--runs N] [--iters N] [--ops comma,list] \
                     [--precisions f32,f64] [--modes comma,list] [--out CSV] \
                     <absolute trig_probe>"
                );
                std::process::exit(0);
            }
            value if value.starts_with('-') => return Err(format!("unknown argument: {value}")),
            value if probe.is_none() => probe = Some(PathBuf::from(value)),
            value => return Err(format!("unexpected argument: {value}")),
        }
        index += 1;
    }
    if runs == 0 || iters == 0 {
        return Err("--runs and --iters must be non-zero".into());
    }
    if ops.is_empty() || ops.iter().any(|op| !OPS.contains(&op.as_str())) {
        return Err(format!("--ops must contain only: {}", OPS.join(",")));
    }
    if precisions.is_empty()
        || precisions
            .iter()
            .any(|precision| !PRECISIONS.contains(&precision.as_str()))
    {
        return Err("--precisions must contain only f32,f64".into());
    }
    if modes.is_empty() || modes.iter().any(|mode| !MODES.contains(&mode.as_str())) {
        return Err(format!("--modes must contain only: {}", MODES.join(",")));
    }
    let probe = probe.ok_or("missing absolute trig_probe path")?;
    if !probe.is_absolute() || !probe.is_file() {
        return Err(format!(
            "trig_probe must be an existing absolute path: {}",
            probe.display()
        ));
    }
    Ok(Args {
        probe,
        out,
        runs,
        iters,
        ops,
        precisions,
        modes,
    })
}

fn write_rows(path: &Path, rows: &[Row]) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|error| error.to_string())?;
    }
    write_csv(path, rows)
}

fn main() -> Result<(), String> {
    let args = parse_args()?;
    let mut rows =
        Vec::with_capacity(args.ops.len() * args.modes.len() * args.precisions.len() * args.runs);
    println!(
        "Direct-kernel perf: {} ops x {} paths x {} precisions x {} runs; {} values/run",
        args.ops.len(),
        args.modes.len(),
        args.precisions.len(),
        args.runs,
        args.iters
    );

    for op in &args.ops {
        for mode in &args.modes {
            for precision in &args.precisions {
                for run in 0..args.runs {
                    let row = run_perf_stat_counters_only(
                        &args.probe,
                        op,
                        mode,
                        precision,
                        args.iters,
                        run,
                        EVENTS,
                    )
                    .map_err(|error| {
                        format!("perf failed for {op}/{mode}/{precision} run {run}: {error}")
                    })?;
                    let ns = row
                        .time_elapsed_s
                        .map(|seconds| seconds * 1.0e9 / args.iters as f64)
                        .unwrap_or(f64::NAN);
                    let cycles = row
                        .cycles
                        .map(|value| value / args.iters as f64)
                        .unwrap_or(f64::NAN);
                    let instructions = row
                        .instructions
                        .map(|value| value / args.iters as f64)
                        .unwrap_or(f64::NAN);
                    println!(
                        "{op:5} {precision} {mode:17} run={run:2} ns/value={ns:9.3} cycles/value={cycles:9.3} ins/value={instructions:9.3}"
                    );
                    rows.push(row);
                    write_rows(&args.out, &rows)?;
                }
            }
        }
    }
    println!("Wrote {}", args.out.display());
    Ok(())
}
