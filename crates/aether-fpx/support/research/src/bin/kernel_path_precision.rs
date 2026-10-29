use std::env;
use std::path::PathBuf;

use aether_fpx_bench::{valid_ulp_dist_f32, valid_ulp_dist_f64};
use aether_fpx::bench_paths::{
    eval_f32, eval_f64, kernel_range_f32, kernel_range_f64, ArithmeticPath, OPS,
};

const DEFAULT_SAMPLES: usize = 65_537;
const DEFAULT_OUT: &str = "target/kernel_path_precision.csv";

struct Args {
    samples: usize,
    out: PathBuf,
    ops: Vec<String>,
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
    let mut samples = DEFAULT_SAMPLES;
    let mut out = PathBuf::from(DEFAULT_OUT);
    let mut ops = OPS.iter().map(|op| (*op).to_string()).collect::<Vec<_>>();
    let args = env::args().skip(1).collect::<Vec<_>>();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--samples" => {
                index += 1;
                samples = args
                    .get(index)
                    .ok_or("missing value for --samples")?
                    .parse()
                    .map_err(|_| "invalid --samples")?;
            }
            "--out" => {
                index += 1;
                out = PathBuf::from(args.get(index).ok_or("missing value for --out")?);
            }
            "--ops" => {
                index += 1;
                ops = parse_list(args.get(index).ok_or("missing value for --ops")?);
            }
            "-h" | "--help" => {
                println!(
                    "usage: kernel_path_precision [--samples N] [--ops comma,list] [--out CSV]"
                );
                std::process::exit(0);
            }
            other => return Err(format!("unknown argument: {other}")),
        }
        index += 1;
    }
    if samples < 2 {
        return Err("--samples must be at least 2".into());
    }
    if ops.is_empty() || ops.iter().any(|op| !OPS.contains(&op.as_str())) {
        return Err(format!("--ops must contain only: {}", OPS.join(",")));
    }
    Ok(Args { samples, out, ops })
}

fn reference_f32(op: &str, x: f32) -> f32 {
    match op {
        "sin" => x.sin(),
        "cos" => x.cos(),
        "tan" => x.tan(),
        "asin" => x.asin(),
        "acos" => x.acos(),
        "atan" => x.atan(),
        "sinh" => x.sinh(),
        "cosh" => x.cosh(),
        "tanh" => x.tanh(),
        "exp" => x.exp(),
        "expm1" => x.exp_m1(),
        "log" => x.ln(),
        "log1p" => x.ln_1p(),
        _ => unreachable!(),
    }
}

fn reference_f64(op: &str, x: f64) -> f64 {
    match op {
        "sin" => x.sin(),
        "cos" => x.cos(),
        "tan" => x.tan(),
        "asin" => x.asin(),
        "acos" => x.acos(),
        "atan" => x.atan(),
        "sinh" => x.sinh(),
        "cosh" => x.cosh(),
        "tanh" => x.tanh(),
        "exp" => x.exp(),
        "expm1" => x.exp_m1(),
        "log" => x.ln(),
        "log1p" => x.ln_1p(),
        _ => unreachable!(),
    }
}

fn write_f32_rows(
    writer: &mut csv::Writer<std::fs::File>,
    op: &str,
    inputs: &[f32],
) -> Result<(), csv::Error> {
    let references = inputs
        .iter()
        .map(|&x| reference_f32(op, x))
        .collect::<Vec<_>>();
    let mut fma = inputs.to_vec();
    eval_f32(op, ArithmeticPath::Fma, &mut fma);

    for path in ArithmeticPath::ALL {
        let mut outputs = inputs.to_vec();
        eval_f32(op, path, &mut outputs);
        let mut max_abs = 0.0_f32;
        let mut max_abs_values = (0.0_f32, 0.0_f32, 0.0_f32);
        let mut max_ulp = 0_u32;
        let mut max_ulp_values = (0.0_f32, 0.0_f32, 0.0_f32);
        let mut sum_abs = 0.0_f64;
        let mut sum_ulp = 0_u128;
        let mut finite = 0_usize;

        for ((&x, &value), &reference) in inputs.iter().zip(&outputs).zip(&references) {
            if !value.is_finite() || !reference.is_finite() {
                continue;
            }
            finite += 1;
            let abs = (value - reference).abs();
            sum_abs += abs as f64;
            if abs > max_abs {
                max_abs = abs;
                max_abs_values = (x, reference, value);
            }
            if let Some(ulp) = valid_ulp_dist_f32(value, reference) {
                sum_ulp += ulp as u128;
                if ulp > max_ulp {
                    max_ulp = ulp;
                    max_ulp_values = (x, reference, value);
                }
            }
        }
        let mismatches = outputs
            .iter()
            .zip(&fma)
            .filter(|(value, baseline)| value.to_bits() != baseline.to_bits())
            .count();
        let mean_abs = sum_abs / finite.max(1) as f64;
        let mean_ulp = sum_ulp as f64 / finite.max(1) as f64;
        writer.write_record([
            op,
            "f32",
            path.name(),
            &inputs.len().to_string(),
            &finite.to_string(),
            &format!("{max_abs:.9e}"),
            &format!("{mean_abs:.9e}"),
            &max_ulp.to_string(),
            &format!("{mean_ulp:.9e}"),
            &mismatches.to_string(),
            &format!("{:.9e}", max_abs_values.0),
            &format!("{:.9e}", max_abs_values.1),
            &format!("{:.9e}", max_abs_values.2),
            &format!("{:.9e}", max_ulp_values.0),
            &format!("{:.9e}", max_ulp_values.1),
            &format!("{:.9e}", max_ulp_values.2),
        ])?;
        println!(
            "{op:5} f32 {:6}: max_abs={max_abs:.3e} max_ulp={max_ulp:<8} mismatch-vs-fma={mismatches}",
            path.name()
        );
    }
    Ok(())
}

fn write_f64_rows(
    writer: &mut csv::Writer<std::fs::File>,
    op: &str,
    inputs: &[f64],
) -> Result<(), csv::Error> {
    let references = inputs
        .iter()
        .map(|&x| reference_f64(op, x))
        .collect::<Vec<_>>();
    let mut fma = inputs.to_vec();
    eval_f64(op, ArithmeticPath::Fma, &mut fma);

    for path in ArithmeticPath::ALL {
        let mut outputs = inputs.to_vec();
        eval_f64(op, path, &mut outputs);
        let mut max_abs = 0.0_f64;
        let mut max_abs_values = (0.0_f64, 0.0_f64, 0.0_f64);
        let mut max_ulp = 0_u64;
        let mut max_ulp_values = (0.0_f64, 0.0_f64, 0.0_f64);
        let mut sum_abs = 0.0_f64;
        let mut sum_ulp = 0_u128;
        let mut finite = 0_usize;

        for ((&x, &value), &reference) in inputs.iter().zip(&outputs).zip(&references) {
            if !value.is_finite() || !reference.is_finite() {
                continue;
            }
            finite += 1;
            let abs = (value - reference).abs();
            sum_abs += abs;
            if abs > max_abs {
                max_abs = abs;
                max_abs_values = (x, reference, value);
            }
            if let Some(ulp) = valid_ulp_dist_f64(value, reference) {
                sum_ulp += ulp as u128;
                if ulp > max_ulp {
                    max_ulp = ulp;
                    max_ulp_values = (x, reference, value);
                }
            }
        }
        let mismatches = outputs
            .iter()
            .zip(&fma)
            .filter(|(value, baseline)| value.to_bits() != baseline.to_bits())
            .count();
        let mean_abs = sum_abs / finite.max(1) as f64;
        let mean_ulp = sum_ulp as f64 / finite.max(1) as f64;
        writer.write_record([
            op,
            "f64",
            path.name(),
            &inputs.len().to_string(),
            &finite.to_string(),
            &format!("{max_abs:.17e}"),
            &format!("{mean_abs:.17e}"),
            &max_ulp.to_string(),
            &format!("{mean_ulp:.17e}"),
            &mismatches.to_string(),
            &format!("{:.17e}", max_abs_values.0),
            &format!("{:.17e}", max_abs_values.1),
            &format!("{:.17e}", max_abs_values.2),
            &format!("{:.17e}", max_ulp_values.0),
            &format!("{:.17e}", max_ulp_values.1),
            &format!("{:.17e}", max_ulp_values.2),
        ])?;
        println!(
            "{op:5} f64 {:6}: max_abs={max_abs:.3e} max_ulp={max_ulp:<8} mismatch-vs-fma={mismatches}",
            path.name()
        );
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = parse_args().map_err(|error| {
        eprintln!("{error}");
        std::io::Error::new(std::io::ErrorKind::InvalidInput, error)
    })?;
    if let Some(parent) = args.out.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut writer = csv::Writer::from_path(&args.out)?;
    writer.write_record([
        "op",
        "precision",
        "path",
        "samples",
        "finite_samples",
        "max_abs",
        "mean_abs",
        "max_ulp",
        "mean_ulp",
        "bit_mismatches_vs_fma",
        "max_abs_x",
        "max_abs_ref",
        "max_abs_value",
        "max_ulp_x",
        "max_ulp_ref",
        "max_ulp_value",
    ])?;

    for op in &args.ops {
        let (lo, hi) = kernel_range_f32(op);
        let inputs = (0..args.samples)
            .map(|index| lo + (hi - lo) * index as f32 / (args.samples - 1) as f32)
            .collect::<Vec<_>>();
        write_f32_rows(&mut writer, op, &inputs)?;

        let (lo, hi) = kernel_range_f64(op);
        let inputs = (0..args.samples)
            .map(|index| lo + (hi - lo) * index as f64 / (args.samples - 1) as f64)
            .collect::<Vec<_>>();
        write_f64_rows(&mut writer, op, &inputs)?;
    }
    writer.flush()?;
    println!("Wrote {}", args.out.display());
    Ok(())
}
