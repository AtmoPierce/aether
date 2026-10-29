use std::fs;

const DEFAULT_INPUT: &str = "target/perf_trig_probe.csv";
const DEFAULT_ITERS: usize = 100_000;

const OPS: &[&str] = &[
    "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "exp", "expm1", "log",
    "log1p",
];
const PRECISIONS: &[&str] = &["f16", "f32", "f64", "f128"];

#[derive(Clone, Debug)]
struct PerfRow {
    op: String,
    mode: String,
    precision: String,
    time_elapsed_s: Option<f64>,
}

#[derive(Clone, Debug)]
struct RuntimeSummaryRow {
    op: String,
    precision: String,
    fpx_ns: Option<f64>,
    libm_ns: Option<f64>,
    nat_ns: Option<f64>,
    x87_ns: Option<f64>,
    fpx_vs_libm: Option<f64>,
    fpx_vs_nat: Option<f64>,
    fpx_vs_x87: Option<f64>,
}

#[derive(Clone, Debug)]
struct KernelSummaryRow {
    op: String,
    precision: String,
    horner_ns: Option<f64>,
    estrin_ns: Option<f64>,
    estrin_vs_horner: Option<f64>,
}

fn parse_opt_f64(text: &str) -> Option<f64> {
    let trimmed = text.trim();
    if trimmed.is_empty() {
        None
    } else {
        trimmed.parse().ok()
    }
}

fn read_perf_rows(path: &str) -> Result<Vec<PerfRow>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("failed to read {path}: {e}"))?;
    let headers: Vec<&str> = text
        .lines()
        .next()
        .unwrap_or_default()
        .split(',')
        .map(str::trim)
        .collect();
    let column = |name: &str| headers.iter().position(|header| *header == name);
    let op_column = column("op").ok_or_else(|| format!("CSV {path} is missing op column"))?;
    let mode_column = column("mode").ok_or_else(|| format!("CSV {path} is missing mode column"))?;
    let precision_column =
        column("precision").ok_or_else(|| format!("CSV {path} is missing precision column"))?;
    let elapsed_column = column("time_elapsed_s")
        .ok_or_else(|| format!("CSV {path} is missing time_elapsed_s column"))?;
    let required_columns = [op_column, mode_column, precision_column, elapsed_column];
    let min_columns = required_columns.into_iter().max().unwrap_or_default() + 1;
    let mut rows = Vec::new();

    for (line_idx, line) in text.lines().enumerate() {
        if line_idx == 0 || line.trim().is_empty() {
            continue;
        }

        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() < min_columns {
            return Err(format!(
                "malformed CSV row {} in {}: expected at least {} columns, got {}",
                line_idx + 1,
                path,
                min_columns,
                cols.len()
            ));
        }

        rows.push(PerfRow {
            op: cols[op_column].to_string(),
            mode: cols[mode_column].to_string(),
            precision: cols[precision_column].to_string(),
            time_elapsed_s: parse_opt_f64(cols[elapsed_column]),
        });
    }

    Ok(rows)
}

fn mean_ns_per_call(rows: &[&PerfRow], iters: usize) -> Option<f64> {
    let mut total_time = 0.0_f64;
    let mut samples = 0usize;

    for row in rows {
        if let Some(time) = row.time_elapsed_s {
            total_time += time;
            samples += 1;
        }
    }

    if samples == 0 {
        None
    } else {
        Some(1e9 * total_time / (iters as f64 * samples as f64))
    }
}

fn speedup(candidate_ns: Option<f64>, baseline_ns: Option<f64>) -> Option<f64> {
    match (candidate_ns, baseline_ns) {
        (Some(candidate), Some(baseline)) if candidate > 0.0 => Some(baseline / candidate),
        _ => None,
    }
}

fn collect_mode_mean(
    rows: &[PerfRow],
    op: &str,
    mode: &str,
    precision: &str,
    iters: usize,
) -> Option<f64> {
    let subset: Vec<&PerfRow> = rows
        .iter()
        .filter(|row| row.op == op && row.mode == mode && row.precision == precision)
        .collect();
    mean_ns_per_call(&subset, iters)
}

fn build_runtime_rows(rows: &[PerfRow], iters: usize) -> Vec<RuntimeSummaryRow> {
    let mut out = Vec::new();

    for &op in OPS {
        for &precision in PRECISIONS {
            let fpx_ns = collect_mode_mean(rows, op, "fpx", precision, iters);
            let libm_ns = collect_mode_mean(rows, op, "libm", precision, iters);
            let nat_ns = collect_mode_mean(rows, op, "nat", precision, iters);
            let x87_ns = collect_mode_mean(rows, op, "x87", precision, iters);

            if [fpx_ns, libm_ns, nat_ns, x87_ns]
                .iter()
                .all(|v| v.is_none())
            {
                continue;
            }

            out.push(RuntimeSummaryRow {
                op: op.to_string(),
                precision: precision.to_string(),
                fpx_ns,
                libm_ns,
                nat_ns,
                x87_ns,
                fpx_vs_libm: speedup(fpx_ns, libm_ns),
                fpx_vs_nat: speedup(fpx_ns, nat_ns),
                fpx_vs_x87: speedup(fpx_ns, x87_ns),
            });
        }
    }

    out
}

fn build_kernel_rows(rows: &[PerfRow], iters: usize) -> Vec<KernelSummaryRow> {
    let mut out = Vec::new();

    for &op in OPS {
        for &precision in PRECISIONS {
            let horner_ns = collect_mode_mean(rows, op, "fpx_kernel_horner", precision, iters);
            let estrin_ns = collect_mode_mean(rows, op, "fpx_kernel_estrin", precision, iters);

            if horner_ns.is_none() && estrin_ns.is_none() {
                continue;
            }

            out.push(KernelSummaryRow {
                op: op.to_string(),
                precision: precision.to_string(),
                horner_ns,
                estrin_ns,
                estrin_vs_horner: speedup(estrin_ns, horner_ns),
            });
        }
    }

    out
}

fn fmt_ns(value: Option<f64>) -> String {
    match value {
        Some(ns) => format!("{ns:.3}"),
        None => "NA".to_string(),
    }
}

fn fmt_speedup(value: Option<f64>) -> String {
    match value {
        Some(ratio) => format!("{ratio:.2}x"),
        None => "NA".to_string(),
    }
}

fn write_runtime_csv(path: &str, rows: &[RuntimeSummaryRow]) -> Result<(), String> {
    let mut out = String::from(
        "op,precision,fpx_ns_per_call,libm_ns_per_call,nat_ns_per_call,x87_ns_per_call,speedup_fpx_vs_libm,speedup_fpx_vs_nat,speedup_fpx_vs_x87\n",
    );

    for row in rows {
        out.push_str(&format!(
            "{},{},{},{},{},{},{},{},{}\n",
            row.op,
            row.precision,
            fmt_ns(row.fpx_ns),
            fmt_ns(row.libm_ns),
            fmt_ns(row.nat_ns),
            fmt_ns(row.x87_ns),
            fmt_speedup(row.fpx_vs_libm),
            fmt_speedup(row.fpx_vs_nat),
            fmt_speedup(row.fpx_vs_x87),
        ));
    }

    fs::write(path, out).map_err(|e| format!("failed to write {path}: {e}"))
}

fn write_kernel_csv(path: &str, rows: &[KernelSummaryRow]) -> Result<(), String> {
    let mut out = String::from(
        "op,precision,horner_ns_per_call,estrin_ns_per_call,speedup_estrin_vs_horner\n",
    );

    for row in rows {
        out.push_str(&format!(
            "{},{},{},{},{}\n",
            row.op,
            row.precision,
            fmt_ns(row.horner_ns),
            fmt_ns(row.estrin_ns),
            fmt_speedup(row.estrin_vs_horner),
        ));
    }

    fs::write(path, out).map_err(|e| format!("failed to write {path}: {e}"))
}

fn build_markdown(
    runtime_rows: &[RuntimeSummaryRow],
    kernel_rows: &[KernelSummaryRow],
    input_path: &str,
) -> String {
    let mut out = String::new();

    out.push_str("# Performance Summary\n\n");
    out.push_str(&format!("Source: `{}`\n\n", input_path));

    out.push_str("## Runtime\n\n");
    out.push_str("| op | prec | fpx ns/call | libm ns/call | nat ns/call | x87 ns/call | fpx vs libm | fpx vs nat | fpx vs x87 |\n");
    out.push_str("|---|---|---:|---:|---:|---:|---:|---:|---:|\n");
    for row in runtime_rows {
        out.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            row.op,
            row.precision,
            fmt_ns(row.fpx_ns),
            fmt_ns(row.libm_ns),
            fmt_ns(row.nat_ns),
            fmt_ns(row.x87_ns),
            fmt_speedup(row.fpx_vs_libm),
            fmt_speedup(row.fpx_vs_nat),
            fmt_speedup(row.fpx_vs_x87),
        ));
    }

    out.push_str("\n## Kernel Strategy\n\n");
    out.push_str("| op | prec | Horner ns/call | Estrin ns/call | Estrin vs Horner |\n");
    out.push_str("|---|---|---:|---:|---:|\n");
    for row in kernel_rows {
        out.push_str(&format!(
            "| {} | {} | {} | {} | {} |\n",
            row.op,
            row.precision,
            fmt_ns(row.horner_ns),
            fmt_ns(row.estrin_ns),
            fmt_speedup(row.estrin_vs_horner),
        ));
    }

    out
}

fn main() {
    let input_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| DEFAULT_INPUT.to_string());
    let rows = read_perf_rows(&input_path).unwrap_or_else(|err| {
        eprintln!("{err}");
        std::process::exit(1);
    });

    let runtime_rows = build_runtime_rows(&rows, DEFAULT_ITERS);
    let kernel_rows = build_kernel_rows(&rows, DEFAULT_ITERS);

    let runtime_csv_path = "target/perf_trig_probe_runtime_summary.csv";
    let kernel_csv_path = "target/perf_trig_probe_kernel_summary.csv";
    let markdown_path = "target/perf_trig_probe_simple.md";

    write_runtime_csv(runtime_csv_path, &runtime_rows).unwrap_or_else(|err| {
        eprintln!("{err}");
        std::process::exit(1);
    });

    write_kernel_csv(kernel_csv_path, &kernel_rows).unwrap_or_else(|err| {
        eprintln!("{err}");
        std::process::exit(1);
    });

    let markdown = build_markdown(&runtime_rows, &kernel_rows, &input_path);
    fs::write(markdown_path, &markdown).unwrap_or_else(|err| {
        eprintln!("failed to write {markdown_path}: {err}");
        std::process::exit(1);
    });

    println!("Wrote {runtime_csv_path}");
    println!("Wrote {kernel_csv_path}");
    println!("Wrote {markdown_path}\n");
    println!("{markdown}");
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_csv(name: &str, text: &str) -> String {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "fpx_simple_table_{name}_{}_{}.csv",
            std::process::id(),
            nanos
        ));
        fs::write(&path, text).unwrap();
        path.to_string_lossy().into_owned()
    }

    #[test]
    fn reads_perf_csv_with_optional_power_columns() {
        let path = temp_csv(
            "power_columns",
            "op,mode,precision,run,instructions,cycles,branches,branch_misses,time_elapsed_s,energy_joules,average_power_watts,power_source\n\
             tan,fpx,f16,0,3029332,4093497,278333,10824,0.053199966,,,\n",
        );

        let rows = read_perf_rows(&path).unwrap();
        let _ = fs::remove_file(path);

        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].op, "tan");
        assert_eq!(rows[0].mode, "fpx");
        assert_eq!(rows[0].precision, "f16");
        assert_eq!(rows[0].time_elapsed_s, Some(0.053199966));
    }
}
