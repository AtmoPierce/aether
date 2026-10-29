use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use aether_fpx_bench::stats::stats;

const DEFAULT_CSV: &str = "target/perf_trig_probe.csv";
const DEFAULT_OUT: &str = "docs/figures/perf_speed_tables.tex";
const DEFAULT_ITERS: f64 = 100_000.0;
const DEFAULT_MODE_ORDER: &[&str] = &[
    "nat",
    "fpx",
    "fpx_kernel_horner",
    "fpx_kernel_estrin",
    "fpx_kernel_scalar",
    "fpx_kernel_fma",
    "fpx_kernel_simd",
    "fpx_kernel_simd_estrin",
    "fpx_kernel_raw_scalar",
    "fpx_kernel_raw_fma",
    "fpx_kernel_raw_simd",
    "fpx_kernel_raw_simd_estrin",
    "libm",
    "x87",
];
const DEFAULT_PRECISION_ORDER: &[&str] = &["f16", "f32", "f64", "f128"];
const MAX_METRIC_ROWS_PER_TABLE: usize = 12;
const CANONICAL_OPS: &[&str] = &[
    "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "exp", "expm1", "log",
    "log1p",
];

#[derive(Debug, Clone)]
struct Row {
    op: String,
    mode: String,
    precision: String,
    iterations: f64,
    ns_per_call: f64,
    cyc_per_call: Option<f64>,
    ins_per_call: Option<f64>,
    energy_nj_per_call: Option<f64>,
    average_power_watts: Option<f64>,
    power_source: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PrecisionRow {
    op: String,
    precision: String,
    max_ulp: Option<u128>,
}

#[derive(Debug, Clone, Copy)]
enum Metric {
    NsPerCall,
    CyclesPerCall,
    InstructionsPerCall,
    EnergyNjPerCall,
    AveragePowerWatts,
}

#[derive(Debug, Clone, Copy)]
enum CoverageKind {
    Speed,
    Precision,
    Power,
}

#[derive(Debug, Clone)]
struct Args {
    csv: PathBuf,
    precision_csv: Option<PathBuf>,
    out: PathBuf,
    iters: f64,
    target: Option<String>,
    host_tag: Option<String>,
    modes: Option<Vec<String>>,
    precisions: Option<Vec<String>>,
    ops: Option<Vec<String>>,
    coverage: bool,
    strict_complete: bool,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = parse_args(env::args().skip(1).collect())?;
    let rows = load_rows(&args.csv, args.iters)?;
    if rows.is_empty() {
        return Err(format!("no usable rows found in {}", args.csv.display()).into());
    }
    let precision_rows = match &args.precision_csv {
        Some(path) => load_precision_rows(path)?,
        None => Vec::new(),
    };
    if args.strict_complete {
        validate_complete_coverage(&rows, &precision_rows)?;
    }

    let latex = render_report(&rows, &precision_rows, &args)?;
    if let Some(parent) = args.out.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(&args.out, format!("{}\n", latex.trim_end()))?;
    println!("Wrote {}", args.out.display());
    Ok(())
}

fn parse_args(args: Vec<String>) -> Result<Args, Box<dyn Error>> {
    let mut csv = PathBuf::from(DEFAULT_CSV);
    let mut precision_csv = None;
    let mut out = PathBuf::from(DEFAULT_OUT);
    let mut iters = DEFAULT_ITERS;
    let mut target = None;
    let mut host_tag = None;
    let mut modes = None;
    let mut precisions = Some(
        DEFAULT_PRECISION_ORDER
            .iter()
            .map(|value| (*value).to_string())
            .collect::<Vec<_>>(),
    );
    let mut ops = None;
    let mut coverage = false;
    let mut strict_complete = false;

    let mut index = 0usize;
    while index < args.len() {
        match args[index].as_str() {
            "--csv" => {
                index += 1;
                csv = PathBuf::from(required_value(&args, index, "--csv")?);
            }
            "--precision-csv" => {
                index += 1;
                precision_csv = Some(PathBuf::from(required_value(
                    &args,
                    index,
                    "--precision-csv",
                )?));
            }
            "--out" => {
                index += 1;
                out = PathBuf::from(required_value(&args, index, "--out")?);
            }
            "--iters" => {
                index += 1;
                iters = required_value(&args, index, "--iters")?.parse()?;
                if !iters.is_finite() || iters <= 0.0 {
                    return Err("--iters must be a positive finite number".into());
                }
            }
            "--target" => {
                index += 1;
                target = Some(required_value(&args, index, "--target")?.to_string());
            }
            "--host-tag" | "--label-prefix" => {
                let flag = args[index].clone();
                index += 1;
                if host_tag.is_some() {
                    return Err("specify only one of --host-tag or --label-prefix".into());
                }
                let tag = slugify(required_value(&args, index, &flag)?);
                if tag.is_empty() {
                    return Err(format!("{flag} must contain at least one letter or number").into());
                }
                host_tag = Some(tag);
            }
            "--modes" => {
                index += 1;
                modes = parse_filter(required_value(&args, index, "--modes")?);
            }
            "--precisions" => {
                index += 1;
                precisions = parse_filter(required_value(&args, index, "--precisions")?);
            }
            "--ops" => {
                index += 1;
                ops = parse_filter(required_value(&args, index, "--ops")?);
            }
            "--coverage" => coverage = true,
            "--strict-complete" => strict_complete = true,
            "-h" | "--help" => {
                print_usage();
                std::process::exit(0);
            }
            other => {
                return Err(format!("unknown argument: {other}").into());
            }
        }
        index += 1;
    }

    Ok(Args {
        csv,
        precision_csv,
        out,
        iters,
        target,
        host_tag,
        modes,
        precisions,
        ops,
        coverage,
        strict_complete,
    })
}

fn required_value<'a>(
    args: &'a [String],
    index: usize,
    flag: &str,
) -> Result<&'a str, Box<dyn Error>> {
    args.get(index)
        .map(|value| value.as_str())
        .ok_or_else(|| format!("missing value for {flag}").into())
}

fn parse_filter(raw: &str) -> Option<Vec<String>> {
    if raw.eq_ignore_ascii_case("all") {
        return None;
    }

    let values = raw
        .split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(ToOwned::to_owned)
        .collect::<Vec<_>>();

    if values.is_empty() {
        None
    } else {
        Some(values)
    }
}

fn print_usage() {
    eprintln!(
        "usage: cargo run -p aether-fpx-research --bin perf_latex_tables -- \\\n  --csv target/perf_trig_probe.csv \\\n  --out docs/figures/perf_speed_tables.tex \\\n  [--iters 100000] \\\n  [--target 'AMD Ryzen 5 3600'] \\\n  [--host-tag x86-ryzen-3600] \\\n  [--modes nat,fpx,fpx_kernel_horner,fpx_kernel_estrin,libm,x87|all] \\\n  [--precisions f16,f32,f64,f128|all] \\\n  [--ops sin,cos,tan,asin,acos,atan,sinh,cosh,tanh,exp,expm1,log,log1p|all]"
    );
}

fn load_rows(csv_path: &Path, iters: f64) -> Result<Vec<Row>, Box<dyn Error>> {
    let text = fs::read_to_string(csv_path)?;
    parse_rows(&text, iters)
}

fn parse_rows(text: &str, iters: f64) -> Result<Vec<Row>, Box<dyn Error>> {
    let headers = text
        .lines()
        .next()
        .unwrap_or_default()
        .split(',')
        .map(str::trim)
        .collect::<Vec<_>>();
    let column = |name: &str| headers.iter().position(|header| *header == name);
    let elapsed_column =
        column("time_elapsed_s").ok_or("CSV is missing required time_elapsed_s column")?;
    let iterations_column = column("iterations");
    let energy_column = column("energy_joules");
    let power_column = column("average_power_watts");
    let power_source_column = column("power_source");
    let mut rows = Vec::new();
    for (line_index, line) in text.lines().enumerate() {
        if line_index == 0 {
            continue;
        }
        let fields = line.split(',').collect::<Vec<_>>();
        if fields.len() <= elapsed_column || fields.len() <= 5 {
            continue;
        }
        let elapsed_s = match fields[elapsed_column].parse::<f64>() {
            Ok(value) if value.is_finite() => value,
            _ => continue,
        };
        let row_iters = iterations_column
            .and_then(|index| fields.get(index))
            .and_then(|value| parse_optional_metric(value))
            .filter(|value| *value > 0.0)
            .unwrap_or(iters);
        let ins_per_call = parse_optional_metric(fields[4]).map(|value| value / row_iters);
        let cyc_per_call = parse_optional_metric(fields[5]).map(|value| value / row_iters);
        let energy_nj_per_call = energy_column
            .and_then(|index| fields.get(index))
            .and_then(|value| parse_optional_metric(value))
            .map(|joules| 1e9 * joules / row_iters);
        let average_power_watts = power_column
            .and_then(|index| fields.get(index))
            .and_then(|value| parse_optional_metric(value));
        let power_source = power_source_column
            .and_then(|index| fields.get(index))
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .map(ToOwned::to_owned);
        rows.push(Row {
            op: fields[0].to_string(),
            mode: fields[1].to_string(),
            precision: fields[2].to_string(),
            iterations: row_iters,
            ns_per_call: 1e9 * elapsed_s / row_iters,
            cyc_per_call,
            ins_per_call,
            energy_nj_per_call,
            average_power_watts,
            power_source,
        });
    }
    Ok(rows)
}

fn load_precision_rows(csv_path: &Path) -> Result<Vec<PrecisionRow>, Box<dyn Error>> {
    let text = fs::read_to_string(csv_path)?;
    parse_precision_rows(&text, csv_path)
}

fn parse_precision_rows(text: &str, csv_path: &Path) -> Result<Vec<PrecisionRow>, Box<dyn Error>> {
    let headers = text
        .lines()
        .next()
        .unwrap_or_default()
        .split(',')
        .map(str::trim)
        .collect::<Vec<_>>();
    let column = |name: &str| headers.iter().position(|header| *header == name);
    let name_column = column("name")
        .ok_or_else(|| format!("{} is missing required name column", csv_path.display()))?;
    let precision_column = column("precision").ok_or_else(|| {
        format!(
            "{} is missing required precision column",
            csv_path.display()
        )
    })?;
    let max_ulp_column = column("max_ulp")
        .ok_or_else(|| format!("{} is missing required max_ulp column", csv_path.display()))?;

    let mut rows = Vec::new();
    for (line_index, line) in text.lines().enumerate().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split(',').map(str::trim).collect::<Vec<_>>();
        let get = |column: usize, name: &str| {
            fields.get(column).copied().ok_or_else(|| {
                format!(
                    "{}:{} is missing the {name} field",
                    csv_path.display(),
                    line_index + 1
                )
            })
        };
        let op = get(name_column, "name")?;
        let precision = get(precision_column, "precision")?;
        if op.is_empty() || precision.is_empty() {
            return Err(format!(
                "{}:{} has an empty name or precision",
                csv_path.display(),
                line_index + 1
            )
            .into());
        }
        let raw_max_ulp = get(max_ulp_column, "max_ulp")?;
        let max_ulp = if raw_max_ulp.is_empty() || raw_max_ulp.eq_ignore_ascii_case("NA") {
            None
        } else {
            Some(raw_max_ulp.parse::<u128>().map_err(|error| {
                format!(
                    "{}:{} has invalid max_ulp {raw_max_ulp:?}: {error}",
                    csv_path.display(),
                    line_index + 1
                )
            })?)
        };
        rows.push(PrecisionRow {
            op: op.to_ascii_lowercase(),
            precision: precision.to_ascii_lowercase(),
            max_ulp,
        });
    }
    Ok(rows)
}

fn parse_optional_metric(raw: &str) -> Option<f64> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    match trimmed.parse::<f64>() {
        Ok(value) if value.is_finite() => Some(value),
        _ => None,
    }
}

fn render_report(
    rows: &[Row],
    precision_rows: &[PrecisionRow],
    args: &Args,
) -> Result<String, Box<dyn Error>> {
    let mut out = render_tables(rows, args)?;

    if args.coverage {
        out.push('\n');
        out.push_str(&render_coverage_table(rows, precision_rows, args));
    }
    if let Some(precision_csv) = &args.precision_csv {
        out.push('\n');
        out.push_str(&format!(
            "% Precision source: {}\n",
            precision_csv.display()
        ));
        out.push_str(&render_max_ulp_table(precision_rows, args));
    }

    Ok(out)
}

fn render_coverage_table(rows: &[Row], precision_rows: &[PrecisionRow], args: &Args) -> String {
    let target_suffix = args
        .target
        .as_deref()
        .map(|value| format!(" on {}", latex_escape(value)))
        .unwrap_or_default();
    let mut out = String::new();
    out.push_str("\\begin{table}[H]\n");
    out.push_str("\\centering\n");
    out.push_str(&format!(
        "\\caption{{Architecture result coverage by function and precision{target_suffix}}}\n"
    ));
    out.push_str(&format!(
        "\\label{{{}}}\n",
        architecture_table_label(args.host_tag.as_deref(), "result-coverage")
    ));
    out.push_str("\\resizebox{\\textwidth}{!}{%\n");
    out.push_str("\\begin{tabular}{llll}\n");
    out.push_str("\\toprule\n");
    out.push_str("Function & Speed & Precision & Power \\\\\n");
    out.push_str("\\midrule\n");
    for op in CANONICAL_OPS {
        out.push_str(&latex_escape(op));
        for kind in [
            CoverageKind::Speed,
            CoverageKind::Precision,
            CoverageKind::Power,
        ] {
            out.push_str(" & ");
            out.push_str(&format_coverage_cell(kind, op, rows, precision_rows));
        }
        out.push_str(" \\\\\n");
    }
    out.push_str("\\bottomrule\n");
    out.push_str("\\end{tabular}%\n");
    out.push_str("}\n");
    out.push_str("\\end{table}\n");
    out
}

fn render_max_ulp_table(precision_rows: &[PrecisionRow], args: &Args) -> String {
    let worst = worst_max_ulp(precision_rows);
    let target_suffix = args
        .target
        .as_deref()
        .map(|value| format!(" on {}", latex_escape(value)))
        .unwrap_or_default();
    let mut out = String::new();
    out.push_str("\\begin{table}[H]\n");
    out.push_str("\\centering\n");
    out.push_str(&format!(
        "\\caption{{Worst FPX ULP error by function and precision across all supplied tests{target_suffix}}}\n"
    ));
    out.push_str(&format!(
        "\\label{{{}}}\n",
        architecture_table_label(args.host_tag.as_deref(), "max-ulp")
    ));
    out.push_str("\\begin{tabular}{lrrrr}\n");
    out.push_str("\\toprule\n");
    out.push_str("Function & f16 & f32 & f64 & f128 \\\\\n");
    out.push_str("\\midrule\n");
    for op in CANONICAL_OPS {
        out.push_str(&latex_escape(op));
        for precision in DEFAULT_PRECISION_ORDER {
            out.push_str(" & ");
            match worst.get(&(op.to_string(), precision.to_string())) {
                Some(value) => out.push_str(&value.to_string()),
                None => out.push_str("Missing"),
            }
        }
        out.push_str(" \\\\\n");
    }
    out.push_str("\\bottomrule\n");
    out.push_str("\\end{tabular}\n");
    out.push_str("\\end{table}\n");
    out
}

fn architecture_table_label(host_tag: Option<&str>, suffix: &str) -> String {
    let host_prefix = host_tag
        .filter(|value| !value.is_empty())
        .map(|value| format!("{}-", slugify(value)))
        .unwrap_or_default();
    format!("tab:{host_prefix}architecture-{suffix}")
}

fn worst_max_ulp(precision_rows: &[PrecisionRow]) -> BTreeMap<(String, String), u128> {
    let mut worst: BTreeMap<(String, String), u128> = BTreeMap::new();
    for row in precision_rows {
        let Some(max_ulp) = row.max_ulp else {
            continue;
        };
        worst
            .entry((row.op.clone(), row.precision.clone()))
            .and_modify(|current| *current = (*current).max(max_ulp))
            .or_insert(max_ulp);
    }
    worst
}

fn format_coverage_cell(
    kind: CoverageKind,
    op: &str,
    rows: &[Row],
    precision_rows: &[PrecisionRow],
) -> String {
    let available = DEFAULT_PRECISION_ORDER
        .iter()
        .copied()
        .filter(|precision| has_coverage(kind, op, precision, rows, precision_rows))
        .collect::<Vec<_>>();
    let missing = DEFAULT_PRECISION_ORDER
        .iter()
        .copied()
        .filter(|precision| !has_coverage(kind, op, precision, rows, precision_rows))
        .collect::<Vec<_>>();

    if missing.is_empty() {
        format!("Available: {}", available.join(", "))
    } else if available.is_empty() {
        format!("Missing: {}", missing.join(", "))
    } else {
        format!(
            "Available: {}; missing: {}",
            available.join(", "),
            missing.join(", ")
        )
    }
}

fn has_coverage(
    kind: CoverageKind,
    op: &str,
    precision: &str,
    rows: &[Row],
    precision_rows: &[PrecisionRow],
) -> bool {
    match kind {
        CoverageKind::Speed => rows
            .iter()
            .any(|row| row.op == op && row.precision == precision),
        CoverageKind::Precision => precision_rows
            .iter()
            .any(|row| row.op == op && row.precision == precision && row.max_ulp.is_some()),
        CoverageKind::Power => rows.iter().any(|row| {
            row.op == op
                && row.precision == precision
                && (row.energy_nj_per_call.is_some() || row.average_power_watts.is_some())
        }),
    }
}

fn validate_complete_coverage(
    rows: &[Row],
    precision_rows: &[PrecisionRow],
) -> Result<(), Box<dyn Error>> {
    let mut gaps = Vec::new();
    for op in CANONICAL_OPS {
        for kind in [
            CoverageKind::Speed,
            CoverageKind::Precision,
            CoverageKind::Power,
        ] {
            let missing = DEFAULT_PRECISION_ORDER
                .iter()
                .copied()
                .filter(|precision| !has_coverage(kind, op, precision, rows, precision_rows))
                .collect::<Vec<_>>();
            if !missing.is_empty() {
                gaps.push(format!(
                    "{op} {} missing {}",
                    kind.label(),
                    missing.join(", ")
                ));
            }
        }
    }
    if gaps.is_empty() {
        Ok(())
    } else {
        Err(format!("incomplete architecture coverage:\n  {}", gaps.join("\n  ")).into())
    }
}

fn render_tables(rows: &[Row], args: &Args) -> Result<String, Box<dyn Error>> {
    let op_order = ordered_values(
        rows.iter().map(|row| row.op.as_str()),
        args.ops.as_deref(),
        &[],
    );
    let precision_order = ordered_values(
        rows.iter().map(|row| row.precision.as_str()),
        args.precisions.as_deref(),
        DEFAULT_PRECISION_ORDER,
    );
    let mode_order = ordered_values(
        rows.iter().map(|row| row.mode.as_str()),
        args.modes.as_deref(),
        DEFAULT_MODE_ORDER,
    );

    let mode_filter = args
        .modes
        .as_ref()
        .map(|values| values.iter().cloned().collect::<BTreeSet<_>>());
    let precision_filter = args
        .precisions
        .as_ref()
        .map(|values| values.iter().cloned().collect::<BTreeSet<_>>());
    let op_filter = args
        .ops
        .as_ref()
        .map(|values| values.iter().cloned().collect::<BTreeSet<_>>());

    let mut grouped: BTreeMap<(String, String, String, &'static str), Vec<f64>> = BTreeMap::new();
    for row in rows {
        if let Some(filter) = &mode_filter {
            if !filter.contains(&row.mode) {
                continue;
            }
        }
        if let Some(filter) = &precision_filter {
            if !filter.contains(&row.precision) {
                continue;
            }
        }
        if let Some(filter) = &op_filter {
            if !filter.contains(&row.op) {
                continue;
            }
        }
        grouped
            .entry((
                row.op.clone(),
                row.precision.clone(),
                row.mode.clone(),
                Metric::NsPerCall.key(),
            ))
            .or_default()
            .push(row.ns_per_call);
        if let Some(value) = row.cyc_per_call {
            grouped
                .entry((
                    row.op.clone(),
                    row.precision.clone(),
                    row.mode.clone(),
                    Metric::CyclesPerCall.key(),
                ))
                .or_default()
                .push(value);
        }
        if let Some(value) = row.ins_per_call {
            grouped
                .entry((
                    row.op.clone(),
                    row.precision.clone(),
                    row.mode.clone(),
                    Metric::InstructionsPerCall.key(),
                ))
                .or_default()
                .push(value);
        }
        if let Some(value) = row.energy_nj_per_call {
            grouped
                .entry((
                    row.op.clone(),
                    row.precision.clone(),
                    row.mode.clone(),
                    Metric::EnergyNjPerCall.key(),
                ))
                .or_default()
                .push(value);
        }
        if let Some(value) = row.average_power_watts {
            grouped
                .entry((
                    row.op.clone(),
                    row.precision.clone(),
                    row.mode.clone(),
                    Metric::AveragePowerWatts.key(),
                ))
                .or_default()
                .push(value);
        }
    }

    let gpu_elements = rows.iter().any(|row| {
        row.power_source
            .as_deref()
            .is_some_and(|source| source.starts_with("nvidia-smi:"))
    });
    let mut out = String::new();
    out.push_str("% Auto-generated by perf_latex_tables.\n");
    if let Some(host_tag) = &args.host_tag {
        out.push_str(&format!("% Host tag: {host_tag}\n"));
    }
    if gpu_elements {
        out.push_str("% Elapsed time and energy are normalized by the number of element evaluations in each sample.\n");
    } else if rows.iter().any(|row| row.iterations != args.iters) {
        out.push_str(&format!(
            "% Counters, elapsed time, and energy use each row's iteration count; missing values fall back to {:.0} calls.\n",
            args.iters
        ));
    } else {
        out.push_str(&format!(
            "% Counters, elapsed time, and energy are normalized by {:.0} calls per sample.\n",
            args.iters
        ));
    }
    out.push_str("% Values are mean $\\pm$ sample standard deviation.\n\n");
    let power_sources = rows
        .iter()
        .filter_map(|row| row.power_source.as_deref())
        .collect::<BTreeSet<_>>();
    if !power_sources.is_empty() {
        out.push_str(&format!(
            "% Power source(s): {}\n\n",
            power_sources.into_iter().collect::<Vec<_>>().join(", ")
        ));
    }
    let pmic_rails = rows
        .iter()
        .any(|row| row.power_source.as_deref() == Some("rpi5-pmic-rails"));

    let mut wrote_any = false;
    for op in op_order {
        let has_rows = precision_order.iter().any(|precision| {
            mode_order.iter().any(|mode| {
                grouped.contains_key(&(
                    op.clone(),
                    precision.clone(),
                    mode.clone(),
                    Metric::NsPerCall.key(),
                ))
            })
        });
        if !has_rows {
            continue;
        }
        wrote_any = true;
        out.push_str(&render_table_for_op(
            &op,
            &precision_order,
            &mode_order,
            &grouped,
            args.target.as_deref(),
            args.host_tag.as_deref(),
            pmic_rails,
            gpu_elements,
        ));
        out.push('\n');
    }

    if !wrote_any {
        return Err("no tables were generated for the requested filters".into());
    }

    Ok(out)
}

fn render_table_for_op(
    op: &str,
    precisions: &[String],
    modes: &[String],
    grouped: &BTreeMap<(String, String, String, &'static str), Vec<f64>>,
    target: Option<&str>,
    host_tag: Option<&str>,
    pmic_rails: bool,
    gpu_elements: bool,
) -> String {
    let mut metrics = vec![Metric::NsPerCall];
    for optional in [
        Metric::CyclesPerCall,
        Metric::InstructionsPerCall,
        Metric::EnergyNjPerCall,
        Metric::AveragePowerWatts,
    ] {
        if grouped
            .keys()
            .any(|(grouped_op, _, _, metric)| grouped_op == op && *metric == optional.key())
        {
            metrics.push(optional);
        }
    }
    let target_suffix = target
        .map(|value| format!(" on {}", latex_escape(value)))
        .unwrap_or_default();
    let caption = format!(
        "{} performance summary by precision and implementation{}",
        latex_escape(&title_case(op)),
        target_suffix,
    );
    let label_prefix = host_tag
        .filter(|value| !value.is_empty())
        .map(|value| format!("{}-", slugify(value)))
        .unwrap_or_default();
    let label = format!("tab:{}perf-{}-speed", label_prefix, slugify(op));
    // GPU tables are commonly widened to the full text width. Limit them to
    // three metric groups so that scaling cannot push the final rows into the
    // footer and page number.
    let max_metric_rows = if gpu_elements {
        9
    } else {
        MAX_METRIC_ROWS_PER_TABLE
    };
    let precisions_per_table = (max_metric_rows / metrics.len()).max(1);
    let part_count = precisions.len().div_ceil(precisions_per_table);
    let mut out = String::new();

    for (part_index, precision_chunk) in precisions.chunks(precisions_per_table).enumerate() {
        // A forced clear after a full-height GPU table can ship the already
        // opened next page with only its footer. Let the following [H] table
        // advance naturally instead; CPU tables retain their explicit split.
        if part_index > 0 && !gpu_elements {
            out.push_str("\\clearpage\n");
        }
        out.push_str("\\begin{table}[H]\n");
        if part_index > 0 {
            out.push_str("\\ContinuedFloat\n");
        }
        out.push_str("\\centering\n");
        if part_index == 0 {
            out.push_str(&format!("\\caption{{{caption}}}\n"));
            out.push_str(&format!("\\label{{{label}}}\n"));
        } else {
            out.push_str(&format!("\\caption[]{{{caption} (continued)}}\n"));
        }
        out.push_str("\\resizebox{\\textwidth}{!}{%\n");
        out.push_str(&format!(
            "\\begin{{tabular}}{{ll{}}}\n",
            "c".repeat(modes.len())
        ));
        out.push_str("\\toprule\n");
        out.push_str("Precision & Metric");
        for mode in modes {
            out.push_str(" & ");
            out.push_str(&mode_header(mode));
        }
        out.push_str(" \\\\\n");
        out.push_str("\\midrule\n");

        for precision in precision_chunk {
            for (metric_index, metric) in metrics.iter().enumerate() {
                if metric_index == 0 {
                    out.push_str(&latex_escape(precision));
                } else {
                    out.push(' ');
                }
                out.push_str(" & ");
                out.push_str(metric.label(pmic_rails, gpu_elements));
                for mode in modes {
                    out.push_str(" & ");
                    if let Some(values) = grouped.get(&(
                        op.to_string(),
                        precision.clone(),
                        mode.clone(),
                        metric.key(),
                    )) {
                        out.push_str(&format_stats(values, *metric, gpu_elements));
                    } else {
                        out.push_str("--");
                    }
                }
                out.push_str(" \\\\\n");
            }
            out.push_str("\\addlinespace\n");
        }

        out.push_str("\\bottomrule\n");
        out.push_str("\\end{tabular}%\n");
        out.push_str("}\n");
        out.push_str("\\end{table}\n");
        if part_index + 1 < part_count {
            out.push('\n');
        }
    }
    out
}

fn ordered_values<'a>(
    present: impl Iterator<Item = &'a str>,
    requested: Option<&[String]>,
    preferred: &[&str],
) -> Vec<String> {
    let present_set = present.map(ToOwned::to_owned).collect::<BTreeSet<_>>();
    let mut ordered = Vec::new();

    if let Some(requested) = requested {
        for value in requested {
            if present_set.contains(value) {
                ordered.push(value.clone());
            }
        }
        return ordered;
    }

    for value in preferred {
        if present_set.contains(*value) {
            ordered.push((*value).to_string());
        }
    }
    for value in present_set {
        if !ordered.iter().any(|existing| existing == &value) {
            ordered.push(value);
        }
    }
    ordered
}

fn format_stats(values: &[f64], metric: Metric, gpu_elements: bool) -> String {
    match stats(values.to_vec()) {
        Some((mean, sd, _n)) if gpu_elements && matches!(metric, Metric::NsPerCall) => {
            format!("${:.4} \\pm {:.4}$", mean, sd)
        }
        Some((mean, sd, _n)) if gpu_elements && matches!(metric, Metric::EnergyNjPerCall) => {
            format!("${:.3} \\pm {:.3}$", mean, sd)
        }
        Some((mean, sd, _n)) => format!("${:.2} \\pm {:.2}$", mean, sd),
        None => "--".to_string(),
    }
}

impl Metric {
    fn key(self) -> &'static str {
        match self {
            Metric::NsPerCall => "ns_per_call",
            Metric::CyclesPerCall => "cyc_per_call",
            Metric::InstructionsPerCall => "ins_per_call",
            Metric::EnergyNjPerCall => "energy_nj_per_call",
            Metric::AveragePowerWatts => "average_power_watts",
        }
    }

    fn label(self, pmic_rails: bool, gpu_elements: bool) -> &'static str {
        match self {
            Metric::NsPerCall if gpu_elements => "ns/element",
            Metric::NsPerCall => "ns/call",
            Metric::CyclesPerCall => "cyc/call",
            Metric::InstructionsPerCall => "ins/call",
            Metric::EnergyNjPerCall if pmic_rails => "PMIC rail nJ/call",
            Metric::AveragePowerWatts if pmic_rails => "PMIC rail mean W",
            Metric::EnergyNjPerCall if gpu_elements => "nJ/element",
            Metric::EnergyNjPerCall => "nJ/call",
            Metric::AveragePowerWatts => "mean W",
        }
    }
}

impl CoverageKind {
    fn label(self) -> &'static str {
        match self {
            CoverageKind::Speed => "speed",
            CoverageKind::Precision => "precision",
            CoverageKind::Power => "power",
        }
    }
}

fn mode_header(mode: &str) -> String {
    match mode {
        "nat" => "native".to_string(),
        "fpx" => "FPX".to_string(),
        "fpx_blocked_horner" => "FPX blocked Horner".to_string(),
        "fpx_horizontal_estrin" => "FPX horizontal Estrin".to_string(),
        "fpx_kernel_horner" => "Horner-fpx".to_string(),
        "nat_kernel" => "native kernel".to_string(),
        "fpx_kernel_estrin" => "Estrin-fpx".to_string(),
        "fpx_kernel_scalar" => "Scalar".to_string(),
        "fpx_kernel_fma" => "Scalar FMA".to_string(),
        "fpx_kernel_simd" => "Packed SIMD FMA".to_string(),
        "fpx_kernel_simd_estrin" => "Packed SIMD Estrin".to_string(),
        "fpx_kernel_raw_scalar" => "Raw scalar".to_string(),
        "fpx_kernel_raw_fma" => "Raw scalar FMA".to_string(),
        "fpx_kernel_raw_simd" => "Raw packed SIMD FMA".to_string(),
        "fpx_kernel_raw_simd_estrin" => "Raw packed SIMD Estrin".to_string(),
        "libm" => "libm".to_string(),
        "x87" => "x87".to_string(),
        other => latex_escape(other),
    }
}

fn title_case(value: &str) -> String {
    let mut chars = value.chars();
    match chars.next() {
        Some(first) => first.to_ascii_uppercase().to_string() + chars.as_str(),
        None => String::new(),
    }
}

fn latex_escape(value: &str) -> String {
    value
        .replace('\\', "\\textbackslash{}")
        .replace('&', "\\&")
        .replace('%', "\\%")
        .replace('_', "\\_")
        .replace('#', "\\#")
        .replace('{', "\\{")
        .replace('}', "\\}")
}

fn slugify(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    let mut last_was_sep = false;
    for ch in value.chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch.to_ascii_lowercase());
            last_was_sep = false;
        } else if !last_was_sep {
            out.push('-');
            last_was_sep = true;
        }
    }
    out.trim_matches('-').to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kernel_headers_match_plot_legend_labels() {
        assert_eq!(mode_header("fpx_kernel_horner"), "Horner-fpx");
        assert_eq!(mode_header("fpx_blocked_horner"), "FPX blocked Horner");
        assert_eq!(
            mode_header("fpx_horizontal_estrin"),
            "FPX horizontal Estrin"
        );
        assert_eq!(mode_header("fpx_kernel_estrin"), "Estrin-fpx");
        assert_eq!(mode_header("fpx_kernel_scalar"), "Scalar");
        assert_eq!(mode_header("fpx_kernel_fma"), "Scalar FMA");
        assert_eq!(mode_header("fpx_kernel_simd"), "Packed SIMD FMA");
        assert_eq!(mode_header("fpx_kernel_raw_scalar"), "Raw scalar");
        assert_eq!(mode_header("fpx_kernel_raw_fma"), "Raw scalar FMA");
        assert_eq!(mode_header("fpx_kernel_raw_simd"), "Raw packed SIMD FMA");
    }

    #[test]
    fn host_tag_is_slugified() {
        let args = parse_args(vec!["--host-tag".into(), "Raspberry Pi 5 / Lab".into()]).unwrap();
        assert_eq!(args.iters, 100_000.0);
        assert_eq!(args.host_tag.as_deref(), Some("raspberry-pi-5-lab"));
        assert!(args.precision_csv.is_none());
        assert!(!args.coverage);
        assert!(!args.strict_complete);
    }

    #[test]
    fn architecture_reporting_flags_are_parsed() {
        let args = parse_args(vec![
            "--precision-csv".into(),
            "target/precision_summary.csv".into(),
            "--coverage".into(),
            "--strict-complete".into(),
        ])
        .unwrap();

        assert_eq!(
            args.precision_csv.as_deref(),
            Some(Path::new("target/precision_summary.csv"))
        );
        assert!(args.coverage);
        assert!(args.strict_complete);
    }

    #[test]
    fn legacy_label_prefix_remains_an_alias() {
        let args = parse_args(vec!["--label-prefix".into(), "Ryzen 3600".into()]).unwrap();
        assert_eq!(args.host_tag.as_deref(), Some("ryzen-3600"));
    }

    #[test]
    fn host_tag_qualifies_generated_table_labels() {
        let args = Args {
            csv: PathBuf::from("unused.csv"),
            precision_csv: None,
            out: PathBuf::from("unused.tex"),
            iters: 1.0,
            target: Some("Raspberry Pi 5".into()),
            host_tag: Some("rpi5-aarch64".into()),
            modes: None,
            precisions: Some(vec!["f64".into()]),
            ops: Some(vec!["sin".into()]),
            coverage: false,
            strict_complete: false,
        };
        let rows = [Row {
            op: "sin".into(),
            mode: "fpx".into(),
            precision: "f64".into(),
            iterations: 1.0,
            ns_per_call: 1.0,
            cyc_per_call: Some(2.0),
            ins_per_call: Some(3.0),
            energy_nj_per_call: Some(4.0),
            average_power_watts: Some(5.0),
            power_source: Some("rpi5-pmic-rails".into()),
        }];

        let latex = render_tables(&rows, &args).unwrap();
        assert!(latex.contains("% Host tag: rpi5-aarch64"));
        assert!(latex.contains("\\label{tab:rpi5-aarch64-perf-sin-speed}"));
        assert!(latex.contains("on Raspberry Pi 5"));
        assert!(latex.contains("PMIC rail nJ/call"));
        assert!(latex.contains("PMIC rail mean W"));
        assert!(latex.contains("% Power source(s): rpi5-pmic-rails"));
    }

    #[test]
    fn power_tables_split_before_exceeding_a_page() {
        let args = Args {
            csv: PathBuf::from("unused.csv"),
            precision_csv: None,
            out: PathBuf::from("unused.tex"),
            iters: 1.0,
            target: Some("Raspberry Pi 5".into()),
            host_tag: Some("rpi5-aarch64".into()),
            modes: Some(vec!["fpx".into()]),
            precisions: Some(
                DEFAULT_PRECISION_ORDER
                    .iter()
                    .map(|precision| (*precision).to_string())
                    .collect(),
            ),
            ops: Some(vec!["sin".into()]),
            coverage: false,
            strict_complete: false,
        };
        let rows = DEFAULT_PRECISION_ORDER
            .iter()
            .map(|precision| Row {
                op: "sin".into(),
                mode: "fpx".into(),
                precision: (*precision).into(),
                iterations: 1.0,
                ns_per_call: 1.0,
                cyc_per_call: Some(2.0),
                ins_per_call: Some(3.0),
                energy_nj_per_call: Some(4.0),
                average_power_watts: Some(5.0),
                power_source: Some("rpi5-pmic-rails".into()),
            })
            .collect::<Vec<_>>();

        let latex = render_tables(&rows, &args).unwrap();
        assert_eq!(latex.matches("\\begin{table}[H]").count(), 2);
        assert_eq!(latex.matches("\\ContinuedFloat").count(), 1);
        assert_eq!(latex.matches("\\clearpage").count(), 1);
        assert_eq!(
            latex
                .matches("\\label{tab:rpi5-aarch64-perf-sin-speed}")
                .count(),
            1
        );
        assert!(latex.contains("(continued)"));
    }

    #[test]
    fn precision_csv_preserves_zero_and_marks_na_missing() {
        let rows = parse_precision_rows(
            "name,precision,test,max_ulp\nsin,f32,a,0\nsin,f64,b,NA\n",
            Path::new("precision.csv"),
        )
        .unwrap();

        assert_eq!(
            rows,
            vec![
                PrecisionRow {
                    op: "sin".into(),
                    precision: "f32".into(),
                    max_ulp: Some(0),
                },
                PrecisionRow {
                    op: "sin".into(),
                    precision: "f64".into(),
                    max_ulp: None,
                },
            ]
        );
    }

    #[test]
    fn row_iterations_override_legacy_normalization_fallback() {
        let rows = parse_rows(
            "op,mode,precision,run,instructions,cycles,branches,branch_misses,time_elapsed_s,iterations,energy_joules\n\
             sin,fpx,f64,0,20,12,0,0,0.000000008,4,0.000000008\n\
             cos,fpx,f64,0,20,10,0,0,0.000000010,,0.000000010\n",
            10.0,
        )
        .unwrap();

        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].iterations, 4.0);
        assert!((rows[0].ns_per_call - 2.0).abs() < 1e-12);
        assert_eq!(rows[0].ins_per_call, Some(5.0));
        assert_eq!(rows[0].cyc_per_call, Some(3.0));
        assert!((rows[0].energy_nj_per_call.unwrap() - 2.0).abs() < 1e-12);

        assert_eq!(rows[1].iterations, 10.0);
        assert!((rows[1].ns_per_call - 1.0).abs() < 1e-12);
        assert_eq!(rows[1].ins_per_call, Some(2.0));
        assert_eq!(rows[1].cyc_per_call, Some(1.0));
        assert!((rows[1].energy_nj_per_call.unwrap() - 1.0).abs() < 1e-12);

        let mut args = reporting_args();
        args.iters = 10.0;
        let latex = render_tables(&rows, &args).unwrap();
        assert!(latex.contains("use each row's iteration count"));
    }

    #[test]
    fn max_ulp_table_uses_worst_test_and_explicit_missing_values() {
        let args = reporting_args();
        let precision_rows = vec![
            PrecisionRow {
                op: "sin".into(),
                precision: "f64".into(),
                max_ulp: Some(2),
            },
            PrecisionRow {
                op: "sin".into(),
                precision: "f64".into(),
                max_ulp: Some(9),
            },
            PrecisionRow {
                op: "sin".into(),
                precision: "f32".into(),
                max_ulp: None,
            },
        ];

        let latex = render_max_ulp_table(&precision_rows, &args);
        assert!(latex.contains("sin & Missing & Missing & 9 & Missing"));
        assert!(latex.contains("\\label{tab:test-host-architecture-max-ulp}"));
        assert_eq!(
            latex
                .lines()
                .filter(|line| CANONICAL_OPS.iter().any(|op| line.starts_with(op)))
                .count(),
            CANONICAL_OPS.len()
        );
    }

    #[test]
    fn coverage_table_has_one_row_per_canonical_function() {
        let args = reporting_args();
        let latex = render_coverage_table(&[], &[], &args);

        for op in CANONICAL_OPS {
            assert!(
                latex
                    .lines()
                    .any(|line| line.starts_with(&format!("{op} &"))),
                "missing coverage row for {op}"
            );
        }
        assert_eq!(latex.matches("Missing: f16, f32, f64, f128").count(), 39);
        assert!(latex.contains("\\label{tab:test-host-architecture-result-coverage}"));
    }

    #[test]
    fn strict_coverage_requires_every_format_for_every_result_kind() {
        let (mut rows, mut precision_rows) = complete_coverage_rows();
        validate_complete_coverage(&rows, &precision_rows).unwrap();

        rows.retain(|row| !(row.op == "sin" && row.precision == "f16"));
        let cos_f32 = rows
            .iter_mut()
            .find(|row| row.op == "cos" && row.precision == "f32")
            .unwrap();
        cos_f32.energy_nj_per_call = None;
        cos_f32.average_power_watts = None;
        precision_rows
            .iter_mut()
            .find(|row| row.op == "log1p" && row.precision == "f128")
            .unwrap()
            .max_ulp = None;

        let error = validate_complete_coverage(&rows, &precision_rows)
            .unwrap_err()
            .to_string();
        assert!(error.contains("sin speed missing f16"));
        assert!(error.contains("sin power missing f16"));
        assert!(error.contains("cos power missing f32"));
        assert!(error.contains("log1p precision missing f128"));
    }

    #[test]
    fn default_report_is_byte_identical_to_existing_tables() {
        let args = reporting_args();
        let rows = vec![coverage_row("sin", "f64")];

        assert_eq!(
            render_report(&rows, &[], &args).unwrap(),
            render_tables(&rows, &args).unwrap()
        );
    }

    fn reporting_args() -> Args {
        Args {
            csv: PathBuf::from("unused.csv"),
            precision_csv: None,
            out: PathBuf::from("unused.tex"),
            iters: 1.0,
            target: Some("Test Architecture".into()),
            host_tag: Some("test-host".into()),
            modes: None,
            precisions: None,
            ops: None,
            coverage: false,
            strict_complete: false,
        }
    }

    fn coverage_row(op: &str, precision: &str) -> Row {
        Row {
            op: op.into(),
            mode: "fpx".into(),
            precision: precision.into(),
            iterations: 1.0,
            ns_per_call: 1.0,
            cyc_per_call: None,
            ins_per_call: None,
            energy_nj_per_call: Some(2.0),
            average_power_watts: Some(3.0),
            power_source: Some("test".into()),
        }
    }

    fn complete_coverage_rows() -> (Vec<Row>, Vec<PrecisionRow>) {
        let rows = CANONICAL_OPS
            .iter()
            .flat_map(|op| {
                DEFAULT_PRECISION_ORDER
                    .iter()
                    .map(move |precision| coverage_row(op, precision))
            })
            .collect();
        let precision_rows = CANONICAL_OPS
            .iter()
            .flat_map(|op| {
                DEFAULT_PRECISION_ORDER
                    .iter()
                    .map(move |precision| PrecisionRow {
                        op: (*op).into(),
                        precision: (*precision).into(),
                        max_ulp: Some(1),
                    })
            })
            .collect();
        (rows, precision_rows)
    }
}
