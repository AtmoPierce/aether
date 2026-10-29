use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use aether_viz::{plot_histograms_with_config, HistogramConfig, HistogramSeries};

const DEFAULT_CSV: &str = "target/perf_trig_probe.csv";
const DEFAULT_OUT: &str = "docs/figures/perf_speed";
const DEFAULT_ITERS: f64 = 10_000_000.0;
const DEFAULT_MODE_ORDER: &[&str] = &[
    "fpx",
    "fpx_kernel_horner",
    "fpx_kernel_estrin",
    "libm",
    "nat",
    "x87",
];
const DEFAULT_PRECISION_ORDER: &[&str] = &["f16", "f32", "f64", "f128"];

#[derive(Debug, Clone)]
struct Row {
    op: String,
    mode: String,
    precision: String,
    ns_per_call: f64,
    energy_nj_per_call: Option<f64>,
    average_power_watts: Option<f64>,
    power_source: Option<String>,
}

#[derive(Clone, Copy, Debug)]
enum Metric {
    Speed,
    Energy,
    Power,
}

#[derive(Debug, Clone)]
struct Args {
    csv: PathBuf,
    out: PathBuf,
    iters: f64,
    target: Option<String>,
    modes: Option<Vec<String>>,
    precisions: Option<Vec<String>>,
    ops: Option<Vec<String>>,
    bins: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = parse_args(env::args().skip(1).collect())?;
    let rows = load_rows(&args.csv, args.iters)?;
    if rows.is_empty() {
        return Err(format!("no usable rows found in {}", args.csv.display()).into());
    }

    fs::create_dir_all(&args.out)?;
    let written = write_plots(&rows, &args)?;
    if written.is_empty() {
        return Err("no plots were generated for the requested filters".into());
    }

    for path in written {
        println!("Wrote {}", path.display());
    }
    Ok(())
}

fn parse_args(args: Vec<String>) -> Result<Args, Box<dyn Error>> {
    let mut csv = PathBuf::from(DEFAULT_CSV);
    let mut out = PathBuf::from(DEFAULT_OUT);
    let mut iters = DEFAULT_ITERS;
    let mut target = None;
    let mut modes = None;
    let mut precisions = Some(
        DEFAULT_PRECISION_ORDER
            .iter()
            .map(|value| (*value).to_string())
            .collect::<Vec<_>>(),
    );
    let mut ops = None;
    let mut bins = 12usize;

    let mut index = 0usize;
    while index < args.len() {
        match args[index].as_str() {
            "--csv" => {
                index += 1;
                csv = PathBuf::from(required_value(&args, index, "--csv")?);
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
            "--bins" => {
                index += 1;
                bins = required_value(&args, index, "--bins")?.parse()?;
            }
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
        out,
        iters,
        target,
        modes,
        precisions,
        ops,
        bins: bins.max(1),
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
        "usage: cargo run -p aether-fpx-research --bin perf_speed_plots -- \\\n  --csv target/perf_trig_probe.csv \\\n  --out docs/figures/perf_speed \\\n  [--iters 10000000] \\\n  [--target 'AMD Ryzen 5 3600'] \\\n  [--modes fpx,fpx_kernel_horner,fpx_kernel_estrin,libm,nat,x87|all] \\\n  [--precisions f32,f64|all] \\\n  [--ops sin,cos,tan|all] \\\n  [--bins 18]"
    );
}

fn load_rows(csv_path: &Path, iters: f64) -> Result<Vec<Row>, Box<dyn Error>> {
    let text = fs::read_to_string(csv_path)?;
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
    let energy_column = column("energy_joules");
    let power_column = column("average_power_watts");
    let source_column = column("power_source");
    let iterations_column = column("iterations");
    let mut rows = Vec::new();
    for (line_index, line) in text.lines().enumerate() {
        if line_index == 0 {
            continue;
        }
        let fields = line.split(',').collect::<Vec<_>>();
        if fields.len() <= elapsed_column || fields.len() < 4 {
            continue;
        }
        let elapsed_s = match fields[elapsed_column].parse::<f64>() {
            Ok(value) if value.is_finite() => value,
            _ => continue,
        };
        let row_iterations = iterations_column
            .and_then(|index| fields.get(index))
            .and_then(|value| value.trim().parse::<f64>().ok())
            .filter(|value| value.is_finite() && *value > 0.0)
            .unwrap_or(iters);
        rows.push(Row {
            op: fields[0].to_string(),
            mode: fields[1].to_string(),
            precision: fields[2].to_string(),
            ns_per_call: 1e9 * elapsed_s / row_iterations,
            energy_nj_per_call: energy_column
                .and_then(|index| fields.get(index))
                .and_then(|value| parse_optional_metric(value))
                .map(|joules| 1e9 * joules / row_iterations),
            average_power_watts: power_column
                .and_then(|index| fields.get(index))
                .and_then(|value| parse_optional_metric(value)),
            power_source: source_column
                .and_then(|index| fields.get(index))
                .map(|value| value.trim())
                .filter(|value| !value.is_empty())
                .map(ToOwned::to_owned),
        });
    }
    Ok(rows)
}

fn parse_optional_metric(raw: &str) -> Option<f64> {
    match raw.trim().parse::<f64>() {
        Ok(value) if value.is_finite() => Some(value),
        _ => None,
    }
}

fn write_plots(rows: &[Row], args: &Args) -> Result<Vec<PathBuf>, Box<dyn Error>> {
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
                Metric::Speed.key(),
            ))
            .or_default()
            .push(row.ns_per_call);
        if let Some(value) = row.energy_nj_per_call {
            grouped
                .entry((
                    row.op.clone(),
                    row.precision.clone(),
                    row.mode.clone(),
                    Metric::Energy.key(),
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
                    Metric::Power.key(),
                ))
                .or_default()
                .push(value);
        }
    }

    let pmic_rails = rows
        .iter()
        .any(|row| row.power_source.as_deref() == Some("rpi5-pmic-rails"));
    let gpu_elements = rows.iter().any(|row| {
        row.power_source
            .as_deref()
            .is_some_and(|source| source.starts_with("nvidia-smi:"))
    });

    let mut written = Vec::new();
    for op in op_order {
        for precision in &precision_order {
            for metric in [Metric::Speed, Metric::Energy, Metric::Power] {
                let values_by_label = mode_order
                    .iter()
                    .filter_map(|mode| {
                        grouped
                            .get(&(op.clone(), precision.clone(), mode.clone(), metric.key()))
                            .filter(|values| !values.is_empty())
                            .map(|values| {
                                (mode_label(mode), mode_palette_index(mode), values.clone())
                            })
                    })
                    .collect::<Vec<_>>();
                if values_by_label.is_empty() {
                    continue;
                }

                let series = values_by_label
                    .iter()
                    .map(|(label, palette_index, values)| {
                        HistogramSeries::new(values)
                            .with_label(label)
                            .with_bins(args.bins)
                            .with_palette(*palette_index)
                    })
                    .collect::<Vec<_>>();
                let path = args.out.join(format!(
                    "{}_{}_{}_hist.svg",
                    slugify(&op),
                    precision,
                    metric.file_component()
                ));
                let title = metric.title(
                    &op,
                    precision,
                    args.target.as_deref(),
                    pmic_rails,
                    gpu_elements,
                );
                let path_string = path.to_string_lossy().into_owned();
                let config = HistogramConfig::publication(&title)
                    .with_x_label(metric.axis_label(pmic_rails, gpu_elements))
                    .with_y_label("samples")
                    .with_save_path(&path_string);
                plot_histograms_with_config(&series, &config)?;
                written.push(path);
            }
        }
    }

    Ok(written)
}

impl Metric {
    fn key(self) -> &'static str {
        match self {
            Self::Speed => "speed",
            Self::Energy => "energy",
            Self::Power => "power",
        }
    }

    fn file_component(self) -> &'static str {
        match self {
            Self::Speed => "speed",
            Self::Energy => "energy",
            Self::Power => "power",
        }
    }

    fn axis_label(self, pmic_rails: bool, gpu_elements: bool) -> &'static str {
        match (self, pmic_rails, gpu_elements) {
            (Self::Speed, _, true) => "ns / element",
            (Self::Speed, _, false) => "ns / call",
            (Self::Energy, true, _) => "PMIC rail nJ / call",
            (Self::Power, true, _) => "PMIC rail mean W",
            (Self::Energy, false, true) => "nJ / element",
            (Self::Energy, false, false) => "nJ / call",
            (Self::Power, false, _) => "mean W",
        }
    }

    fn title(
        self,
        op: &str,
        precision: &str,
        target: Option<&str>,
        pmic_rails: bool,
        _gpu_elements: bool,
    ) -> String {
        let subject = match (self, pmic_rails) {
            (Self::Speed, _) => "speed",
            (Self::Energy, true) => "PMIC rail energy",
            (Self::Power, true) => "PMIC rail power",
            (Self::Energy, false) => "energy",
            (Self::Power, false) => "power",
        };
        match target {
            Some(target) => format!(
                "{op} {precision} {subject} — {}",
                compact_target_label(target)
            ),
            None => format!("{op} {precision} {subject} histogram"),
        }
    }
}

fn compact_target_label(target: &str) -> &str {
    match target {
        "Raspberry Pi 5, Cortex-A76 (AArch64)" => "RPi 5 AArch64",
        _ => target,
    }
}

fn mode_label(mode: &str) -> String {
    match mode {
        "fpx" => "FPX".to_string(),
        "fpx_kernel_horner" => "Horner-fpx".to_string(),
        "nat_kernel" => "native kernel".to_string(),
        "fpx_kernel_estrin" => "Estrin-fpx".to_string(),
        "libm" => "libm".to_string(),
        "nat" => "native".to_string(),
        "x87" => "x87".to_string(),
        other => other.to_string(),
    }
}

fn mode_palette_index(mode: &str) -> usize {
    match mode {
        "fpx" => 0,
        "fpx_kernel_horner" => 1,
        "fpx_kernel_estrin" => 2,
        "libm" => 3,
        "nat" => 4,
        "x87" => 5,
        other => {
            let hash = other.bytes().fold(2_166_136_261_u32, |hash, byte| {
                (hash ^ u32::from(byte)).wrapping_mul(16_777_619)
            });
            6 + hash as usize % 94
        }
    }
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

fn slugify(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    let mut last_was_sep = false;
    for ch in value.chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch.to_ascii_lowercase());
            last_was_sep = false;
        } else if !last_was_sep {
            out.push('_');
            last_was_sep = true;
        }
    }
    out.trim_matches('_').to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_iterations_override_the_legacy_fallback() {
        let path = env::temp_dir().join(format!(
            "fpx-perf-iterations-{}-{}.csv",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::write(
            &path,
            "operation,mode,precision,run,instructions,cycles,time_elapsed_s,energy_joules,average_power_watts,power_source,iterations\n\
             sin,fpx,f64,0,20,30,1,1,2,fixture,2\n",
        )
        .unwrap();

        let rows = load_rows(&path, 10.0).unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].ns_per_call, 500_000_000.0);
        assert_eq!(rows[0].energy_nj_per_call, Some(500_000_000.0));
    }

    #[test]
    fn writes_speed_energy_and_power_plots_for_pmic_rows() {
        let out = env::temp_dir().join(format!(
            "fpx-power-plots-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&out).unwrap();
        let args = Args {
            csv: PathBuf::from("unused.csv"),
            out: out.clone(),
            iters: 10_000_000.0,
            target: Some("Raspberry Pi 5".into()),
            modes: None,
            precisions: Some(vec!["f64".into()]),
            ops: Some(vec!["sin".into()]),
            bins: 4,
        };
        let rows = vec![
            Row {
                op: "sin".into(),
                mode: "fpx".into(),
                precision: "f64".into(),
                ns_per_call: 10.0,
                energy_nj_per_call: Some(20.0),
                average_power_watts: Some(2.0),
                power_source: Some("rpi5-pmic-rails".into()),
            },
            Row {
                op: "sin".into(),
                mode: "libm".into(),
                precision: "f64".into(),
                ns_per_call: 12.0,
                energy_nj_per_call: Some(25.0),
                average_power_watts: Some(2.1),
                power_source: Some("rpi5-pmic-rails".into()),
            },
        ];

        let written = write_plots(&rows, &args).unwrap();
        assert_eq!(written.len(), 3);
        assert!(out.join("sin_f64_speed_hist.svg").exists());
        assert!(out.join("sin_f64_energy_hist.svg").exists());
        assert!(out.join("sin_f64_power_hist.svg").exists());
        let svg = fs::read_to_string(out.join("sin_f64_energy_hist.svg")).unwrap();
        assert!(svg.contains("<svg width=\"1200\" height=\"675\""));

        fs::remove_dir_all(out).unwrap();
    }

    #[test]
    fn compacts_the_aarch64_plot_target() {
        assert_eq!(
            compact_target_label("Raspberry Pi 5, Cortex-A76 (AArch64)"),
            "RPi 5 AArch64"
        );
    }

    #[test]
    fn uses_compact_kernel_legend_labels() {
        assert_eq!(mode_label("fpx_kernel_horner"), "Horner-fpx");
        assert_eq!(mode_label("fpx_kernel_estrin"), "Estrin-fpx");
    }
}
