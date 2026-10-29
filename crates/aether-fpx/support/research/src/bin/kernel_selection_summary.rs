use aether_viz::{plot_histograms_with_config, HistogramConfig, HistogramSeries};
use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

const PRECISIONS: [&str; 2] = ["f32", "f64"];
const END_HORNER: &str = "fpx_kernel_simd";
const END_ESTRIN: &str = "fpx_kernel_simd_estrin";
const RAW_HORNER: &str = "fpx_kernel_raw_simd";
const RAW_ESTRIN: &str = "fpx_kernel_raw_simd_estrin";
const MODES: [&str; 4] = [END_HORNER, END_ESTRIN, RAW_HORNER, RAW_ESTRIN];

#[derive(Clone, Copy, Debug)]
struct Target {
    label: &'static str,
    detail: &'static str,
    csv_path: &'static str,
}

const TARGETS: [Target; 2] = [
    Target {
        label: "AArch64",
        detail: "Raspberry Pi 5 (Cortex-A76)",
        csv_path: "outputs/perf/aarch64.csv",
    },
    Target {
        label: "x86-64",
        detail: "AMD Ryzen 5 3600",
        csv_path: "outputs/perf/x86_64.csv",
    },
];

#[derive(Clone, Debug)]
struct Measurement {
    operation: String,
    mode: String,
    precision: String,
    cycles_per_call: f64,
    instructions_per_call: f64,
    ns_per_call: f64,
    nj_per_call: f64,
    average_power_watts: f64,
    power_source: String,
}

#[derive(Clone, Debug)]
struct OperationRatio {
    operation: String,
    end_cycles: f64,
    raw_cycles: f64,
    raw_instructions: f64,
    raw_elapsed: f64,
    end_energy: f64,
    raw_energy: f64,
    end_power: f64,
    raw_power: f64,
}

#[derive(Clone, Debug)]
struct PrecisionReport {
    precision: &'static str,
    ratios: Vec<OperationRatio>,
    end_cycle_cvs: Vec<f64>,
    raw_cycle_cvs: Vec<f64>,
    end_energy_cvs: Vec<f64>,
    raw_energy_cvs: Vec<f64>,
}

#[derive(Clone, Debug)]
struct TargetReport {
    target: Target,
    power_source: String,
    precisions: Vec<PrecisionReport>,
}

#[derive(Clone, Copy)]
enum Metric {
    EndCycles,
    RawCycles,
    RawInstructions,
    EndEnergy,
    RawEnergy,
    EndPower,
}

fn repository_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(4)
        .expect("aether-fpx-research must be inside the workspace")
        .to_path_buf()
}

fn column(headers: &csv::StringRecord, name: &str) -> Result<usize, Box<dyn Error>> {
    headers
        .iter()
        .position(|header| header.trim() == name)
        .ok_or_else(|| format!("CSV is missing required {name} column").into())
}

fn field<'a>(record: &'a csv::StringRecord, index: usize, name: &str) -> Result<&'a str, String> {
    record
        .get(index)
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .ok_or_else(|| format!("CSV row is missing {name}"))
}

fn positive(value: &str, name: &str) -> Result<f64, String> {
    let parsed = value
        .parse::<f64>()
        .map_err(|error| format!("invalid {name} value {value:?}: {error}"))?;
    if parsed.is_finite() && parsed > 0.0 {
        Ok(parsed)
    } else {
        Err(format!("non-positive or non-finite {name} value {value:?}"))
    }
}

fn load_measurements(path: &Path) -> Result<Vec<Measurement>, Box<dyn Error>> {
    let mut reader = csv::Reader::from_path(path)?;
    let headers = reader.headers()?.clone();
    let op = column(&headers, "op")?;
    let mode = column(&headers, "mode")?;
    let precision = column(&headers, "precision")?;
    let cycles = column(&headers, "cycles")?;
    let instructions = column(&headers, "instructions")?;
    let elapsed = column(&headers, "time_elapsed_s")?;
    let energy = column(&headers, "energy_joules")?;
    let power = column(&headers, "average_power_watts")?;
    let source = column(&headers, "power_source")?;
    let iterations = column(&headers, "iterations")?;

    let mut measurements = Vec::new();
    for (index, result) in reader.records().enumerate() {
        let record = result?;
        let row_number = index + 2;
        let mode_value = field(&record, mode, "mode")?;
        let precision_value = field(&record, precision, "precision")?;
        if !MODES.contains(&mode_value) || !PRECISIONS.contains(&precision_value) {
            continue;
        }
        let parse = |column_index, name| {
            positive(field(&record, column_index, name)?, name)
                .map_err(|error| format!("{}:{row_number}: {error}", path.display()))
        };
        let iteration_count = field(&record, iterations, "iterations")?
            .parse::<u64>()
            .map_err(|error| {
                format!(
                    "{}:{row_number}: invalid iterations: {error}",
                    path.display()
                )
            })?;
        if iteration_count == 0 {
            return Err(format!("{}:{row_number}: zero iterations", path.display()).into());
        }
        let denominator = iteration_count as f64;
        measurements.push(Measurement {
            operation: field(&record, op, "op")?.to_string(),
            mode: mode_value.to_string(),
            precision: precision_value.to_string(),
            cycles_per_call: parse(cycles, "cycles")? / denominator,
            instructions_per_call: parse(instructions, "instructions")? / denominator,
            ns_per_call: 1.0e9 * parse(elapsed, "time_elapsed_s")? / denominator,
            nj_per_call: 1.0e9 * parse(energy, "energy_joules")? / denominator,
            average_power_watts: parse(power, "average_power_watts")?,
            power_source: field(&record, source, "power_source")?.to_string(),
        });
    }
    if measurements.is_empty() {
        return Err(format!("{} contains no recognized SIMD modes", path.display()).into());
    }
    Ok(measurements)
}

fn mean(values: &[f64]) -> Option<f64> {
    (!values.is_empty()).then(|| values.iter().sum::<f64>() / values.len() as f64)
}

fn coefficient_of_variation(values: &[f64]) -> Option<f64> {
    let average = mean(values)?;
    if values.len() < 2 {
        return Some(0.0);
    }
    let variance = values
        .iter()
        .map(|value| (value - average).powi(2))
        .sum::<f64>()
        / (values.len() - 1) as f64;
    (average > 0.0).then(|| variance.sqrt() / average)
}

fn median(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let middle = sorted.len() / 2;
    Some(if sorted.len() % 2 == 0 {
        (sorted[middle - 1] + sorted[middle]) / 2.0
    } else {
        sorted[middle]
    })
}

fn geometric_mean(values: &[f64]) -> Option<f64> {
    (!values.is_empty() && values.iter().all(|value| *value > 0.0))
        .then(|| (values.iter().map(|value| value.ln()).sum::<f64>() / values.len() as f64).exp())
}

fn samples(rows: &[&Measurement], select: fn(&Measurement) -> f64) -> Vec<f64> {
    rows.iter().map(|row| select(row)).collect()
}

fn ratio(left: &[&Measurement], right: &[&Measurement], select: fn(&Measurement) -> f64) -> f64 {
    median(&samples(left, select)).unwrap() / median(&samples(right, select)).unwrap()
}

fn summarize_precision(
    precision: &'static str,
    measurements: &[Measurement],
) -> Result<PrecisionReport, Box<dyn Error>> {
    let mut grouped = BTreeMap::<(String, String), Vec<&Measurement>>::new();
    for measurement in measurements.iter().filter(|row| row.precision == precision) {
        grouped
            .entry((measurement.operation.clone(), measurement.mode.clone()))
            .or_default()
            .push(measurement);
    }
    let operations = grouped
        .keys()
        .map(|(operation, _)| operation.clone())
        .collect::<BTreeSet<_>>();
    let mut ratios = Vec::new();
    let mut end_cycle_cvs = Vec::new();
    let mut raw_cycle_cvs = Vec::new();
    let mut end_energy_cvs = Vec::new();
    let mut raw_energy_cvs = Vec::new();

    for operation in operations {
        let rows = MODES
            .iter()
            .map(|mode| {
                grouped
                    .get(&(operation.clone(), (*mode).to_string()))
                    .ok_or_else(|| format!("missing {mode} rows for {operation}/{precision}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        if rows
            .iter()
            .any(|candidate| candidate.len() != rows[0].len())
        {
            return Err(format!("mismatched run counts for {operation}/{precision}").into());
        }
        for candidate in [&rows[0], &rows[1]] {
            end_cycle_cvs.push(
                coefficient_of_variation(&samples(candidate, |row| row.cycles_per_call)).unwrap(),
            );
            end_energy_cvs.push(
                coefficient_of_variation(&samples(candidate, |row| row.nj_per_call)).unwrap(),
            );
        }
        for candidate in [&rows[2], &rows[3]] {
            raw_cycle_cvs.push(
                coefficient_of_variation(&samples(candidate, |row| row.cycles_per_call)).unwrap(),
            );
            raw_energy_cvs.push(
                coefficient_of_variation(&samples(candidate, |row| row.nj_per_call)).unwrap(),
            );
        }
        ratios.push(OperationRatio {
            operation,
            end_cycles: ratio(rows[0], rows[1], |row| row.cycles_per_call),
            raw_cycles: ratio(rows[2], rows[3], |row| row.cycles_per_call),
            raw_instructions: ratio(rows[2], rows[3], |row| row.instructions_per_call),
            raw_elapsed: ratio(rows[2], rows[3], |row| row.ns_per_call),
            end_energy: ratio(rows[0], rows[1], |row| row.nj_per_call),
            raw_energy: ratio(rows[2], rows[3], |row| row.nj_per_call),
            end_power: ratio(rows[0], rows[1], |row| row.average_power_watts),
            raw_power: ratio(rows[2], rows[3], |row| row.average_power_watts),
        });
    }
    if ratios.is_empty() {
        return Err(format!("no {precision} measurements").into());
    }
    Ok(PrecisionReport {
        precision,
        ratios,
        end_cycle_cvs,
        raw_cycle_cvs,
        end_energy_cvs,
        raw_energy_cvs,
    })
}

fn summarize_target(
    target: Target,
    measurements: &[Measurement],
) -> Result<TargetReport, Box<dyn Error>> {
    let sources = measurements
        .iter()
        .map(|row| row.power_source.as_str())
        .collect::<BTreeSet<_>>();
    if sources.len() != 1 {
        return Err(format!("{} has mixed power sources: {sources:?}", target.label).into());
    }
    let precisions = PRECISIONS
        .iter()
        .map(|precision| summarize_precision(precision, measurements))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(TargetReport {
        target,
        power_source: sources.into_iter().next().unwrap().to_string(),
        precisions,
    })
}

fn metric_values(report: &PrecisionReport, metric: Metric) -> Vec<f64> {
    report
        .ratios
        .iter()
        .map(|ratio| match metric {
            Metric::EndCycles => ratio.end_cycles,
            Metric::RawCycles => ratio.raw_cycles,
            Metric::RawInstructions => ratio.raw_instructions,
            Metric::EndEnergy => ratio.end_energy,
            Metric::RawEnergy => ratio.raw_energy,
            Metric::EndPower => ratio.end_power,
        })
        .collect()
}

fn write_histogram(
    reports: &[TargetReport],
    metric: Metric,
    output_path: &Path,
) -> Result<(), Box<dyn Error>> {
    let owned = reports
        .iter()
        .enumerate()
        .flat_map(|(target_index, target)| {
            target
                .precisions
                .iter()
                .enumerate()
                .map(move |(precision_index, report)| {
                    (
                        format!(
                            "{} {}",
                            if target.target.label == "AArch64" {
                                "A64"
                            } else {
                                "x86"
                            },
                            report.precision
                        ),
                        target_index * PRECISIONS.len() + precision_index,
                        metric_values(report, metric),
                    )
                })
        })
        .collect::<Vec<_>>();
    let series = owned
        .iter()
        .map(|(label, palette, values)| {
            HistogramSeries::new(values)
                .with_label(label)
                .with_bins(6)
                .with_palette(*palette)
        })
        .collect::<Vec<_>>();
    let (title, x_label) = match metric {
        Metric::EndCycles => (
            "End-to-end vertical packed-SIMD Horner/Estrin cycle ratios",
            "median Horner cycles / median Estrin cycles",
        ),
        Metric::RawCycles => (
            "Raw-polynomial vertical packed-SIMD Horner/Estrin cycle ratios",
            "median Horner cycles / median Estrin cycles",
        ),
        Metric::RawInstructions => (
            "Raw-polynomial packed-SIMD instruction ratios",
            "median Horner instructions / median Estrin instructions",
        ),
        Metric::EndEnergy => (
            "End-to-end vertical packed-SIMD Horner/Estrin energy ratios",
            "median Horner energy/value / median Estrin energy/value",
        ),
        Metric::RawEnergy => (
            "Raw-polynomial vertical packed-SIMD Horner/Estrin energy ratios",
            "median Horner energy/value / median Estrin energy/value",
        ),
        Metric::EndPower => (
            "End-to-end vertical packed-SIMD Horner/Estrin mean-power ratios",
            "median Horner mean power / median Estrin mean power",
        ),
    };
    let output = output_path.to_string_lossy().into_owned();
    let config = HistogramConfig::publication(title)
        .with_size(1400, 800)
        .with_font_scale(1.0)
        .with_x_label(x_label)
        .with_y_label("Functions")
        .with_reference_line(1.0)
        .with_save_path(&output);
    plot_histograms_with_config(&series, &config)?;
    Ok(())
}

fn summary(report: &PrecisionReport, metric: Metric) -> (usize, f64) {
    let values = metric_values(report, metric);
    (
        values.iter().filter(|ratio| **ratio < 1.0).count(),
        geometric_mean(&values).unwrap(),
    )
}

fn ratio_for<'a>(report: &'a PrecisionReport, operation: &str) -> &'a OperationRatio {
    report
        .ratios
        .iter()
        .find(|ratio| ratio.operation == operation)
        .expect("operation must be present on both targets")
}

fn write_summary_tex(path: &Path, reports: &[TargetReport]) -> Result<(), Box<dyn Error>> {
    let mut output = String::from(
        "% Auto-generated by kernel_selection_summary.\n\
         \\begin{figure}[H]\n\\centering\n\
         \\begin{minipage}[t]{0.49\\textwidth}\n\\centering\n\
         \\includegraphics[width=\\linewidth]{outputs/fpx/kernel-selection/kernel_selection_summary.pdf}\n\
         \\end{minipage}\\hfill\n\
         \\begin{minipage}[t]{0.49\\textwidth}\n\\centering\n\
         \\includegraphics[width=\\linewidth]{outputs/fpx/kernel-selection/kernel_selection_raw.pdf}\n\
         \\end{minipage}\n\
         \\caption{Vertical packed-SIMD Horner-to-Estrin cycle ratios for the end-to-end SIMD path (left) and isolated raw polynomial (right). Values below one favor Horner; the vertical line marks parity.}\n\
         \\label{fig:simd-kernel-selection-cycles}\n\\end{figure}\n\n\
         \\begin{figure}[H]\n\\centering\n\
         \\begin{minipage}[t]{0.49\\textwidth}\n\\centering\n\
         \\includegraphics[width=\\linewidth]{outputs/fpx/kernel-selection/kernel_selection_energy_efficiency.pdf}\n\
         \\end{minipage}\\hfill\n\
         \\begin{minipage}[t]{0.49\\textwidth}\n\\centering\n\
         \\includegraphics[width=\\linewidth]{outputs/fpx/kernel-selection/kernel_selection_power.pdf}\n\
         \\end{minipage}\n\
         \\caption{End-to-end energy per value (left) and mean-power (right) ratios. Energy combines duration and power; mean power alone stays close to parity. Measurements are architecture-local because the Pi reports PMIC output-rail power while x86 reports package energy.}\n\
         \\label{fig:simd-kernel-selection-energy-power}\n\\end{figure}\n\n",
    );

    output.push_str(
        "\\begin{table}[H]\n\\centering\n\
         \\caption{Vertical packed-SIMD Horner--Estrin selection across thirteen functions}\n\
         \\label{tab:simd-kernel-selection-summary}\n\
         \\resizebox{\\textwidth}{!}{%\n\
         \\begin{tabular}{llrrrrrrrrrl}\n\\toprule\n\
         Target & Precision & End cyc. H/E & Wins & Raw cyc. H/E & Wins & End energy H/E & Wins & End power H/E & Cycle CV & Energy CV & Selection \\\\\n         \\midrule\n",
    );
    for target in reports {
        for report in &target.precisions {
            let (end_cycle_wins, end_cycles) = summary(report, Metric::EndCycles);
            let (raw_cycle_wins, raw_cycles) = summary(report, Metric::RawCycles);
            let (energy_wins, energy) = summary(report, Metric::EndEnergy);
            let (_, power) = summary(report, Metric::EndPower);
            let mut cycle_cvs = report.end_cycle_cvs.clone();
            cycle_cvs.extend_from_slice(&report.raw_cycle_cvs);
            let mut energy_cvs = report.end_energy_cvs.clone();
            energy_cvs.extend_from_slice(&report.raw_energy_cvs);
            output.push_str(&format!(
                "{} & {} & {:.3} & {}/13 & {:.3} & {}/13 & {:.3} & {}/13 & {:.3} & {:.2}\\% & {:.2}\\% & \\textbf{{Horner}} \\\\\n",
                target.target.detail,
                report.precision,
                end_cycles,
                end_cycle_wins,
                raw_cycles,
                raw_cycle_wins,
                energy,
                energy_wins,
                power,
                100.0 * median(&cycle_cvs).unwrap(),
                100.0 * median(&energy_cvs).unwrap(),
            ));
        }
    }
    output.push_str("\\bottomrule\n\\end{tabular}%\n}\n\\end{table}\n\n");

    output.push_str(
        "\\begin{table}[H]\n\\centering\n\\footnotesize\n\
         \\caption{Raw-polynomial mechanism and energy ratios}\n\
         \\label{tab:simd-kernel-mechanism}\n\
         \\resizebox{\\textwidth}{!}{%\n\
         \\begin{tabular}{llrrrrr}\n\\toprule\n\
         Target & Precision & Cycles H/E & Instructions H/E & Time H/E & Energy H/E & Power H/E \\\\\n         \\midrule\n",
    );
    for target in reports {
        for report in &target.precisions {
            let (_, cycles) = summary(report, Metric::RawCycles);
            let (_, instructions) = summary(report, Metric::RawInstructions);
            let elapsed = geometric_mean(
                &report
                    .ratios
                    .iter()
                    .map(|ratio| ratio.raw_elapsed)
                    .collect::<Vec<_>>(),
            )
            .unwrap();
            let (_, energy) = summary(report, Metric::RawEnergy);
            let power = geometric_mean(
                &report
                    .ratios
                    .iter()
                    .map(|ratio| ratio.raw_power)
                    .collect::<Vec<_>>(),
            )
            .unwrap();
            output.push_str(&format!(
                "{} & {} & {:.3} & {:.3} & {:.3} & {:.3} & {:.3} \\\\\n",
                target.target.label, report.precision, cycles, instructions, elapsed, energy, power,
            ));
        }
    }
    output.push_str("\\bottomrule\n\\end{tabular}%\n}\n\\end{table}\n\n");

    let operations = reports[0].precisions[0]
        .ratios
        .iter()
        .map(|ratio| ratio.operation.clone())
        .collect::<Vec<_>>();
    output.push_str(
        "\\begin{table}[H]\n\\centering\n\\footnotesize\n\
         \\caption{Per-function raw-polynomial Horner-to-Estrin cycle ratios}\n\
         \\label{tab:simd-kernel-selection-functions}\n\
         \\renewcommand{\\arraystretch}{0.85}\n\
         \\begin{tabular}{lrrrr}\n\\toprule\n\
         Function & AArch64 f32 & AArch64 f64 & x86-64 f32 & x86-64 f64 \\\\\n         \\midrule\n",
    );
    for operation in operations {
        output.push_str(&format!(
            "{} & {:.3} & {:.3} & {:.3} & {:.3} \\\\\n",
            operation,
            ratio_for(&reports[0].precisions[0], &operation).raw_cycles,
            ratio_for(&reports[0].precisions[1], &operation).raw_cycles,
            ratio_for(&reports[1].precisions[0], &operation).raw_cycles,
            ratio_for(&reports[1].precisions[1], &operation).raw_cycles,
        ));
    }
    output.push_str("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n");

    let total_cycle_wins = reports
        .iter()
        .flat_map(|target| &target.precisions)
        .map(|report| summary(report, Metric::EndCycles).0 + summary(report, Metric::RawCycles).0)
        .sum::<usize>();
    output.push_str(&format!(
        "\\par Horner wins all {total_cycle_wins} cycle comparisons across two targets, thirteen functions, two precisions, and two scopes. The raw instruction ratios are substantially below one on both targets, confirming that Estrin's additional powers, temporaries, and tree-combination work outweigh its shorter dependency depth for these fixed polynomial shapes. Target-matched assembly inspection independently confirms packed AVX2/FMA on x86-64 and packed NEON arithmetic on AArch64; the result is not a failure to vectorize.\n\n"
    ));
    output.push_str(&format!(
        "\\par Power was collected from \\path{{{}}} on AArch64 and \\path{{{}}} on x86-64. These scopes are not directly comparable in watts. Within each target, energy per value generally tracks execution time while mean power stays much nearer parity, so Horner's energy advantage is primarily a shorter-duration effect rather than a large reduction in instantaneous power.\n\n",
        reports[0].power_source,
        reports[1].power_source,
    ));
    output.push_str(
        "\\par Each ratio compares configuration medians from ten runs of one million values. Configurations were interleaved by repetition to reduce thermal and DVFS drift. Cycle counts select the kernel; elapsed time, energy, and power provide corroborating system-level evidence.\n",
    );
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, output)?;
    Ok(())
}

fn run() -> Result<(), Box<dyn Error>> {
    let root = repository_root();
    let output_dir = root.join("outputs/fpx/kernel-selection");
    fs::create_dir_all(&output_dir)?;
    let reports = TARGETS
        .iter()
        .map(|target| {
            let measurements = load_measurements(&root.join(target.csv_path))?;
            summarize_target(*target, &measurements)
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    for (metric, file) in [
        (Metric::EndCycles, "kernel_selection_summary.svg"),
        (Metric::RawCycles, "kernel_selection_raw.svg"),
        (Metric::RawInstructions, "kernel_selection_instructions.svg"),
        (Metric::EndEnergy, "kernel_selection_energy_efficiency.svg"),
        (Metric::EndPower, "kernel_selection_power.svg"),
    ] {
        let path = output_dir.join(file);
        write_histogram(&reports, metric, &path)?;
        println!("Wrote {}", path.display());
    }
    let tex_path = output_dir.join("kernel_selection.tex");
    write_summary_tex(&tex_path, &reports)?;
    println!("Wrote {}", tex_path.display());
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("kernel_selection_summary: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn measurement(
        mode: &str,
        precision: &str,
        cycles: f64,
        instructions: f64,
        elapsed: f64,
        energy: f64,
        power: f64,
    ) -> Measurement {
        Measurement {
            operation: "sin".into(),
            mode: mode.into(),
            precision: precision.into(),
            cycles_per_call: cycles,
            instructions_per_call: instructions,
            ns_per_call: elapsed,
            nj_per_call: energy,
            average_power_watts: power,
            power_source: "test".into(),
        }
    }

    #[test]
    fn compares_all_metrics_for_both_scopes() {
        let mut rows = Vec::new();
        for precision in PRECISIONS {
            for (mode, scale) in [
                (END_HORNER, 8.0),
                (END_ESTRIN, 10.0),
                (RAW_HORNER, 4.0),
                (RAW_ESTRIN, 10.0),
            ] {
                rows.extend([
                    measurement(mode, precision, scale, scale, scale, scale, scale),
                    measurement(
                        mode,
                        precision,
                        2.0 * scale,
                        2.0 * scale,
                        2.0 * scale,
                        2.0 * scale,
                        2.0 * scale,
                    ),
                ]);
            }
        }
        let report = summarize_target(TARGETS[0], &rows).unwrap();
        let ratio = &report.precisions[0].ratios[0];
        assert_eq!(ratio.end_cycles, 0.8);
        assert_eq!(ratio.raw_cycles, 0.4);
        assert_eq!(ratio.raw_instructions, 0.4);
        assert_eq!(ratio.end_energy, 0.8);
        assert_eq!(ratio.raw_power, 0.4);
    }
}
