use std::env;
use std::fs::{self, OpenOptions};
#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

use aether_fpx_bench::{
    export::write_csv,
    format::{fmt_f64, fmt_i64},
    perf::{parse_perf_stat_csv, run_perf_stat},
    pi_pmic::{read_pmic_rail_power_watts, run_pi_pmic_power},
    report::{print_comparison, print_summary},
    sysfs_power::SysfsPower,
    types::Row,
};

// trig_probe argv: <op> <impl> <precision> [iters]
const OPS: &[&str] = &[
    "tan", "sin", "cos", "asin", "acos", "atan", "sinh", "cosh", "tanh", "exp", "log", "log1p",
    "expm1",
];
const MODES: &[&str] = &[
    "fpx",
    "fpx_blocked_horner",
    "fpx_horizontal_estrin",
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
    "nat",
    "x87",
];
const PRECISIONS: &[&str] = &["f16", "f32", "f64", "f128"];
const X87_AVAILABLE: bool = cfg!(fpx_x87);

// perf events to collect
const EVENTS: &[&str] = &["instructions", "cycles", "branches", "branch-misses"];

const RUNS: usize = 10;
const ITERS: usize = 100_000;
const PMIC_RUNS: usize = 10;
const PMIC_ITERS: usize = 100_000;
const PMIC_SAMPLE_MS: u64 = 20;
const PERF_CSV: &str = "perf_trig_probe.csv";
const PMIC_CSV: &str = "target/perf_trig_probe_pi_pmic.csv";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum PowerSource {
    Perf,
    PiPmic,
}

#[derive(Debug)]
struct Args {
    trig_probe: PathBuf,
    power_source: PowerSource,
    runs: usize,
    iters: usize,
    power_runs: usize,
    power_iters: usize,
    sample_ms: u64,
    ops: Vec<String>,
    modes: Vec<String>,
    precisions: Vec<String>,
    require_power: bool,
    out: PathBuf,
}

struct StagedElf {
    path: PathBuf,
}

impl StagedElf {
    fn as_path(&self) -> &Path {
        &self.path
    }
}

impl Drop for StagedElf {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

fn parse_args() -> Result<Args, String> {
    let mut power_source = PowerSource::Perf;
    let mut runs = RUNS;
    let mut iters = ITERS;
    let mut power_runs = PMIC_RUNS;
    let mut power_iters = PMIC_ITERS;
    let mut sample_ms = PMIC_SAMPLE_MS;
    let mut ops = OPS.iter().map(|value| (*value).to_string()).collect();
    // Horizontal modes remain available through --modes, but the thesis-facing
    // default capture compares the vertical packed-input paths.
    let mut modes = MODES
        .iter()
        .filter(|&&mode| !matches!(mode, "fpx_blocked_horner" | "fpx_horizontal_estrin"))
        .map(|value| (*value).to_string())
        .collect();
    let mut precisions = PRECISIONS
        .iter()
        .map(|value| (*value).to_string())
        .collect();
    let mut require_power = false;
    let mut out = None;
    let mut trig_probe = None;
    let args = env::args().skip(1).collect::<Vec<_>>();
    let mut index = 0usize;
    while index < args.len() {
        let flag = args[index].as_str();
        let value = |index: usize| {
            args.get(index)
                .map(String::as_str)
                .ok_or_else(|| format!("missing value for {flag}"))
        };
        match flag {
            "--power-source" => {
                index += 1;
                power_source = match value(index)? {
                    "perf" => PowerSource::Perf,
                    "pi-pmic" => PowerSource::PiPmic,
                    other => return Err(format!("unknown power source: {other}")),
                };
            }
            "--power-runs" => {
                index += 1;
                power_runs = value(index)?.parse().map_err(|_| "invalid --power-runs")?;
            }
            "--runs" => {
                index += 1;
                runs = value(index)?.parse().map_err(|_| "invalid --runs")?;
            }
            "--iters" => {
                index += 1;
                iters = value(index)?.parse().map_err(|_| "invalid --iters")?;
            }
            "--power-iters" => {
                index += 1;
                power_iters = value(index)?.parse().map_err(|_| "invalid --power-iters")?;
            }
            "--sample-ms" => {
                index += 1;
                sample_ms = value(index)?.parse().map_err(|_| "invalid --sample-ms")?;
            }
            "--ops" => {
                index += 1;
                ops = parse_list(value(index)?);
            }
            "--modes" => {
                index += 1;
                modes = parse_list(value(index)?);
            }
            "--precisions" => {
                index += 1;
                precisions = parse_list(value(index)?);
            }
            "--require-power" => require_power = true,
            "--out" => {
                index += 1;
                out = Some(PathBuf::from(value(index)?));
            }
            "-h" | "--help" => {
                println!(
                    r#"usage: perf [--power-source perf|pi-pmic] [--runs N] [--iters N] \
  [--power-runs N] [--power-iters N] [--sample-ms N] \
  [--ops comma,list] [--modes comma,list] [--precisions comma,list] [--require-power] \
  [--out CSV] <absolute trig_probe>

perf default: target/{}
pi-pmic default: {}
PBS defaults include $PBS_JOBID to avoid job collisions. --out takes precedence.
pi-pmic settings: runs={}, iters={}, sample-ms={}"#,
                    PERF_CSV, PMIC_CSV, PMIC_RUNS, PMIC_ITERS, PMIC_SAMPLE_MS,
                );
                std::process::exit(0);
            }
            path if path.starts_with('-') => return Err(format!("unknown argument: {path}")),
            path if trig_probe.is_some() => {
                return Err(format!("unexpected trig_probe path: {path}"));
            }
            path => trig_probe = Some(PathBuf::from(path)),
        }
        index += 1;
    }
    if runs == 0 || iters == 0 || power_runs == 0 || power_iters == 0 || sample_ms == 0 {
        return Err("runs, iterations, and sample interval must be non-zero".into());
    }
    validate_selection("--ops", &ops, OPS)?;
    validate_selection("--modes", &modes, MODES)?;
    validate_selection("--precisions", &precisions, PRECISIONS)?;
    let pbs_job_id = env::var("PBS_JOBID").ok();
    let default_out = default_output_path(power_source, pbs_job_id.as_deref());
    let trig_probe = trig_probe.ok_or("missing absolute path to trig_probe")?;
    if !trig_probe.is_absolute() {
        return Err(format!(
            "trig_probe path must be absolute: {}",
            trig_probe.display()
        ));
    }
    if !trig_probe.is_file() {
        return Err(format!(
            "trig_probe is not a file: {}",
            trig_probe.display()
        ));
    }
    Ok(Args {
        trig_probe,
        power_source,
        runs,
        iters,
        power_runs,
        power_iters,
        sample_ms,
        ops,
        modes,
        precisions,
        require_power,
        out: out.unwrap_or(default_out),
    })
}

fn parse_list(value: &str) -> Vec<String> {
    value
        .split(',')
        .map(str::trim)
        .filter(|item| !item.is_empty())
        .map(str::to_owned)
        .collect()
}

fn validate_selection(flag: &str, values: &[String], allowed: &[&str]) -> Result<(), String> {
    if values.is_empty()
        || values
            .iter()
            .any(|value| !allowed.contains(&value.as_str()))
    {
        return Err(format!("{flag} must contain only: {}", allowed.join(",")));
    }
    Ok(())
}

fn default_output_path(power_source: PowerSource, pbs_job_id: Option<&str>) -> PathBuf {
    let default = match power_source {
        PowerSource::Perf => PERF_CSV,
        PowerSource::PiPmic => "perf_trig_probe_pi_pmic.csv",
    };
    let output = PathBuf::from("target").join(default);
    let Some(job_id) = pbs_job_id
        .map(sanitize_output_component)
        .filter(|job_id| !job_id.is_empty())
    else {
        return output;
    };

    let stem = output
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("perf_trig_probe");
    output.with_file_name(format!("{stem}_{job_id}.csv"))
}

fn sanitize_output_component(value: &str) -> String {
    value
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || matches!(character, '.' | '_' | '-') {
                character
            } else {
                '_'
            }
        })
        .collect()
}

fn prepare_output_path(out: &Path) -> Result<(), String> {
    if let Some(parent) = out.parent().filter(|parent| !parent.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|error| {
            format!(
                "failed to create CSV output directory {}: {error}",
                parent.display()
            )
        })?;
    }

    OpenOptions::new()
        .create(true)
        .append(true)
        .open(out)
        .map(|_| ())
        .map_err(|error| format!("cannot write CSV output {}: {error}", out.display()))
}

fn write_results_csv(out: &Path, rows: &[Row]) -> Result<(), String> {
    prepare_output_path(out)?;
    write_csv(out, rows)
        .map_err(|error| format!("failed to write CSV output {}: {error}", out.display()))
}

fn write_results_csv_or_exit(out: &Path, rows: &[Row]) {
    write_results_csv(out, rows).unwrap_or_else(|error| {
        eprintln!("{error}");
        std::process::exit(1);
    });
}

fn stage_elf_for_execution(source: &Path) -> Result<StagedElf, String> {
    let temp_root = env::var_os("PBS_TMPDIR")
        .or_else(|| env::var_os("TMPDIR"))
        .map(PathBuf::from)
        .unwrap_or_else(env::temp_dir);
    fs::create_dir_all(&temp_root).map_err(|error| {
        format!(
            "failed to create staging directory {}: {error}",
            temp_root.display()
        )
    })?;

    let file_name = source
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("trig_probe");
    let job = env::var("PBS_JOBID")
        .map(|job| sanitize_output_component(&job))
        .unwrap_or_else(|_| "local".to_string());
    let staged = temp_root.join(format!("fpx-{job}-{}-{file_name}", std::process::id()));
    fs::copy(source, &staged).map_err(|error| {
        format!(
            "failed to stage benchmark ELF {} to {}: {error}",
            source.display(),
            staged.display()
        )
    })?;
    make_executable(source, &staged)?;
    Ok(StagedElf { path: staged })
}

#[cfg(unix)]
fn make_executable(source: &Path, staged: &Path) -> Result<(), String> {
    let mut permissions = fs::metadata(source)
        .map_err(|error| {
            format!(
                "failed to read ELF permissions {}: {error}",
                source.display()
            )
        })?
        .permissions();
    permissions.set_mode(permissions.mode() | 0o500);
    fs::set_permissions(staged, permissions).map_err(|error| {
        format!(
            "failed to set executable permissions on {}: {error}",
            staged.display()
        )
    })
}

#[cfg(not(unix))]
fn make_executable(_source: &Path, _staged: &Path) -> Result<(), String> {
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        default_output_path, parse_energy_preflight, parse_perf_list_energy_event,
        sanitize_output_component, supports_mode, PowerSource, OPS,
    };
    use std::path::PathBuf;

    #[test]
    fn keeps_the_local_default_output_name() {
        assert_eq!(
            default_output_path(PowerSource::Perf, None),
            PathBuf::from("target/perf_trig_probe.csv")
        );
    }

    #[test]
    fn scopes_pbs_output_by_job_id() {
        assert_eq!(
            default_output_path(PowerSource::PiPmic, Some("52.mercury")),
            PathBuf::from("target/perf_trig_probe_pi_pmic_52.mercury.csv")
        );
    }

    #[test]
    fn sanitizes_job_ids_for_staged_paths() {
        assert_eq!(sanitize_output_component("52.mercury/0"), "52.mercury_0");
    }

    #[test]
    fn detects_package_energy_from_perf_list() {
        let event = parse_perf_list_energy_event(
            r#"
  power/energy-pkg/                 [Kernel PMU event]
  power/event=0..255/modifier       [Raw event descriptor]
"#,
        );
        assert_eq!(event.as_deref(), Some("power/energy-pkg/"));
    }

    #[test]
    fn prefers_package_energy_from_perf_list() {
        let event = parse_perf_list_energy_event(
            r#"
  power/energy-cores/               [Kernel PMU event]
  power/energy-pkg/                 [Kernel PMU event]
"#,
        );
        assert_eq!(event.as_deref(), Some("power/energy-pkg/"));
    }

    #[test]
    fn diagnoses_non_numeric_perf_energy_preflight() {
        let error = parse_energy_preflight(
            "power/energy-pkg/",
            "<not counted>,Joules,power/energy-pkg/,0,0.00,,",
        )
        .unwrap_err();

        assert!(error.contains("<not counted>"));
    }

    #[test]
    fn all_ops_expose_kernel_paths_for_f32_and_f64() {
        for op in OPS {
            for precision in ["f16", "f32", "f64", "f128"] {
                assert!(supports_mode(op, "fpx", precision));
                assert_eq!(
                    supports_mode(op, "fpx_blocked_horner", precision),
                    matches!(precision, "f32" | "f64")
                );
                assert_eq!(
                    supports_mode(op, "fpx_horizontal_estrin", precision),
                    matches!(precision, "f32" | "f64")
                );
                assert!(supports_mode(op, "nat", precision));
                assert_eq!(supports_mode(op, "libm", precision), precision != "f128");
                assert_eq!(
                    supports_mode(op, "fpx_kernel_horner", precision),
                    matches!(precision, "f32" | "f64")
                );
                assert_eq!(
                    supports_mode(op, "fpx_kernel_estrin", precision),
                    matches!(precision, "f32" | "f64")
                );
                for mode in [
                    "fpx_kernel_scalar",
                    "fpx_kernel_fma",
                    "fpx_kernel_simd",
                    "fpx_kernel_simd_estrin",
                    "fpx_kernel_raw_scalar",
                    "fpx_kernel_raw_fma",
                    "fpx_kernel_raw_simd",
                    "fpx_kernel_raw_simd_estrin",
                ] {
                    assert_eq!(
                        supports_mode(op, mode, precision),
                        matches!(precision, "f32" | "f64")
                    );
                }
            }
        }
    }
}

fn discover_energy_event() -> Option<String> {
    if let Ok(configured) = env::var("FPX_PERF_ENERGY_EVENT") {
        let configured = configured.trim();
        if configured.is_empty()
            || configured.eq_ignore_ascii_case("none")
            || configured.eq_ignore_ascii_case("off")
        {
            return None;
        }
        return Some(configured.to_string());
    }

    discover_energy_event_from_perf_list().or_else(discover_energy_event_from_sysfs)
}

fn discover_energy_event_from_perf_list() -> Option<String> {
    let output = Command::new("perf")
        .env("LC_ALL", "C")
        .arg("list")
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }

    parse_perf_list_energy_event(&String::from_utf8_lossy(&output.stdout))
}

fn parse_perf_list_energy_event(text: &str) -> Option<String> {
    let mut candidates = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let Some(event) = line.split_whitespace().find_map(clean_perf_list_event) else {
            continue;
        };
        candidates.push((energy_event_priority(event), event.to_string()));
    }
    candidates.sort();
    candidates.into_iter().next().map(|(_, event)| event)
}

fn clean_perf_list_event(token: &str) -> Option<&str> {
    let token = token.trim_end_matches(',');
    let lower = token.to_ascii_lowercase();
    (lower.contains("/energy") && lower.ends_with('/') && !lower.contains("event="))
        .then_some(token)
}

fn discover_energy_event_from_sysfs() -> Option<String> {
    let devices = fs::read_dir("/sys/bus/event_source/devices").ok()?;
    let mut candidates = Vec::new();
    for device in devices.flatten() {
        let device_name = device.file_name().to_string_lossy().into_owned();
        let events_path = device.path().join("events");
        let Ok(events) = fs::read_dir(&events_path) else {
            continue;
        };
        for event in events.flatten() {
            let event_name = event.file_name().to_string_lossy().into_owned();
            if !event_name.to_ascii_lowercase().contains("energy") || event_name.contains('.') {
                continue;
            }

            let unit_path = events_path.join(format!("{event_name}.unit"));
            if let Ok(unit) = fs::read_to_string(unit_path) {
                if !unit.trim().eq_ignore_ascii_case("joules") {
                    continue;
                }
            }

            let event = format!("{device_name}/{event_name}/");
            candidates.push((energy_event_priority(&event), event));
        }
    }
    candidates.sort();
    candidates.into_iter().next().map(|(_, event)| event)
}

fn energy_event_priority(event: &str) -> u8 {
    let event = event.to_ascii_lowercase();
    if event.contains("energy-pkg") {
        0
    } else if event.contains("energy-psys") {
        1
    } else if event.contains("energy-cores") {
        2
    } else if event.contains("energy-ram") || event.contains("energy-dram") {
        3
    } else {
        4
    }
}

fn validate_perf_energy_event(event: &str) -> Result<f64, String> {
    let output = Command::new("perf")
        .env("LC_ALL", "C")
        .arg("stat")
        .arg("-x")
        .arg(",")
        .arg("--no-big-num")
        .arg("-e")
        .arg(event)
        .arg("--")
        .arg("sleep")
        .arg("0.05")
        .output()
        .map_err(|error| format!("failed to run perf energy preflight: {error}"))?;
    let stderr = String::from_utf8_lossy(&output.stderr);

    if !output.status.success() {
        return Err(first_relevant_perf_line(&stderr, event)
            .unwrap_or_else(|| stderr.trim().to_string())
            .trim()
            .to_string());
    }

    parse_energy_preflight(event, &stderr)
}

fn parse_energy_preflight(event: &str, stderr: &str) -> Result<f64, String> {
    parse_perf_stat_csv(stderr)
        .iter()
        .find_map(|(name, value)| {
            (name == event || name.to_ascii_lowercase().contains("energy")).then_some(*value)
        })
        .ok_or_else(|| {
            first_relevant_perf_line(stderr, event)
                .unwrap_or_else(|| "event did not appear in perf stat output".to_string())
        })
}

fn first_relevant_perf_line(stderr: &str, event: &str) -> Option<String> {
    stderr
        .lines()
        .map(str::trim)
        .find(|line| line.contains(event) || line.to_ascii_lowercase().contains("energy"))
        .map(ToOwned::to_owned)
}

fn supports_mode(op: &str, mode: &str, prec: &str) -> bool {
    match mode {
        "x87" => X87_AVAILABLE && prec == "f64" && matches!(op, "sin" | "cos" | "tan"),
        "libm" => prec != "f128",
        "fpx_blocked_horner"
        | "fpx_horizontal_estrin"
        | "fpx_kernel_horner"
        | "fpx_kernel_estrin"
        | "fpx_kernel_scalar"
        | "fpx_kernel_fma"
        | "fpx_kernel_simd"
        | "fpx_kernel_simd_estrin"
        | "fpx_kernel_raw_scalar"
        | "fpx_kernel_raw_fma"
        | "fpx_kernel_raw_simd"
        | "fpx_kernel_raw_simd_estrin" => matches!(prec, "f32" | "f64"),
        _ => true,
    }
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(4)
        .expect("aether-fpx-research should live under the workspace root")
        .to_path_buf()
}

fn find_release_elf_in(workspace_root: &Path, bin_name: &str) -> Option<PathBuf> {
    let plain = workspace_root.join("target").join("release").join(bin_name);
    if plain.exists() {
        return Some(plain);
    }

    let deps = workspace_root.join("target").join("release").join("deps");
    let rd = fs::read_dir(&deps).ok()?;

    for ent in rd.flatten() {
        let p = ent.path();
        if !p.is_file() {
            continue;
        }
        if let Some(fname) = p.file_name().and_then(|s| s.to_str()) {
            if fname.starts_with(&format!("{bin_name}-")) && p.extension().is_none() {
                return Some(p);
            }
        }
    }

    None
}

fn build_release_elf(bin_name: &str) -> Result<PathBuf, String> {
    let root = workspace_root();
    if let Some(elf) = find_release_elf_in(&root, bin_name) {
        return Ok(elf);
    }

    let mut command = Command::new("cargo");
    command.current_dir(&root).args([
        "build",
        "--release",
        "-p",
        "aether-fpx-research",
        "--bin",
        bin_name,
    ]);
    if !X87_AVAILABLE {
        command.arg("--no-default-features");
    }
    let status = command.status().map_err(|error| {
        format!(
            "{bin_name} is not built and Cargo could not be invoked: {error}. \
                 Build it before using sudo with: cargo build --release -p aether-fpx-research \
                 --bin perf --bin {bin_name}"
        )
    })?;

    if !status.success() {
        return Err(format!(
            "cargo build --release for {bin_name} failed with status {status}"
        ));
    }

    find_release_elf_in(&root, bin_name)
        .ok_or_else(|| format!("built {bin_name}, but could not locate the resulting ELF"))
}

fn run_pi_pmic_collection(elf: &Path, args: &Args) -> Result<(), String> {
    println!(
        "Pi PMIC rail collection: {} runs, {} calls/run, {} ms sampling",
        args.power_runs, args.power_iters, args.sample_ms
    );
    println!("Libraries are interleaved by repetition to reduce thermal and DVFS drift.");
    let interval = Duration::from_millis(args.sample_ms);
    let mut rows = Vec::new();

    for run in 0..args.power_runs {
        for op in args.ops.iter().map(String::as_str) {
            for prec in args.precisions.iter().map(String::as_str) {
                for mode in args.modes.iter().map(String::as_str) {
                    if !supports_mode(op, mode, prec) {
                        continue;
                    }
                    let row =
                        run_pi_pmic_power(elf, op, mode, prec, args.power_iters, run, interval)?;
                    let elapsed = row.time_elapsed_s.unwrap_or_default();
                    if elapsed < 3.0 * interval.as_secs_f64() {
                        eprintln!(
                            "warning: {op}/{mode}/{prec} ran for {elapsed:.3}s; increase --power-iters for reliable PMIC sampling"
                        );
                    }
                    println!(
                        "run {run:2}: op={op} mode={mode} prec={prec} elapsed={elapsed:.3}s rail-energy(J)={} rail-power(W)={}",
                        fmt_f64(row.energy_joules),
                        fmt_f64(row.average_power_watts),
                    );
                    rows.push(row);
                }
            }
        }
    }

    write_results_csv(&args.out, &rows)?;
    println!("Wrote {}", args.out.display());
    println!(
        "Generate tables with --iters {} so energy is normalized by the PMIC workload size.",
        args.power_iters
    );
    Ok(())
}

fn main() {
    let args = parse_args().unwrap_or_else(|error| {
        eprintln!("{error}");
        std::process::exit(2);
    });
    let staged_elf = stage_elf_for_execution(&args.trig_probe).unwrap_or_else(|error| {
        eprintln!("{error}");
        std::process::exit(2);
    });
    let elf = staged_elf.as_path();

    println!("ELF: {}", args.trig_probe.display());
    println!("Staged ELF: {}", elf.display());
    prepare_output_path(&args.out).unwrap_or_else(|error| {
        eprintln!("{error}");
        eprintln!("hint: select a writable location with --out <CSV>");
        std::process::exit(2);
    });
    println!("CSV output: {}", args.out.display());

    if args.power_source == PowerSource::PiPmic {
        run_pi_pmic_collection(&elf, &args).unwrap_or_else(|error| {
            eprintln!("Pi PMIC collection failed: {error}");
            std::process::exit(1);
        });
        return;
    }

    let mut rows: Vec<Row> = Vec::new();
    let discovered_energy_event = discover_energy_event();
    let energy_event =
        discovered_energy_event.and_then(|event| match validate_perf_energy_event(&event) {
            Ok(joules) => {
                println!("Energy counter: {event} (preflight {joules:.6} J)");
                Some(event)
            }
            Err(error) => {
                println!("Energy counter: {event} (unusable: {error})");
                None
            }
        });
    let mut events = EVENTS
        .iter()
        .map(|event| (*event).to_string())
        .collect::<Vec<_>>();
    if let Some(event) = &energy_event {
        events.push(event.clone());
    } else {
        println!(
            "Energy counter: unavailable (set FPX_PERF_ENERGY_EVENT to override or 'off' to disable)"
        );
        match SysfsPower::autodetect() {
            Some(source) => println!("Sysfs power fallback: {}", source.source()),
            None => {
                println!(
                    "Sysfs power fallback: unavailable (no supported /sys/class/powercap or /sys/class/hwmon CPU power source found)"
                );
                match read_pmic_rail_power_watts() {
                    Ok(watts) => println!("Pi PMIC fallback: available ({watts:.3} W)"),
                    Err(error) => println!("Pi PMIC fallback: unavailable ({error})"),
                }
            }
        }
    }
    let event_refs = events.iter().map(String::as_str).collect::<Vec<_>>();

    println!("Configurations are interleaved by repetition to reduce thermal and DVFS drift.");
    for run in 0..args.runs {
        for op in args.ops.iter().map(String::as_str) {
            for prec in args.precisions.iter().map(String::as_str) {
                for mode in args.modes.iter().map(String::as_str) {
                    if !supports_mode(op, mode, prec) {
                        continue;
                    }

                    println!(
                        "\n== perf stat run={run} op={op} mode={mode} prec={prec} (iters={}) ==",
                        args.iters
                    );

                    let row = run_perf_stat(
                        elf, op, mode, prec, args.iters, run, &event_refs,
                    ).unwrap_or_else(|e| {
                        eprintln!("perf failed: {e}");
                        if e.contains("No permission to enable")
                            || e.contains("Permission denied")
                            || e.contains("perf_event_paranoid")
                        {
                            eprintln!("hint: run crates/aether-fpx/support/research/src/bin/perf.rs under sudo, or lower kernel.perf_event_paranoid");
                        } else {
                            eprintln!("hint: trig_probe failed before perf completed; the benchmark stderr above is the relevant failure.");
                        }
                        std::process::exit(1);
                    });

                    let cpi = match (row.cycles, row.instructions) {
                        (Some(c), Some(i)) if i > 0.0 => Some(c / i),
                        _ => None,
                    };

                    let ns_per_call = row.time_elapsed_s.map(|t| 1e9 * t / args.iters as f64);

                    println!(
                        "run {:2}: ins={} cyc={} CPI={} br={} miss={} ns/call={} energy(J)={} avg-power(W)={} power-source={}",
                        run,
                        fmt_i64(row.instructions),
                        fmt_i64(row.cycles),
                        fmt_f64(cpi),
                        fmt_i64(row.branches),
                        fmt_i64(row.branch_misses),
                        fmt_f64(ns_per_call),
                        fmt_f64(row.energy_joules),
                        fmt_f64(row.average_power_watts),
                        row.power_source.as_deref().unwrap_or("-"),
                    );

                    if args.require_power
                        && row.energy_joules.is_none()
                        && row.average_power_watts.is_none()
                    {
                        eprintln!(
                            "power collection failed for {op}/{mode}/{prec} run {run}: no readable energy PMU, sysfs sensor, or Pi PMIC source"
                        );
                        eprintln!(
                            "hint: grant access to the package energy counter or use --power-source pi-pmic on Raspberry Pi 5"
                        );
                        std::process::exit(1);
                    }

                    rows.push(row);
                    write_results_csv_or_exit(&args.out, &rows);
                }
            }
        }
    }

    write_results_csv_or_exit(&args.out, &rows);
    println!("\nWrote {}", args.out.display());
    let rows_with_power = rows
        .iter()
        .filter(|row| row.energy_joules.is_some() || row.average_power_watts.is_some())
        .count();
    println!(
        "CSV rows with power metrics: {rows_with_power}/{}",
        rows.len()
    );

    println!("\nSummary (mean (+/-) stdev):");
    for op in args.ops.iter().map(String::as_str) {
        for mode in args.modes.iter().map(String::as_str) {
            for prec in args.precisions.iter().map(String::as_str) {
                let subset: Vec<&Row> = rows
                    .iter()
                    .filter(|r| r.op == op && r.mode == mode && r.precision == prec)
                    .collect();

                if !subset.is_empty() {
                    print_summary(op, mode, prec, args.iters, &subset);
                }
            }
        }
    }

    println!("\nDirect comparison (fpx vs nat):");
    for op in args.ops.iter().map(String::as_str) {
        for prec in args.precisions.iter().map(String::as_str) {
            let fpx: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "fpx" && r.precision == prec)
                .collect();

            let nat: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "nat" && r.precision == prec)
                .collect();

            if !fpx.is_empty() && !nat.is_empty() {
                print_comparison(op, prec, args.iters, "fpx", &fpx, "nat", &nat);
            }
        }
    }

    println!("\nDirect comparison (fpx vs libm):");
    for op in args.ops.iter().map(String::as_str) {
        for prec in args.precisions.iter().map(String::as_str) {
            let fpx: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "fpx" && r.precision == prec)
                .collect();

            let libm: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "libm" && r.precision == prec)
                .collect();

            if !fpx.is_empty() && !libm.is_empty() {
                print_comparison(op, prec, args.iters, "fpx", &fpx, "libm", &libm);
            }
        }
    }

    println!("\nDirect comparison (fpx vs x87):");
    for op in args.ops.iter().map(String::as_str) {
        for prec in args.precisions.iter().map(String::as_str) {
            let fpx: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "fpx" && r.precision == prec)
                .collect();

            let x87: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "x87" && r.precision == prec)
                .collect();

            if !fpx.is_empty() && !x87.is_empty() {
                print_comparison(op, prec, args.iters, "fpx", &fpx, "x87", &x87);
            }
        }
    }

    println!("\nDirect comparison (kernel Horner vs kernel Estrin):");
    for op in args.ops.iter().map(String::as_str) {
        for prec in args.precisions.iter().map(String::as_str) {
            let horner: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "fpx_kernel_horner" && r.precision == prec)
                .collect();

            let estrin: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "fpx_kernel_estrin" && r.precision == prec)
                .collect();

            if !horner.is_empty() && !estrin.is_empty() {
                print_comparison(
                    op,
                    prec,
                    args.iters,
                    "fpx_kernel_horner",
                    &horner,
                    "fpx_kernel_estrin",
                    &estrin,
                );
            }
        }
    }

    println!("\nDirect comparison (full fpx vs kernel Horner):");
    for op in args.ops.iter().map(String::as_str) {
        for prec in args.precisions.iter().map(String::as_str) {
            let full: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "fpx" && r.precision == prec)
                .collect();

            let kernel: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "fpx_kernel_horner" && r.precision == prec)
                .collect();

            if !full.is_empty() && !kernel.is_empty() {
                print_comparison(
                    op,
                    prec,
                    args.iters,
                    "fpx",
                    &full,
                    "fpx_kernel_horner",
                    &kernel,
                );
            }
        }
    }

    println!("\nDirect comparison (nat vs x87):");
    for op in args.ops.iter().map(String::as_str) {
        for prec in args.precisions.iter().map(String::as_str) {
            let nat: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "nat" && r.precision == prec)
                .collect();

            let x87: Vec<&Row> = rows
                .iter()
                .filter(|r| r.op == op && r.mode == "x87" && r.precision == prec)
                .collect();

            if !nat.is_empty() && !x87.is_empty() {
                print_comparison(op, prec, args.iters, "nat", &nat, "x87", &x87);
            }
        }
    }
}
