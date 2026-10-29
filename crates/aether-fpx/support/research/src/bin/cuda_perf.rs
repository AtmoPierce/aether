use std::env;
use std::fs;
use std::io::Write;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

use aether_fpx_bench::{export::write_csv, types::Row};

const CUDA_MAKE_DIR: &str = "crates/aether-fpx/src/gpu/cuda";
const CUDA_BENCH_TARGET: &str = "build/fpx_cuda_perf";
const RUNNER_SCRIPT: &str = "./fpx-runner.sh";
const CUDA_PERF_DEFAULT_GENCODE: &str =
    "-gencode arch=compute_90,code=sm_90 -gencode arch=compute_90a,code=sm_90a";
const EXPORT_COMPATIBLE_CSV_FLAG: &str = "--export-compatible-csv";
const RAW_CUDA_CSV_HEADER: &str =
    "op,mode,precision,run,elements,launches,time_elapsed_s,ns_per_launch,ns_per_element,elements_per_second,sample_output,average_power_watts,power_source";
const DEFAULT_POWER_SAMPLE_MS: &str = "200";

struct SubmitPlan {
    forwarded_args: Vec<String>,
    csv_rel: PathBuf,
    raw_csv_rel: PathBuf,
    log_rel: PathBuf,
    job_script_rel: PathBuf,
    queue: String,
    cpus: String,
    memory: String,
    walltime: String,
    gpus: String,
    job_name: Option<String>,
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(4)
        .expect("aether-fpx-research should live under the workspace root")
        .to_path_buf()
}

fn timestamp_tag() -> String {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time should be after the Unix epoch");
    format!("{}_{}", now.as_secs(), now.subsec_millis())
}

fn shell_quote(value: &str) -> String {
    if value.is_empty() {
        return "''".to_string();
    }
    format!("'{}'", value.replace('\'', "'\"'\"'"))
}

fn raw_csv_path_for(csv_path: &Path) -> PathBuf {
    let stem = csv_path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("cuda_perf");
    let file_name = match csv_path.extension().and_then(|ext| ext.to_str()) {
        Some(ext) if !ext.is_empty() => format!("{stem}_raw.{ext}"),
        _ => format!("{stem}_raw"),
    };

    match csv_path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        Some(parent) => parent.join(file_name),
        None => PathBuf::from(file_name),
    }
}

fn prepare_output_path(out: &Path) -> Result<(), String> {
    if let Some(parent) = out.parent().filter(|parent| !parent.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|err| {
            format!(
                "failed to create CSV output directory {}: {err}",
                parent.display()
            )
        })?;
    }
    Ok(())
}

fn parse_usize(value: &str, field: &str, line_number: usize, path: &Path) -> Result<usize, String> {
    value.parse::<usize>().map_err(|_| {
        format!(
            "invalid {field} at line {line_number} in {}: {value}",
            path.display()
        )
    })
}

fn parse_f64(value: &str, field: &str, line_number: usize, path: &Path) -> Result<f64, String> {
    value.parse::<f64>().map_err(|_| {
        format!(
            "invalid {field} at line {line_number} in {}: {value}",
            path.display()
        )
    })
}

fn parse_optional_f64(
    value: &str,
    field: &str,
    line_number: usize,
    path: &Path,
) -> Result<Option<f64>, String> {
    if value.trim().is_empty() {
        return Ok(None);
    }
    let parsed = parse_f64(value, field, line_number, path)?;
    Ok(parsed.is_finite().then_some(parsed))
}

fn parse_raw_csv_rows(text: &str, path: &Path) -> Result<Vec<Row>, String> {
    let mut lines = text.lines();
    let header = lines.next().unwrap_or_default().trim();
    if header != RAW_CUDA_CSV_HEADER {
        return Err(format!(
            "unexpected CUDA raw CSV header in {}: {header}",
            path.display()
        ));
    }

    let mut rows = Vec::new();
    for (line_index, line) in lines.enumerate() {
        let line_number = line_index + 2;
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() != 13 {
            return Err(format!(
                "expected 13 CSV fields at line {line_number} in {}, found {}",
                path.display(),
                fields.len()
            ));
        }

        let elements = parse_usize(fields[4], "elements", line_number, path)?;
        let launches = parse_usize(fields[5], "launches", line_number, path)?;
        let iterations = elements.checked_mul(launches).ok_or_else(|| {
            format!(
                "iterations overflow at line {line_number} in {}: {elements} * {launches}",
                path.display()
            )
        })?;
        let time_elapsed_s = parse_f64(fields[6], "time_elapsed_s", line_number, path)?;
        let average_power_watts =
            parse_optional_f64(fields[11], "average_power_watts", line_number, path)?;

        rows.push(Row {
            op: fields[0].to_string(),
            mode: fields[1].to_string(),
            precision: fields[2].to_string(),
            run: parse_usize(fields[3], "run", line_number, path)?,
            iterations,
            instructions: None,
            cycles: None,
            branches: None,
            branch_misses: None,
            time_elapsed_s: Some(time_elapsed_s),
            energy_joules: average_power_watts.map(|watts| watts * time_elapsed_s),
            average_power_watts,
            power_source: if average_power_watts.is_some() && !fields[12].trim().is_empty() {
                Some(fields[12].trim().to_string())
            } else {
                None
            },
        });
    }

    if rows.is_empty() {
        return Err(format!(
            "{} did not contain any benchmark rows",
            path.display()
        ));
    }

    Ok(rows)
}

fn load_compatible_rows(raw_csv: &Path) -> Result<Vec<Row>, String> {
    let text = fs::read_to_string(raw_csv)
        .map_err(|err| format!("failed to read {}: {err}", raw_csv.display()))?;
    parse_raw_csv_rows(&text, raw_csv)
}

fn export_compatible_csv(raw_csv: &Path, out: &Path) -> Result<(), String> {
    let rows = load_compatible_rows(raw_csv)?;
    prepare_output_path(out)?;
    write_csv(out, &rows)
        .map_err(|err| format!("failed to write CSV output {}: {err}", out.display()))
}

struct ExportArgs {
    raw_csv: PathBuf,
    out: PathBuf,
}

fn parse_export_args(cli_args: &[String]) -> Result<ExportArgs, String> {
    let mut raw_csv = None;
    let mut out = None;
    let mut index = 0usize;

    while index < cli_args.len() {
        match cli_args[index].as_str() {
            "--raw-csv" => {
                if index + 1 >= cli_args.len() {
                    return Err("missing value for --raw-csv".into());
                }
                raw_csv = Some(PathBuf::from(&cli_args[index + 1]));
                index += 2;
            }
            "--out" => {
                if index + 1 >= cli_args.len() {
                    return Err("missing value for --out".into());
                }
                out = Some(PathBuf::from(&cli_args[index + 1]));
                index += 2;
            }
            other => return Err(format!("unknown argument for exporter mode: {other}")),
        }
    }

    Ok(ExportArgs {
        raw_csv: raw_csv.ok_or("--raw-csv is required in exporter mode")?,
        out: out.ok_or("--out is required in exporter mode")?,
    })
}

fn maybe_run_export_mode(cli_args: &[String]) -> Result<bool, String> {
    if cli_args.first().map(String::as_str) != Some(EXPORT_COMPATIBLE_CSV_FLAG) {
        return Ok(false);
    }

    let args = parse_export_args(&cli_args[1..])?;
    export_compatible_csv(&args.raw_csv, &args.out)?;
    println!(
        "Wrote compatible CSV {} from raw CUDA CSV {}",
        args.out.display(),
        args.raw_csv.display()
    );
    Ok(true)
}

fn parse_walltime_seconds(value: &str) -> Result<u64, String> {
    let mut parts = value.split(':');
    let hours = parts
        .next()
        .ok_or_else(|| "--walltime must use HH:MM:SS".to_string())?
        .parse::<u64>()
        .map_err(|_| "--walltime must use HH:MM:SS".to_string())?;
    let minutes = parts
        .next()
        .ok_or_else(|| "--walltime must use HH:MM:SS".to_string())?
        .parse::<u64>()
        .map_err(|_| "--walltime must use HH:MM:SS".to_string())?;
    let seconds = parts
        .next()
        .ok_or_else(|| "--walltime must use HH:MM:SS".to_string())?
        .parse::<u64>()
        .map_err(|_| "--walltime must use HH:MM:SS".to_string())?;

    if parts.next().is_some() || minutes >= 60 || seconds >= 60 {
        return Err("--walltime must use HH:MM:SS".into());
    }

    Ok(hours * 3600 + minutes * 60 + seconds)
}

fn looks_like_job_id(line: &str) -> bool {
    let trimmed = line.trim();
    let mut parts = trimmed.split('.');
    let job_number = parts.next().unwrap_or("");
    let suffix = parts.next().unwrap_or("");

    !job_number.is_empty()
        && !suffix.is_empty()
        && parts.next().is_none()
        && job_number.chars().all(|ch| ch.is_ascii_digit())
}

fn print_usage() {
    println!(
        "usage: cargo cuda-perf -- [cuda-benchmark-options]\n\n\
This command builds the standalone CUDA benchmark, generates a batch job script,\n\
and submits it via fpx-runner.sh / run_gpu_hopper94. The job writes a CPU-compatible\n\
CSV at --csv, keeps the raw CUDA metrics in a sibling *_raw.csv sidecar, and fills\n\
the optional power columns from nvidia-smi when it is available on the GPU node.\n\n\
Submission options handled by this wrapper:\n\
    --queue <name>         default: gpu\n\
    --cpus <n>             default: 1\n\
    --memory <amount>      default: 1gb\n\
    --walltime <HH:MM:SS>  default: 01:00:00\n\
    --gpus <n>             default: 1\n\
    --job-name <name>      default: helper-generated name\n\n\
Benchmark options forwarded to the CUDA executable:\n\
  --only <op1,op2,...>\n\
    --modes <m1,m2,...>\n\
  --precisions <p1,p2,...>\n\
  --runs <n>\n\
  --launches <n>\n\
  --elements <n>\n\
  --block-size <n>\n\
    --csv <path>                  Output CPU-compatible CSV\n\n\
Examples:\n\
  cargo cuda-perf\n\
    cargo cuda-perf -- --memory 4gb --walltime 01:00:00\n\
  cargo cuda-perf -- --only sin,cos,tan --precisions f32,f64 --runs 5\n\
                cargo cuda-perf -- --memory 4gb --walltime 01:00:00 --csv target/cuda_perf_sin.csv --only sin --precisions f16,fp8_e4m3\n"
    );
}

fn normalize_forwarded_args(cli_args: &[String]) -> Result<Vec<String>, String> {
    let mut normalized = Vec::new();
    let mut index = 0usize;

    while index < cli_args.len() {
        let arg = &cli_args[index];
        if arg == "--only" || arg == "--modes" || arg == "--precisions" {
            normalized.push(arg.clone());
            index += 1;

            if index >= cli_args.len() || cli_args[index].starts_with('-') {
                return Err(format!("missing value for {arg}"));
            }

            let mut value = String::new();
            while index < cli_args.len() && !cli_args[index].starts_with('-') {
                value.push_str(&cli_args[index]);
                index += 1;
            }
            normalized.push(value);
            continue;
        }

        normalized.push(arg.clone());
        index += 1;
    }

    Ok(normalized)
}

fn parse_plan(cli_args: &[String]) -> Result<SubmitPlan, String> {
    if cli_args.iter().any(|arg| arg == "--help" || arg == "-h") {
        print_usage();
        std::process::exit(0);
    }

    let stamp = timestamp_tag();
    let normalized_args = normalize_forwarded_args(cli_args)?;
    let mut csv_rel = PathBuf::from(format!("target/cuda_perf_{stamp}.csv"));
    let log_rel = PathBuf::from(format!("target/cuda_perf_job_{stamp}.log"));
    let job_script_rel = PathBuf::from(format!("cuda_perf_job_{stamp}.sh"));
    let mut queue = "gpu".to_string();
    let mut cpus = "1".to_string();
    let mut memory = "1gb".to_string();
    let mut walltime = "01:00:00".to_string();
    let mut gpus = "1".to_string();
    let mut job_name = None;
    let mut forwarded_args = Vec::new();

    let mut index = 0usize;
    while index < normalized_args.len() {
        match normalized_args[index].as_str() {
            "--queue" => {
                if index + 1 >= normalized_args.len() {
                    return Err("missing value for --queue".into());
                }
                queue = normalized_args[index + 1].clone();
                index += 2;
            }
            "--cpus" => {
                if index + 1 >= normalized_args.len() {
                    return Err("missing value for --cpus".into());
                }
                cpus = normalized_args[index + 1].clone();
                index += 2;
            }
            "--memory" => {
                if index + 1 >= normalized_args.len() {
                    return Err("missing value for --memory".into());
                }
                memory = normalized_args[index + 1].clone();
                index += 2;
            }
            "--walltime" => {
                if index + 1 >= normalized_args.len() {
                    return Err("missing value for --walltime".into());
                }
                walltime = normalized_args[index + 1].clone();
                index += 2;
            }
            "--gpus" => {
                if index + 1 >= normalized_args.len() {
                    return Err("missing value for --gpus".into());
                }
                gpus = normalized_args[index + 1].clone();
                index += 2;
            }
            "--job-name" => {
                if index + 1 >= normalized_args.len() {
                    return Err("missing value for --job-name".into());
                }
                job_name = Some(normalized_args[index + 1].clone());
                index += 2;
            }
            "--csv" => {
                if index + 1 >= normalized_args.len() {
                    return Err("missing value for --csv".into());
                }
                csv_rel = PathBuf::from(&normalized_args[index + 1]);
                index += 2;
            }
            _ => {
                forwarded_args.push(normalized_args[index].clone());
                index += 1;
            }
        }
    }

    let raw_csv_rel = raw_csv_path_for(&csv_rel);
    if !normalized_args.iter().any(|arg| arg == "--power-sample-ms") {
        forwarded_args.push("--power-sample-ms".into());
        forwarded_args.push(DEFAULT_POWER_SAMPLE_MS.into());
    }
    forwarded_args.push("--csv".into());
    forwarded_args.push(raw_csv_rel.to_string_lossy().into_owned());

    if !memory.chars().any(|ch| ch.is_ascii_alphabetic()) {
        return Err("--memory must include units, for example 1gb or 200mb".into());
    }

    if parse_walltime_seconds(&walltime)? < 3600 {
        return Err("--walltime must be at least 01:00:00 for run_gpu_hopper94".into());
    }

    Ok(SubmitPlan {
        forwarded_args,
        csv_rel,
        raw_csv_rel,
        log_rel,
        job_script_rel,
        queue,
        cpus,
        memory,
        walltime,
        gpus,
        job_name,
    })
}

fn build_cuda_bench(root: &Path) -> Result<PathBuf, String> {
    let make_dir = root.join(CUDA_MAKE_DIR);
    let mut command = Command::new("make");
    command.current_dir(&make_dir).arg("-B").arg("bench");

    if env::var_os("CUDA_GENCODE").is_none() {
        command.env("CUDA_GENCODE", CUDA_PERF_DEFAULT_GENCODE);
    }

    let status = command
        .status()
        .map_err(|err| format!("failed to invoke make bench: {err}"))?;

    if !status.success() {
        return Err(format!("make bench failed with status {status}"));
    }

    let bench_path = make_dir.join(CUDA_BENCH_TARGET);
    if !bench_path.exists() {
        return Err(format!(
            "make bench reported success, but {} does not exist",
            bench_path.display()
        ));
    }

    Ok(bench_path)
}

fn write_job_script(
    root: &Path,
    bench_path: &Path,
    converter_path: &Path,
    plan: &SubmitPlan,
) -> Result<(), String> {
    let script_path = root.join(&plan.job_script_rel);
    let log_path = root.join(&plan.log_rel);

    if let Some(parent) = script_path.parent() {
        fs::create_dir_all(parent).map_err(|err| {
            format!(
                "failed to create script directory {}: {err}",
                parent.display()
            )
        })?;
    }
    if let Some(parent) = log_path.parent() {
        fs::create_dir_all(parent)
            .map_err(|err| format!("failed to create log directory {}: {err}", parent.display()))?;
    }

    let mut command_parts = Vec::new();
    command_parts.push(shell_quote(&bench_path.to_string_lossy()));
    for arg in &plan.forwarded_args {
        command_parts.push(shell_quote(arg));
    }
    let benchmark_command = command_parts.join(" ");
    let export_command = format!(
        "{} {} --raw-csv {} --out {}",
        shell_quote(&converter_path.to_string_lossy()),
        shell_quote(EXPORT_COMPATIBLE_CSV_FLAG),
        shell_quote(&plan.raw_csv_rel.to_string_lossy()),
        shell_quote(&plan.csv_rel.to_string_lossy()),
    );

    let script = format!(
        "#!/bin/bash\n\
set -uo pipefail\n\
source /apps/profiles/modules_asax.sh.dyn\n\
module load cuda/12.8.1\n\
cd {cwd}\n\
exec > {log} 2>&1\n\
echo \"[cuda-perf] started $(date -Is)\"\n\
echo \"[cuda-perf] working directory: $(pwd)\"\n\
echo \"[cuda-perf] benchmark command: {cmd}\"\n\
{cmd}\n\
status=$?\n\
if [[ $status -eq 0 ]]; then\n\
echo \"[cuda-perf] export command: {export_cmd}\"\n\
{export_cmd}\n\
status=$?\n\
fi\n\
echo \"[cuda-perf] finished $(date -Is) status=$status\"\n\
exit $status\n",
        cwd = shell_quote(&root.to_string_lossy()),
        log = shell_quote(&log_path.to_string_lossy()),
        cmd = benchmark_command,
        export_cmd = export_command,
    );

    fs::write(&script_path, script)
        .map_err(|err| format!("failed to write {}: {err}", script_path.display()))?;

    let mut perms = fs::metadata(&script_path)
        .map_err(|err| format!("failed to stat {}: {err}", script_path.display()))?
        .permissions();
    perms.set_mode(0o755);
    fs::set_permissions(&script_path, perms)
        .map_err(|err| format!("failed to chmod {}: {err}", script_path.display()))?;

    Ok(())
}

fn submit_job(root: &Path, plan: &SubmitPlan) -> Result<(), String> {
    let mut child = Command::new("bash")
        .current_dir(root)
        .arg(RUNNER_SCRIPT)
        .arg(&plan.job_script_rel)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|err| format!("failed to invoke {RUNNER_SCRIPT}: {err}"))?;

    let answers = format!(
        "{}\n{}\n{}\n{}\n{}\n{}\n",
        plan.queue,
        plan.cpus,
        plan.walltime,
        plan.memory,
        plan.gpus,
        plan.job_name.as_deref().unwrap_or("")
    );

    child
        .stdin
        .as_mut()
        .ok_or_else(|| format!("failed to open stdin for {RUNNER_SCRIPT}"))?
        .write_all(answers.as_bytes())
        .map_err(|err| format!("failed to provide run details to {RUNNER_SCRIPT}: {err}"))?;

    let output = child
        .wait_with_output()
        .map_err(|err| format!("failed while waiting on {RUNNER_SCRIPT}: {err}"))?;

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    print!("{}", stdout);
    eprint!("{}", stderr);

    if !output.status.success() {
        return Err(format!(
            "{} failed with status {}",
            RUNNER_SCRIPT, output.status
        ));
    }

    if !stdout.lines().any(looks_like_job_id) {
        return Err(format!(
            "{} did not report a PBS job id, so the batch submission likely failed",
            RUNNER_SCRIPT
        ));
    }

    Ok(())
}

fn main() {
    let cli_args: Vec<String> = env::args().skip(1).collect();

    if maybe_run_export_mode(&cli_args).unwrap_or_else(|err| {
        eprintln!("cuda_perf: {err}");
        std::process::exit(1);
    }) {
        return;
    }

    let root = workspace_root();

    let plan = parse_plan(&cli_args).unwrap_or_else(|err| {
        eprintln!("cuda_perf: {err}");
        std::process::exit(1);
    });

    let bench_path = build_cuda_bench(&root).unwrap_or_else(|err| {
        eprintln!("cuda_perf: {err}");
        std::process::exit(1);
    });

    let converter_path = env::current_exe().unwrap_or_else(|err| {
        eprintln!("cuda_perf: failed to locate exporter executable: {err}");
        std::process::exit(1);
    });

    write_job_script(&root, &bench_path, &converter_path, &plan).unwrap_or_else(|err| {
        eprintln!("cuda_perf: {err}");
        std::process::exit(1);
    });

    println!("Submitting CUDA benchmark batch job.");
    println!("  benchmark: {}", bench_path.display());
    println!("  job script: {}", plan.job_script_rel.display());
    println!("  log file: {}", plan.log_rel.display());
    println!("  perf csv: {}", plan.csv_rel.display());
    println!("  raw perf csv: {}", plan.raw_csv_rel.display());
    println!(
        "  report hint: cargo cuda-perf-report -- --csv {}",
        plan.raw_csv_rel.display()
    );

    submit_job(&root, &plan).unwrap_or_else(|err| {
        eprintln!("cuda_perf: {err}");
        std::process::exit(1);
    });

    println!("Batch job submitted. Results will be written asynchronously.");
}

#[cfg(test)]
mod tests {
    use super::{parse_raw_csv_rows, raw_csv_path_for};
    use std::path::{Path, PathBuf};

    #[test]
    fn derives_raw_sidecar_path_next_to_requested_csv() {
        assert_eq!(
            raw_csv_path_for(Path::new("target/cuda_perf.csv")),
            PathBuf::from("target/cuda_perf_raw.csv")
        );
    }

    #[test]
    fn converts_raw_cuda_rows_to_the_shared_perf_schema() {
        let text = concat!(
            "op,mode,precision,run,elements,launches,time_elapsed_s,ns_per_launch,ns_per_element,elements_per_second,sample_output,average_power_watts,power_source\n",
            "sin,nat,f32,2,8,4,0.25,1.0,2.0,3.0,4.0,250.0,nvidia-smi:power.draw\n",
        );

        let rows = parse_raw_csv_rows(text, Path::new("target/cuda_perf_raw.csv")).unwrap();
        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        assert_eq!(row.op, "sin");
        assert_eq!(row.mode, "nat");
        assert_eq!(row.precision, "f32");
        assert_eq!(row.run, 2);
        assert_eq!(row.iterations, 32);
        assert_eq!(row.time_elapsed_s, Some(0.25));
        assert_eq!(row.instructions, None);
        assert_eq!(row.average_power_watts, Some(250.0));
        assert_eq!(row.energy_joules, Some(62.5));
        assert_eq!(row.power_source.as_deref(), Some("nvidia-smi:power.draw"));
    }
}
