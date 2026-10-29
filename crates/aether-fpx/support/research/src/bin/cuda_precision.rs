use std::env;
use std::fs;
use std::io::Write;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

const CUDA_MAKE_DIR: &str = "crates/aether-fpx/src/gpu/cuda";
const CUDA_PRECISION_TARGET: &str = "build/fpx_cuda_precision";
const RUNNER_SCRIPT: &str = "./fpx-runner.sh";
const CUDA_PRECISION_DEFAULT_GENCODE: &str =
    "-gencode arch=compute_90,code=sm_90 -gencode arch=compute_90a,code=sm_90a";

struct SubmitPlan {
    forwarded_args: Vec<String>,
    csv_rel: PathBuf,
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
        "usage: cargo cuda-precision -- [cuda-precision-options]\n\n\
This command builds the standalone CUDA precision benchmark, generates a batch\n\
job script, and submits it via fpx-runner.sh / run_gpu_hopper94.\n\n\
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
  --samples <n>\n\
  --block-size <n>\n\
  --csv <path>\n\n\
Examples:\n\
  cargo cuda-precision\n\
  cargo cuda-precision -- --only sin,cos,tan --precisions f16,f32,f64\n\
  cargo cuda-precision -- --modes fpx,nat --samples 262144 --csv target/cuda_precision_summary_sin.csv --only sin\n"
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
    let mut csv_rel = PathBuf::from(format!("target/cuda_precision_summary_{stamp}.csv"));
    let log_rel = PathBuf::from(format!("target/cuda_precision_job_{stamp}.log"));
    let job_script_rel = PathBuf::from(format!("cuda_precision_job_{stamp}.sh"));
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
                forwarded_args.push(normalized_args[index].clone());
                forwarded_args.push(normalized_args[index + 1].clone());
                index += 2;
            }
            _ => {
                forwarded_args.push(normalized_args[index].clone());
                index += 1;
            }
        }
    }

    if !forwarded_args.iter().any(|arg| arg == "--csv") {
        forwarded_args.push("--csv".into());
        forwarded_args.push(csv_rel.to_string_lossy().into_owned());
    }

    if !memory.chars().any(|ch| ch.is_ascii_alphabetic()) {
        return Err("--memory must include units, for example 1gb or 200mb".into());
    }

    if parse_walltime_seconds(&walltime)? < 3600 {
        return Err("--walltime must be at least 01:00:00 for run_gpu_hopper94".into());
    }

    Ok(SubmitPlan {
        forwarded_args,
        csv_rel,
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

fn build_cuda_precision(root: &Path) -> Result<PathBuf, String> {
    let make_dir = root.join(CUDA_MAKE_DIR);
    let mut command = Command::new("make");
    command.current_dir(&make_dir).arg("-B").arg("precision");

    if env::var_os("CUDA_GENCODE").is_none() {
        command.env("CUDA_GENCODE", CUDA_PRECISION_DEFAULT_GENCODE);
    }

    let status = command
        .status()
        .map_err(|err| format!("failed to invoke make precision: {err}"))?;

    if !status.success() {
        return Err(format!("make precision failed with status {status}"));
    }

    let bench_path = make_dir.join(CUDA_PRECISION_TARGET);
    if !bench_path.exists() {
        return Err(format!(
            "make precision reported success, but {} does not exist",
            bench_path.display()
        ));
    }

    Ok(bench_path)
}

fn write_job_script(root: &Path, bench_path: &Path, plan: &SubmitPlan) -> Result<(), String> {
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

    let script = format!(
        "#!/bin/bash\n\
set -uo pipefail\n\
source /apps/profiles/modules_asax.sh.dyn\n\
module load cuda/12.8.1\n\
cd {cwd}\n\
exec > {log} 2>&1\n\
echo \"[cuda-precision] started $(date -Is)\"\n\
echo \"[cuda-precision] working directory: $(pwd)\"\n\
echo \"[cuda-precision] benchmark command: {cmd}\"\n\
{cmd}\n\
status=$?\n\
echo \"[cuda-precision] finished $(date -Is) status=$status\"\n\
exit $status\n",
        cwd = shell_quote(&root.to_string_lossy()),
        log = shell_quote(&log_path.to_string_lossy()),
        cmd = benchmark_command,
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
    let root = workspace_root();

    let plan = parse_plan(&cli_args).unwrap_or_else(|err| {
        eprintln!("cuda_precision: {err}");
        std::process::exit(1);
    });

    let bench_path = build_cuda_precision(&root).unwrap_or_else(|err| {
        eprintln!("cuda_precision: {err}");
        std::process::exit(1);
    });

    write_job_script(&root, &bench_path, &plan).unwrap_or_else(|err| {
        eprintln!("cuda_precision: {err}");
        std::process::exit(1);
    });

    println!("Submitting CUDA precision batch job.");
    println!("  benchmark: {}", bench_path.display());
    println!("  job script: {}", plan.job_script_rel.display());
    println!("  log file: {}", plan.log_rel.display());
    println!("  precision csv: {}", plan.csv_rel.display());

    submit_job(&root, &plan).unwrap_or_else(|err| {
        eprintln!("cuda_precision: {err}");
        std::process::exit(1);
    });

    println!("Batch job submitted. Results will be written asynchronously.");
}
