// aether-fpx-bench/src/runner.rs
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::Instant;

use crate::perf::parse_perf_stat_csv;
use crate::types::Row;

pub fn run_perf(
    elf: &Path,
    op: &str,
    mode: &str,
    prec: &str,
    iters: usize,
    run: usize,
    events: &[String],
) -> Result<Row, String> {
    let events_joined = events.join(",");

    let mut cmd = Command::new("perf");
    cmd.env("LC_ALL", "C")
        .arg("stat")
        .arg("-x")
        .arg(",")
        .arg("--no-big-num")
        .arg("-e")
        .arg(events_joined)
        .arg("--")
        .arg(elf)
        .arg(op)
        .arg(mode)
        .arg(prec)
        .arg(iters.to_string())
        .stdout(Stdio::null())
        .stderr(Stdio::piped());

    let start = Instant::now();
    let out = cmd
        .output()
        .map_err(|e| format!("perf launch failed: {e}"))?;
    let elapsed = start.elapsed().as_secs_f64();

    if !out.status.success() {
        return Err(String::from_utf8_lossy(&out.stderr).to_string());
    }

    let parsed = parse_perf_stat_csv(&String::from_utf8_lossy(&out.stderr));

    Ok(Row {
        op: op.into(),
        mode: mode.into(),
        precision: prec.into(),
        run,
        iterations: iters,
        instructions: parsed.get("instructions").copied(),
        cycles: parsed.get("cycles").copied(),
        branches: parsed.get("branches").copied(),
        branch_misses: parsed.get("branch-misses").copied(),
        time_elapsed_s: Some(elapsed),
        energy_joules: None,
        average_power_watts: None,
        power_source: None,
    })
}

use std::fs;
use std::path::PathBuf;

/// Locate a release-built ELF for a binary name in a Cargo workspace.
///
/// Tries:
///   target/release/<bin>
///   target/release/deps/<bin>-<hash>
pub fn find_release_elf(bin_name: &str) -> Option<PathBuf> {
    let plain = PathBuf::from(format!("target/release/{bin_name}"));
    if plain.exists() {
        return Some(plain);
    }

    let deps = PathBuf::from("target/release/deps");
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
