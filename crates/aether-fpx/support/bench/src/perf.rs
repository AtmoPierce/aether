// aether-fpx-bench/src/perf.rs
use std::collections::HashMap;
use std::path::Path;
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;

use crate::pi_pmic::{start_pi_pmic_power, PiPmicRun};
use crate::sysfs_power::PowerRun as SysfsPowerRun;
use crate::sysfs_power::{PowerReading, SysfsPower};
use crate::types::Row;

const POWER_SAMPLE_INTERVAL: Duration = Duration::from_millis(20);
static REPORTED_MISSING_PERF_ENERGY: AtomicBool = AtomicBool::new(false);

enum FallbackPowerRun {
    Sysfs(SysfsPowerRun),
    PiPmic(PiPmicRun),
}

impl FallbackPowerRun {
    fn start() -> Option<Self> {
        SysfsPower::autodetect()
            .map(|source| Self::Sysfs(source.start()))
            .or_else(|| {
                start_pi_pmic_power(POWER_SAMPLE_INTERVAL)
                    .ok()
                    .map(Self::PiPmic)
            })
    }

    fn finish(self, elapsed: Duration) -> Option<PowerReading> {
        match self {
            Self::Sysfs(run) => run.finish(elapsed),
            Self::PiPmic(run) => run.finish(elapsed).ok(),
        }
    }
}

/// Run `perf stat` on a benchmark binary and return a populated Row.
pub fn run_perf_stat(
    elf: &Path,
    op: &str,
    mode: &str,
    prec: &str,
    iters: usize,
    run: usize,
    events: &[&str],
) -> Result<Row, String> {
    run_perf_stat_impl(elf, op, mode, prec, iters, run, events, true)
}

/// Run `perf stat` without probing or sampling fallback power sources.
///
/// This keeps short counter-only microbenchmarks from being dominated by
/// sysfs or PMIC discovery and sampling latency.
pub fn run_perf_stat_counters_only(
    elf: &Path,
    op: &str,
    mode: &str,
    prec: &str,
    iters: usize,
    run: usize,
    events: &[&str],
) -> Result<Row, String> {
    run_perf_stat_impl(elf, op, mode, prec, iters, run, events, false)
}

fn run_perf_stat_impl(
    elf: &Path,
    op: &str,
    mode: &str,
    prec: &str,
    iters: usize,
    run: usize,
    events: &[&str],
    collect_fallback_power: bool,
) -> Result<Row, String> {
    let (out, elapsed, sysfs_reading) =
        match run_perf_command(elf, op, mode, prec, iters, events, collect_fallback_power) {
            Ok(result) => result,
            Err(error) => return Err(error),
        };

    let mut out = out;
    let mut elapsed = elapsed;
    let mut sysfs_reading = sysfs_reading;
    if !out.status.success()
        && should_retry_without_energy(events, &String::from_utf8_lossy(&out.stderr))
    {
        let retry_events = events
            .iter()
            .copied()
            .filter(|event| !is_energy_event(event))
            .collect::<Vec<_>>();
        let retry = run_perf_command(
            elf,
            op,
            mode,
            prec,
            iters,
            &retry_events,
            collect_fallback_power,
        )?;
        out = retry.0;
        elapsed = retry.1;
        sysfs_reading = retry.2;
    }

    if !out.status.success() {
        return Err(String::from_utf8_lossy(&out.stderr).into());
    }

    let elapsed_s = elapsed.as_secs_f64();
    let parsed = parse_perf_stat_csv(&String::from_utf8_lossy(&out.stderr));
    let measured_elapsed_s = parsed
        .get("seconds time elapsed")
        .copied()
        .filter(|value| *value > 0.0)
        .unwrap_or(elapsed_s);
    let energy_joules = parsed
        .iter()
        .find_map(|(event, value)| is_energy_event(event).then_some(*value));
    if energy_joules.is_none() && events.iter().any(|event| is_energy_event(event)) {
        report_missing_perf_energy_once(events, &String::from_utf8_lossy(&out.stderr));
    }
    let average_power_watts = energy_joules
        .filter(|_| measured_elapsed_s > 0.0)
        .map(|energy| energy / measured_elapsed_s);
    let power_reading =
        power_reading_from_perf(&parsed, energy_joules, average_power_watts).or(sysfs_reading);

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

        time_elapsed_s: Some(measured_elapsed_s),
        energy_joules: power_reading
            .as_ref()
            .and_then(|reading| reading.energy_joules),
        average_power_watts: power_reading
            .as_ref()
            .and_then(|reading| reading.average_power_watts),
        power_source: power_reading.map(|reading| reading.source),
    })
}

fn run_perf_command(
    elf: &Path,
    op: &str,
    mode: &str,
    prec: &str,
    iters: usize,
    events: &[&str],
    collect_fallback_power: bool,
) -> Result<(Output, Duration, Option<PowerReading>), String> {
    let mut cmd = Command::new("perf");
    cmd.env("LC_ALL", "C")
        .arg("stat")
        .arg("-x")
        .arg(",")
        .arg("--no-big-num");
    for event in events {
        cmd.arg("-e").arg(event);
    }
    cmd.arg("--")
        .arg(elf)
        .arg(op)
        .arg(mode)
        .arg(prec)
        .arg(iters.to_string())
        .stdout(Stdio::null())
        .stderr(Stdio::piped());

    let fallback_power_run = collect_fallback_power
        .then(FallbackPowerRun::start)
        .flatten();

    let start = std::time::Instant::now();
    let out = cmd.output().map_err(|e| e.to_string())?;
    let elapsed = start.elapsed();
    let fallback_power_reading = fallback_power_run.and_then(|run| run.finish(elapsed));

    Ok((out, elapsed, fallback_power_reading))
}

fn is_energy_event(event: &str) -> bool {
    event.to_ascii_lowercase().contains("energy")
}

fn should_retry_without_energy(events: &[&str], stderr: &str) -> bool {
    if !events.iter().any(|event| is_energy_event(event)) {
        return false;
    }

    let stderr = stderr.to_ascii_lowercase();
    stderr.contains("energy")
        || stderr.contains("no permission to enable")
        || stderr.contains("permission denied")
        || stderr.contains("not supported")
        || stderr.contains("event syntax error")
}

fn power_reading_from_perf(
    parsed: &HashMap<String, f64>,
    energy_joules: Option<f64>,
    average_power_watts: Option<f64>,
) -> Option<PowerReading> {
    let energy_joules = energy_joules?;
    Some(PowerReading {
        source: parsed
            .keys()
            .find(|event| is_energy_event(event))
            .map(|event| format!("perf:{event}"))
            .unwrap_or_else(|| "perf:energy".to_string()),
        energy_joules: Some(energy_joules),
        average_power_watts,
    })
}

fn report_missing_perf_energy_once(events: &[&str], stderr: &str) {
    if REPORTED_MISSING_PERF_ENERGY.swap(true, Ordering::Relaxed) {
        return;
    }

    let requested = events
        .iter()
        .copied()
        .find(|event| is_energy_event(event))
        .unwrap_or("energy");
    if let Some(line) = perf_energy_diagnostic_line(stderr, requested) {
        eprintln!(
            "warning: requested perf energy event `{requested}` did not produce a numeric Joules value: {line}"
        );
    } else {
        eprintln!(
            "warning: requested perf energy event `{requested}` did not appear in perf stat output"
        );
    }
}

fn perf_energy_diagnostic_line<'a>(stderr: &'a str, requested: &str) -> Option<&'a str> {
    stderr
        .lines()
        .map(str::trim)
        .find(|line| line.contains(requested) || line.to_ascii_lowercase().contains("energy"))
}

/// Parse `perf stat -x,` CSV output (stderr).
pub fn parse_perf_stat_csv(text: &str) -> std::collections::HashMap<String, f64> {
    let mut m = std::collections::HashMap::new();

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let mut it = line.splitn(3, ',');
        let raw = it.next().unwrap_or("").trim();
        let unit = it.next().unwrap_or("").trim();
        let rest = it.next().unwrap_or("").trim();

        if raw.is_empty() || raw.starts_with('<') {
            continue;
        }

        let Ok(v) = raw.replace([' ', '\t'], "").parse::<f64>() else {
            continue;
        };

        let event = rest
            .split(',')
            .next()
            .filter(|event| !event.trim().is_empty())
            .unwrap_or(unit)
            .trim();
        if !event.is_empty() {
            m.insert(event.to_string(), v);
        }
    }

    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_energy_event_with_joule_unit() {
        let parsed = parse_perf_stat_csv(
            "0.125000,Joules,power/energy-pkg/,1,100.00,,\n0.010000,,seconds time elapsed\n",
        );
        assert_eq!(parsed.get("power/energy-pkg/"), Some(&0.125));
        assert_eq!(parsed.get("seconds time elapsed"), Some(&0.01));
        assert!(is_energy_event("power/energy-pkg/"));
    }

    #[test]
    fn parses_elapsed_time_when_perf_omits_empty_unit_column() {
        let parsed = parse_perf_stat_csv("0.010000,seconds time elapsed\n");
        assert_eq!(parsed.get("seconds time elapsed"), Some(&0.01));
    }

    #[test]
    fn retries_failed_energy_events_without_losing_perf_counters() {
        assert!(should_retry_without_energy(
            &["instructions", "power/energy-pkg/"],
            "No permission to enable power/energy-pkg/ event."
        ));
        assert!(!should_retry_without_energy(
            &["instructions", "cycles"],
            "No permission to enable cycles event."
        ));
    }

    #[test]
    fn reports_non_numeric_energy_lines() {
        assert_eq!(
            perf_energy_diagnostic_line(
                "<not counted>,Joules,power/energy-pkg/,0,0.00,,",
                "power/energy-pkg/"
            ),
            Some("<not counted>,Joules,power/energy-pkg/,0,0.00,,")
        );
    }
}
