use std::collections::HashMap;
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

use crate::sysfs_power::PowerReading;
use crate::types::Row;

pub const SOURCE: &str = "rpi5-pmic-rails";

pub struct PiPmicRun {
    initial_watts: f64,
    stop: Arc<AtomicBool>,
    sampler: thread::JoinHandle<Vec<(Duration, f64)>>,
}

/// Sum the output power of every PMIC rail with matching current and voltage
/// readings in `vcgencmd pmic_read_adc` output.
pub fn parse_pmic_rail_power_watts(text: &str) -> Result<f64, String> {
    let mut currents = HashMap::new();
    let mut voltages = HashMap::new();

    for line in text.lines() {
        let Some((left, right)) = line.trim().split_once('=') else {
            continue;
        };
        let Some(name) = left.split_whitespace().next() else {
            continue;
        };
        let Some(value) = right
            .trim()
            .strip_suffix(['A', 'V'])
            .and_then(|raw| raw.parse::<f64>().ok())
        else {
            continue;
        };

        if let Some(rail) = name.strip_suffix("_A") {
            currents.insert(rail.to_string(), value);
        } else if let Some(rail) = name.strip_suffix("_V") {
            voltages.insert(rail.to_string(), value);
        }
    }

    let mut matched = 0usize;
    let watts = currents
        .iter()
        .filter_map(|(rail, current)| {
            voltages.get(rail).map(|voltage| {
                matched += 1;
                current * voltage
            })
        })
        .sum::<f64>();

    if matched == 0 || !watts.is_finite() {
        Err("vcgencmd returned no matching PMIC current/voltage rail pairs".into())
    } else {
        Ok(watts)
    }
}

pub fn read_pmic_rail_power_watts() -> Result<f64, String> {
    let output = Command::new("vcgencmd")
        .arg("pmic_read_adc")
        .output()
        .map_err(|error| format!("failed to run vcgencmd pmic_read_adc: {error}"))?;
    if !output.status.success() {
        return Err(format!(
            "vcgencmd pmic_read_adc failed: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    parse_pmic_rail_power_watts(&String::from_utf8_lossy(&output.stdout))
}

pub fn start_pi_pmic_power(sample_interval: Duration) -> Result<PiPmicRun, String> {
    let initial_watts = read_pmic_rail_power_watts()?;
    let stop = Arc::new(AtomicBool::new(false));
    let sampler_stop = Arc::clone(&stop);
    let sampler = thread::spawn(move || {
        let start = Instant::now();
        let mut samples = Vec::new();
        while !sampler_stop.load(Ordering::Relaxed) {
            if let Ok(watts) = read_pmic_rail_power_watts() {
                samples.push((start.elapsed(), watts));
            }
            thread::sleep(sample_interval);
        }
        samples
    });

    Ok(PiPmicRun {
        initial_watts,
        stop,
        sampler,
    })
}

impl PiPmicRun {
    pub fn finish(self, elapsed: Duration) -> Result<PowerReading, String> {
        self.stop.store(true, Ordering::Relaxed);
        let samples = self
            .sampler
            .join()
            .map_err(|_| "PMIC sampling thread panicked".to_string())?;
        let final_watts = read_pmic_rail_power_watts()?;
        let elapsed_s = elapsed.as_secs_f64();
        let mut within_run = vec![(0.0, self.initial_watts)];
        within_run.extend(
            samples
                .into_iter()
                .map(|(time, watts)| (time.as_secs_f64(), watts))
                .filter(|(time, _)| *time > 0.0 && *time < elapsed_s),
        );
        within_run.push((elapsed_s, final_watts));
        within_run.sort_by(|lhs, rhs| lhs.0.total_cmp(&rhs.0));

        let energy_joules = within_run
            .windows(2)
            .map(|pair| 0.5 * (pair[0].1 + pair[1].1) * (pair[1].0 - pair[0].0))
            .sum::<f64>();
        let average_power_watts = (elapsed_s > 0.0).then_some(energy_joules / elapsed_s);

        Ok(PowerReading {
            source: SOURCE.into(),
            energy_joules: Some(energy_joules),
            average_power_watts,
        })
    }
}

pub fn run_pi_pmic_power(
    elf: &Path,
    op: &str,
    mode: &str,
    prec: &str,
    iters: usize,
    run: usize,
    sample_interval: Duration,
) -> Result<Row, String> {
    let power_run = start_pi_pmic_power(sample_interval)?;
    let start = Instant::now();
    let output = Command::new(elf)
        .arg(op)
        .arg(mode)
        .arg(prec)
        .arg(iters.to_string())
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .output()
        .map_err(|error| format!("benchmark launch failed: {error}"));
    let elapsed = start.elapsed();
    let power = power_run.finish(elapsed)?;
    let output = output?;
    if !output.status.success() {
        return Err(String::from_utf8_lossy(&output.stderr).into_owned());
    }

    let elapsed_s = elapsed.as_secs_f64();

    Ok(Row {
        op: op.into(),
        mode: mode.into(),
        precision: prec.into(),
        run,
        iterations: iters,
        instructions: None,
        cycles: None,
        branches: None,
        branch_misses: None,
        time_elapsed_s: Some(elapsed_s),
        energy_joules: power.energy_joules,
        average_power_watts: power.average_power_watts,
        power_source: Some(power.source),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sums_only_matched_pmic_rails() {
        let text = "VDD_CORE_A current(7)=2.00000000A\n\
                    3V3_SYS_A current(1)=0.50000000A\n\
                    VDD_CORE_V volt(15)=0.80000000V\n\
                    3V3_SYS_V volt(9)=3.30000000V\n\
                    EXT5V_V volt(24)=5.10000000V\n";
        let watts = parse_pmic_rail_power_watts(text).unwrap();
        assert!((watts - 3.25).abs() < 1e-12);
    }
}
