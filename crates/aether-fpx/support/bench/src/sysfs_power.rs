use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

const SAMPLE_INTERVAL: Duration = Duration::from_millis(20);

#[derive(Clone, Debug)]
pub enum SysfsPower {
    Energy { source: String, energy_uj: PathBuf },
    Power { source: String, power_uw: PathBuf },
}

#[derive(Clone, Debug)]
pub struct PowerReading {
    pub source: String,
    pub energy_joules: Option<f64>,
    pub average_power_watts: Option<f64>,
}

pub enum PowerRun {
    Energy {
        source: String,
        energy_uj: PathBuf,
        start_joules: Option<f64>,
    },
    Power {
        source: String,
        stop: Arc<AtomicBool>,
        handle: thread::JoinHandle<Vec<(Duration, f64)>>,
    },
}

impl SysfsPower {
    pub fn autodetect() -> Option<Self> {
        Self::autodetect_from_roots(
            Path::new("/sys/class/powercap"),
            Path::new("/sys/class/hwmon"),
        )
        .or_else(|| {
            Self::autodetect_from_roots(
                Path::new("/sys/devices/virtual/powercap"),
                Path::new("/sys/class/hwmon"),
            )
        })
    }

    pub fn source(&self) -> &str {
        match self {
            SysfsPower::Energy { source, .. } | SysfsPower::Power { source, .. } => source,
        }
    }

    fn autodetect_from_roots(powercap_root: &Path, hwmon_root: &Path) -> Option<Self> {
        find_powercap_energy(powercap_root)
            .filter(|path| read_energy_joules(path).is_some())
            .map(|path| SysfsPower::Energy {
                source: format!("sysfs-powercap:{}", path.display()),
                energy_uj: path,
            })
            .or_else(|| {
                find_hwmon_energy(hwmon_root)
                    .filter(|path| read_energy_joules(path).is_some())
                    .map(|path| SysfsPower::Energy {
                        source: format!("sysfs-hwmon-energy:{}", path.display()),
                        energy_uj: path,
                    })
            })
            .or_else(|| {
                find_hwmon_power(hwmon_root)
                    .filter(|path| read_power_watts(path).is_some())
                    .map(|path| SysfsPower::Power {
                        source: format!("sysfs-hwmon-power:{}", path.display()),
                        power_uw: path,
                    })
            })
    }

    pub fn start(&self) -> PowerRun {
        match self {
            SysfsPower::Energy { source, energy_uj } => PowerRun::Energy {
                source: source.clone(),
                energy_uj: energy_uj.clone(),
                start_joules: read_energy_joules(energy_uj),
            },
            SysfsPower::Power { source, power_uw } => {
                let source = source.clone();
                let power_uw = power_uw.clone();
                let stop = Arc::new(AtomicBool::new(false));
                let sampler_stop = Arc::clone(&stop);
                let handle = thread::spawn(move || {
                    let start = Instant::now();
                    let mut samples = Vec::new();
                    while !sampler_stop.load(Ordering::Relaxed) {
                        if let Some(watts) = read_power_watts(&power_uw) {
                            samples.push((start.elapsed(), watts));
                        }
                        thread::sleep(SAMPLE_INTERVAL);
                    }
                    samples
                });
                PowerRun::Power {
                    source,
                    stop,
                    handle,
                }
            }
        }
    }
}

impl PowerRun {
    pub fn finish(self, elapsed: Duration) -> Option<PowerReading> {
        match self {
            PowerRun::Energy {
                source,
                energy_uj,
                start_joules,
            } => {
                let start_joules = start_joules?;
                let end_joules = read_energy_joules(&energy_uj)?;
                let elapsed_s = elapsed.as_secs_f64();
                let energy_joules = end_joules - start_joules;
                if energy_joules < 0.0 || !energy_joules.is_finite() {
                    return None;
                }
                Some(PowerReading {
                    source,
                    energy_joules: Some(energy_joules),
                    average_power_watts: (elapsed_s > 0.0).then_some(energy_joules / elapsed_s),
                })
            }
            PowerRun::Power {
                source,
                stop,
                handle,
            } => {
                stop.store(true, Ordering::Relaxed);
                let mut samples = handle.join().ok()?;
                let elapsed_s = elapsed.as_secs_f64();
                if elapsed_s <= 0.0 {
                    return None;
                }
                samples.retain(|(time, watts)| {
                    let seconds = time.as_secs_f64();
                    seconds > 0.0 && seconds < elapsed_s && watts.is_finite()
                });
                if samples.is_empty() {
                    return None;
                }
                let average_power_watts =
                    samples.iter().map(|(_, watts)| *watts).sum::<f64>() / samples.len() as f64;
                Some(PowerReading {
                    source,
                    energy_joules: Some(average_power_watts * elapsed_s),
                    average_power_watts: Some(average_power_watts),
                })
            }
        }
    }
}

fn read_energy_joules(path: &Path) -> Option<f64> {
    let microjoules: f64 = fs::read_to_string(path).ok()?.trim().parse().ok()?;
    Some(microjoules / 1_000_000.0)
}

fn read_power_watts(path: &Path) -> Option<f64> {
    let microwatts: f64 = fs::read_to_string(path).ok()?.trim().parse().ok()?;
    Some(microwatts / 1_000_000.0)
}

fn find_powercap_energy(root: &Path) -> Option<PathBuf> {
    let mut candidates = find_all(root, "energy_uj");
    candidates.sort_by_key(|path| powercap_score(path));
    candidates.into_iter().next()
}

fn find_hwmon_energy(root: &Path) -> Option<PathBuf> {
    find_hwmon_sensor(root, SensorKind::Energy)
}

fn find_hwmon_power(root: &Path) -> Option<PathBuf> {
    find_hwmon_sensor(root, SensorKind::Power)
}

fn find_all(root: &Path, filename: &str) -> Vec<PathBuf> {
    let mut matches = Vec::new();
    let mut stack = vec![root.to_path_buf()];
    let mut seen = 0usize;

    while let Some(dir) = stack.pop() {
        if seen > 4096 {
            break;
        }
        seen += 1;

        if let Ok(entries) = fs::read_dir(&dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_dir() {
                    stack.push(path);
                } else if path
                    .file_name()
                    .map(|name| name == filename)
                    .unwrap_or(false)
                {
                    matches.push(path);
                }
            }
        }
    }

    matches
}

fn powercap_score(path: &Path) -> (u8, String) {
    let parent = path.parent().unwrap_or(path);
    let name = read_lower(parent.join("name"));
    let dir = parent
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    let haystack = format!("{name} {dir}");

    let rank = if haystack.contains("package") || dir.ends_with(":0") {
        0
    } else if haystack.contains("core") {
        1
    } else if haystack.contains("psys") {
        2
    } else if haystack.contains("dram") {
        3
    } else {
        4
    };

    (rank, path.to_string_lossy().into_owned())
}

#[derive(Clone, Copy)]
enum SensorKind {
    Energy,
    Power,
}

struct SensorCandidate {
    score: (u8, u8, String),
    path: PathBuf,
}

fn find_hwmon_sensor(root: &Path, kind: SensorKind) -> Option<PathBuf> {
    let entries = fs::read_dir(root).ok()?;
    let mut candidates = Vec::new();

    for entry in entries.flatten() {
        let base = entry.path();
        let chip_name = read_lower(base.join("name"));
        for path in hwmon_input_files(&base, kind) {
            let Some(file_name) = path.file_name().and_then(|name| name.to_str()) else {
                continue;
            };
            let label = hwmon_label(&base, file_name);
            let Some(sensor_rank) = hwmon_sensor_rank(&chip_name, &label) else {
                continue;
            };
            candidates.push(SensorCandidate {
                score: (
                    sensor_rank,
                    hwmon_input_rank(file_name),
                    path.to_string_lossy().into_owned(),
                ),
                path,
            });
        }
    }

    candidates.sort_by(|left, right| left.score.cmp(&right.score));
    candidates
        .into_iter()
        .map(|candidate| candidate.path)
        .next()
}

fn hwmon_input_files(base: &Path, kind: SensorKind) -> Vec<PathBuf> {
    let mut files = Vec::new();
    let Ok(entries) = fs::read_dir(base) else {
        return files;
    };

    for entry in entries.flatten() {
        let path = entry.path();
        let Some(file_name) = path.file_name().and_then(|name| name.to_str()) else {
            continue;
        };

        match kind {
            SensorKind::Energy if is_numbered_input(file_name, "energy") => files.push(path),
            SensorKind::Power if is_numbered_power(file_name) => files.push(path),
            _ => {}
        }
    }

    files
}

fn is_numbered_input(file_name: &str, prefix: &str) -> bool {
    numbered_channel(file_name, prefix)
        .map(|(_, suffix)| suffix == "input")
        .unwrap_or(false)
}

fn is_numbered_power(file_name: &str) -> bool {
    numbered_channel(file_name, "power")
        .map(|(_, suffix)| suffix == "average" || suffix == "input")
        .unwrap_or(false)
}

fn numbered_channel<'a>(file_name: &'a str, prefix: &str) -> Option<(&'a str, &'a str)> {
    let rest = file_name.strip_prefix(prefix)?;
    let split_at = rest.find('_')?;
    let (channel, suffix) = rest.split_at(split_at);
    if channel.is_empty() || !channel.chars().all(|ch| ch.is_ascii_digit()) {
        return None;
    }
    Some((channel, &suffix[1..]))
}

fn hwmon_label(base: &Path, input_name: &str) -> String {
    let Some((channel, _)) =
        numbered_channel(input_name, "power").or_else(|| numbered_channel(input_name, "energy"))
    else {
        return String::new();
    };

    for prefix in ["power", "energy"] {
        let label = read_lower(base.join(format!("{prefix}{channel}_label")));
        if !label.is_empty() {
            return label;
        }
    }

    String::new()
}

fn hwmon_sensor_rank(chip_name: &str, label: &str) -> Option<u8> {
    let haystack = format!("{chip_name} {label}");
    if haystack.contains("gpu") && !haystack.contains("cpu") {
        return None;
    }
    if contains_any(
        &haystack,
        &[
            "cpu",
            "package",
            "core",
            "rapl",
            "ppt",
            "svi2_p_core",
            "vdd_cpu",
        ],
    ) {
        return Some(0);
    }
    if contains_any(
        &haystack,
        &["soc", "vcore", "vdd", "processor", "scmi", "arm"],
    ) {
        return Some(1);
    }
    if contains_any(
        chip_name,
        &[
            "amd_energy",
            "zenpower",
            "k10temp",
            "ina3221",
            "ina2",
            "pmbus",
            "power_meter",
            "ltc",
            "adm127",
        ],
    ) {
        return Some(2);
    }
    None
}

fn hwmon_input_rank(file_name: &str) -> u8 {
    if file_name.contains("_average") {
        0
    } else {
        1
    }
}

fn contains_any(haystack: &str, needles: &[&str]) -> bool {
    needles.iter().any(|needle| haystack.contains(needle))
}

fn read_lower(path: impl AsRef<Path>) -> String {
    fs::read_to_string(path)
        .map(|content| content.trim().to_ascii_lowercase())
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_root(name: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "fpx_sysfs_power_{name}_{}_{}",
            std::process::id(),
            nanos
        ));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        root
    }

    #[test]
    fn detects_x86_package_powercap_energy() {
        let root = temp_root("powercap");
        let powercap = root.join("powercap");
        let package = powercap.join("intel-rapl:0");
        let core = powercap.join("intel-rapl:0:0");
        fs::create_dir_all(&package).unwrap();
        fs::create_dir_all(&core).unwrap();
        fs::write(package.join("name"), "package-0\n").unwrap();
        fs::write(package.join("energy_uj"), "1000000\n").unwrap();
        fs::write(core.join("name"), "core\n").unwrap();
        fs::write(core.join("energy_uj"), "2000000\n").unwrap();

        match SysfsPower::autodetect_from_roots(&powercap, &root.join("hwmon")) {
            Some(SysfsPower::Energy { energy_uj, .. }) => {
                assert_eq!(energy_uj, package.join("energy_uj"));
            }
            other => panic!("unexpected power source: {other:?}"),
        }
    }

    #[test]
    fn detects_amd_hwmon_energy_without_powercap() {
        let root = temp_root("amd_hwmon");
        let powercap = root.join("powercap");
        let hwmon = root.join("hwmon");
        let chip = hwmon.join("hwmon0");
        fs::create_dir_all(&powercap).unwrap();
        fs::create_dir_all(&chip).unwrap();
        fs::write(chip.join("name"), "amd_energy\n").unwrap();
        fs::write(chip.join("energy1_label"), "package-0\n").unwrap();
        fs::write(chip.join("energy1_input"), "42000000\n").unwrap();

        match SysfsPower::autodetect_from_roots(&powercap, &hwmon) {
            Some(SysfsPower::Energy { energy_uj, .. }) => {
                assert_eq!(energy_uj, chip.join("energy1_input"));
            }
            other => panic!("unexpected power source: {other:?}"),
        }
    }

    #[test]
    fn detects_aarch64_cpu_power_hwmon() {
        let root = temp_root("aarch64_hwmon");
        let powercap = root.join("powercap");
        let hwmon = root.join("hwmon");
        let chip = hwmon.join("hwmon0");
        fs::create_dir_all(&powercap).unwrap();
        fs::create_dir_all(&chip).unwrap();
        fs::write(chip.join("name"), "ina3221x\n").unwrap();
        fs::write(chip.join("power1_label"), "VDD_CPU\n").unwrap();
        fs::write(chip.join("power1_input"), "1500000\n").unwrap();

        match SysfsPower::autodetect_from_roots(&powercap, &hwmon) {
            Some(SysfsPower::Power { power_uw, .. }) => {
                assert_eq!(power_uw, chip.join("power1_input"));
            }
            other => panic!("unexpected power source: {other:?}"),
        }
    }

    #[test]
    fn ignores_gpu_only_hwmon_power() {
        let root = temp_root("gpu_hwmon");
        let powercap = root.join("powercap");
        let hwmon = root.join("hwmon");
        let chip = hwmon.join("hwmon0");
        fs::create_dir_all(&powercap).unwrap();
        fs::create_dir_all(&chip).unwrap();
        fs::write(chip.join("name"), "amdgpu\n").unwrap();
        fs::write(chip.join("power1_label"), "GPU\n").unwrap();
        fs::write(chip.join("power1_average"), "99000000\n").unwrap();

        assert!(SysfsPower::autodetect_from_roots(&powercap, &hwmon).is_none());
    }
}
