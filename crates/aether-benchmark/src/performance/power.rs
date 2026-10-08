// src/performance/power.rs
#![cfg(feature = "std")]
use std::{
    fs,
    path::{Path, PathBuf},
};

#[derive(Clone, Debug)]
pub enum CPUEnergy {
    Powercap { energy_uj: PathBuf },    // cumulative energy
    HwmonEnergy { energy_uj: PathBuf }, // cumulative energy
    Hwmon { power_uw: PathBuf },        // instantaneous µW
    None,
}

impl CPUEnergy {
    pub fn autodetect() -> Self {
        Self::autodetect_from_roots(
            Path::new("/sys/class/powercap"),
            Path::new("/sys/class/hwmon"),
        )
    }

    fn autodetect_from_roots(powercap_root: &Path, hwmon_root: &Path) -> Self {
        if let Some(p) = find_powercap_energy(powercap_root) {
            return CPUEnergy::Powercap { energy_uj: p };
        }
        if let Some(p) = find_hwmon_cpu_energy(hwmon_root) {
            return CPUEnergy::HwmonEnergy { energy_uj: p };
        }
        if let Some(p) = find_hwmon_cpu_power(hwmon_root) {
            return CPUEnergy::Hwmon { power_uw: p };
        }
        CPUEnergy::None
    }

    pub fn read_total_j(&self) -> Option<f64> {
        match self {
            CPUEnergy::Powercap { energy_uj } | CPUEnergy::HwmonEnergy { energy_uj } => {
                let s = fs::read_to_string(energy_uj).ok()?;
                let uj: u64 = s.trim().parse().ok()?;
                Some(uj as f64 / 1_000_000.0)
            }
            _ => None,
        }
    }
    pub fn read_watts(&self) -> Option<f64> {
        match self {
            CPUEnergy::Hwmon { power_uw } => {
                let s = fs::read_to_string(power_uw).ok()?;
                let uw: f64 = s.trim().parse().ok()?;
                Some(uw / 1_000_000.0)
            }
            _ => None,
        }
    }
    pub fn debug_source_path(&self) -> String {
        match self {
            CPUEnergy::Powercap { energy_uj } => format!("powercap:{}", energy_uj.display()),
            CPUEnergy::HwmonEnergy { energy_uj } => format!("hwmon-energy:{}", energy_uj.display()),
            CPUEnergy::Hwmon { power_uw } => format!("hwmon:{}", power_uw.display()),
            CPUEnergy::None => "none".to_string(),
        }
    }
}

fn find_all(root: &Path, filename: &str) -> Vec<PathBuf> {
    let mut matches = Vec::new();
    let mut stack = vec![root.to_path_buf()];
    let max_nodes = 4096;
    let mut seen = 0usize;
    while let Some(dir) = stack.pop() {
        if seen > max_nodes {
            break;
        }
        seen += 1;
        if let Ok(read) = fs::read_dir(&dir) {
            for e in read.flatten() {
                let p = e.path();
                if p.is_dir() {
                    stack.push(p);
                } else if p.file_name().map(|n| n == filename).unwrap_or(false) {
                    matches.push(p);
                }
            }
        }
    }
    matches
}

fn find_powercap_energy(root: &Path) -> Option<PathBuf> {
    let mut candidates = find_all(root, "energy_uj");
    candidates.sort_by_key(|path| powercap_score(path));
    candidates.into_iter().next()
}

fn powercap_score(path: &Path) -> (u8, String) {
    let parent = path.parent().unwrap_or(path);
    let name = read_trimmed_lower(parent.join("name"));
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

fn find_hwmon_cpu_energy(hwmon_root: &Path) -> Option<PathBuf> {
    find_hwmon_sensor(hwmon_root, SensorKind::Energy)
}

fn find_hwmon_cpu_power(hwmon_root: &Path) -> Option<PathBuf> {
    find_hwmon_sensor(hwmon_root, SensorKind::Power)
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

fn find_hwmon_sensor(hwmon_root: &Path, kind: SensorKind) -> Option<PathBuf> {
    let dir = fs::read_dir(hwmon_root).ok()?;
    let mut candidates = Vec::new();

    for ent in dir.flatten() {
        let base = ent.path();
        let chip_name = read_trimmed_lower(base.join("name"));
        for path in hwmon_input_files(&base, kind) {
            let Some(file_name) = path.file_name().and_then(|name| name.to_str()) else {
                continue;
            };
            let label = hwmon_label(&base, file_name);
            let Some(sensor_rank) = hwmon_sensor_rank(&chip_name, &label) else {
                continue;
            };
            let input_rank = hwmon_input_rank(file_name);
            candidates.push(SensorCandidate {
                score: (sensor_rank, input_rank, path.to_string_lossy().into_owned()),
                path,
            });
        }
    }

    candidates.sort_by(|a, b| a.score.cmp(&b.score));
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
            SensorKind::Energy if is_numbered_hwmon_input(file_name, "energy") => files.push(path),
            SensorKind::Power if is_numbered_hwmon_power(file_name) => files.push(path),
            _ => {}
        }
    }

    files
}

fn is_numbered_hwmon_input(file_name: &str, prefix: &str) -> bool {
    numbered_hwmon_channel(file_name, prefix)
        .map(|(_, suffix)| suffix == "input")
        .unwrap_or(false)
}

fn is_numbered_hwmon_power(file_name: &str) -> bool {
    numbered_hwmon_channel(file_name, "power")
        .map(|(_, suffix)| suffix == "average" || suffix == "input")
        .unwrap_or(false)
}

fn numbered_hwmon_channel<'a>(file_name: &'a str, prefix: &str) -> Option<(&'a str, &'a str)> {
    let rest = file_name.strip_prefix(prefix)?;
    let split_at = rest.find('_')?;
    let (channel, suffix) = rest.split_at(split_at);
    if channel.is_empty() || !channel.chars().all(|ch| ch.is_ascii_digit()) {
        return None;
    }
    Some((channel, &suffix[1..]))
}

fn hwmon_label(base: &Path, input_name: &str) -> String {
    let Some((channel, _)) = numbered_hwmon_channel(input_name, "power")
        .or_else(|| numbered_hwmon_channel(input_name, "energy"))
    else {
        return String::new();
    };

    for prefix in ["power", "energy"] {
        let label = read_trimmed_lower(base.join(format!("{prefix}{channel}_label")));
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

fn read_trimmed_lower(path: impl AsRef<Path>) -> String {
    fs::read_to_string(path)
        .map(|content| content.trim().to_ascii_lowercase())
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(name: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "aether_benchmark_power_{name}_{}_{}",
            std::process::id(),
            nanos
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        path
    }

    #[test]
    fn autodetect_prefers_package_powercap_energy_on_x86_layouts() {
        let root = temp_dir("powercap");
        let powercap = root.join("powercap");
        let hwmon = root.join("hwmon");
        let package = powercap.join("intel-rapl:0");
        let core = powercap.join("intel-rapl:0:0");
        fs::create_dir_all(&package).unwrap();
        fs::create_dir_all(&core).unwrap();
        fs::write(package.join("name"), "package-0\n").unwrap();
        fs::write(package.join("energy_uj"), "125000000\n").unwrap();
        fs::write(core.join("name"), "core\n").unwrap();
        fs::write(core.join("energy_uj"), "5000000\n").unwrap();

        match CPUEnergy::autodetect_from_roots(&powercap, &hwmon) {
            CPUEnergy::Powercap { energy_uj } => {
                assert_eq!(energy_uj, package.join("energy_uj"));
                assert_eq!(
                    CPUEnergy::Powercap { energy_uj }.read_total_j(),
                    Some(125.0)
                );
            }
            other => panic!("unexpected power source: {other:?}"),
        }
    }

    #[test]
    fn autodetect_uses_amd_hwmon_energy_when_powercap_is_absent() {
        let root = temp_dir("amd_hwmon_energy");
        let powercap = root.join("powercap");
        let chip = root.join("hwmon").join("hwmon0");
        fs::create_dir_all(&powercap).unwrap();
        fs::create_dir_all(&chip).unwrap();
        fs::write(chip.join("name"), "amd_energy\n").unwrap();
        fs::write(chip.join("energy1_label"), "package-0\n").unwrap();
        fs::write(chip.join("energy1_input"), "42000000\n").unwrap();

        match CPUEnergy::autodetect_from_roots(&powercap, &root.join("hwmon")) {
            CPUEnergy::HwmonEnergy { energy_uj } => {
                assert_eq!(energy_uj, chip.join("energy1_input"));
                assert_eq!(
                    CPUEnergy::HwmonEnergy { energy_uj }.read_total_j(),
                    Some(42.0)
                );
            }
            other => panic!("unexpected power source: {other:?}"),
        }
    }

    #[test]
    fn autodetect_uses_aarch64_board_power_monitors() {
        let root = temp_dir("aarch64_hwmon_power");
        let powercap = root.join("powercap");
        let chip = root.join("hwmon").join("hwmon0");
        fs::create_dir_all(&powercap).unwrap();
        fs::create_dir_all(&chip).unwrap();
        fs::write(chip.join("name"), "ina3221x\n").unwrap();
        fs::write(chip.join("power1_label"), "VDD_CPU\n").unwrap();
        fs::write(chip.join("power1_input"), "1500000\n").unwrap();

        match CPUEnergy::autodetect_from_roots(&powercap, &root.join("hwmon")) {
            CPUEnergy::Hwmon { power_uw } => {
                assert_eq!(power_uw, chip.join("power1_input"));
                assert_eq!(CPUEnergy::Hwmon { power_uw }.read_watts(), Some(1.5));
            }
            other => panic!("unexpected power source: {other:?}"),
        }
    }

    #[test]
    fn autodetect_does_not_report_gpu_only_hwmon_as_cpu_power() {
        let root = temp_dir("gpu_hwmon_power");
        let powercap = root.join("powercap");
        let chip = root.join("hwmon").join("hwmon0");
        fs::create_dir_all(&powercap).unwrap();
        fs::create_dir_all(&chip).unwrap();
        fs::write(chip.join("name"), "amdgpu\n").unwrap();
        fs::write(chip.join("power1_label"), "GPU\n").unwrap();
        fs::write(chip.join("power1_average"), "99000000\n").unwrap();

        assert!(matches!(
            CPUEnergy::autodetect_from_roots(&powercap, &root.join("hwmon")),
            CPUEnergy::None
        ));
    }
}
