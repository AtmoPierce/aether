// src/performance/power.rs
#![cfg(feature = "std")]
use std::{
    fs,
    path::{Path, PathBuf},
};

#[derive(Clone, Debug)]
pub enum CPUEnergy {
    Powercap { energy_uj: PathBuf }, // cumulative energy
    Hwmon { power_uw: PathBuf },     // instantaneous µW
    None,
}

impl CPUEnergy {
    pub fn autodetect() -> Self {
        if let Some(p) = find_first(&PathBuf::from("/sys/class/powercap"), "energy_uj") {
            return CPUEnergy::Powercap { energy_uj: p };
        }
        if let Some(p) = find_hwmon_cpu_power() {
            return CPUEnergy::Hwmon { power_uw: p };
        }
        CPUEnergy::None
    }

    pub fn read_total_j(&self) -> Option<f64> {
        match self {
            CPUEnergy::Powercap { energy_uj } => {
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
            CPUEnergy::Hwmon { power_uw } => format!("hwmon:{}", power_uw.display()),
            CPUEnergy::None => "none".to_string(),
        }
    }
}

// --- helpers ---
fn find_first(root: &Path, filename: &str) -> Option<PathBuf> {
    // small DFS without external crates
    let mut stack = vec![root.to_path_buf()];
    let max_nodes = 4096; // safety
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
                    return Some(p);
                }
            }
        }
    }
    None
}

fn find_hwmon_cpu_power() -> Option<PathBuf> {
    let hw = PathBuf::from("/sys/class/hwmon");
    let dir = fs::read_dir(&hw).ok()?;
    for ent in dir.flatten() {
        let base = ent.path();
        let name = fs::read_to_string(base.join("name")).unwrap_or_default();
        // Try CPU sensors first
        let is_cpu = name.contains("k10temp") || name.contains("zenpower");
        if !is_cpu {
            continue;
        }
        for n in 1..=8 {
            for key in [
                "power1_average",
                "power1_input",
                "power2_average",
                "power2_input",
            ] {
                let p = base.join(key.replace('1', &n.to_string()));
                if p.exists() {
                    return Some(p);
                }
            }
        }
    }
    None
}
