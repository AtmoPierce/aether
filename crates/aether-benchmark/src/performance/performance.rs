#![cfg(feature = "std")]
// src/performance/performance.rs
use crate::performance::clock::unix_nanos_now;
use std::{io::Write, sync::Mutex};

/// Timing CSV: unix_ns, stage, elapsed_ms
pub struct TimingCsv(pub Mutex<std::fs::File>);
impl TimingCsv {
    pub fn new(path: &str) -> Self {
        let mut f = std::fs::File::create(path).expect("timing csv");
        let _ = writeln!(f, "unix_ns,stage,elapsed_ms");
        Self(Mutex::new(f))
    }
    #[inline]
    pub fn timing_ms(&self, stage: &str, ms: f64) {
        let _ = writeln!(
            self.0.lock().unwrap(),
            "{},{},{}",
            unix_nanos_now(),
            stage,
            ms
        );
    }
}

/// Power CSV: unix_ns, source, power_W, energy_J_total
pub struct PowerCsv(pub Mutex<std::fs::File>);
impl PowerCsv {
    pub fn new(path: &str) -> Self {
        let mut f = std::fs::File::create(path).expect("power csv");
        let _ = writeln!(f, "unix_ns,source,power_W,energy_J_total");
        Self(Mutex::new(f))
    }
    /// Write one sample; pass `None` for columns you can't provide
    #[inline]
    pub fn sample(&self, source: &str, power_w: Option<f64>, energy_j_total: Option<f64>) {
        let w = power_w.map(|v| v.to_string()).unwrap_or_else(String::new);
        let j = energy_j_total
            .map(|v| v.to_string())
            .unwrap_or_else(String::new);
        let _ = writeln!(
            self.0.lock().unwrap(),
            "{},{},{},{}",
            unix_nanos_now(),
            source,
            w,
            j
        );
    }
}
