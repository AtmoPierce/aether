#![cfg(feature = "std")]

// src/performance/tests/power_and_timing.rs
use crate::performance::clock::now;
use crate::performance::performance::*;
use crate::performance::power::*;
use crate::performance::timer::*;
use std::{fs, thread, time::Duration};

// periodic sampler -> writes to PowerCsv
struct Sampler {
    stop: std::sync::Arc<std::sync::atomic::AtomicBool>,
    handle: Option<std::thread::JoinHandle<()>>,
}
impl Sampler {
    fn start(cpu: CPUEnergy, sink: &'static PowerCsv, period_ms: u64) -> Self {
        use std::sync::{
            atomic::{AtomicBool, Ordering},
            Arc,
        };
        let stop = Arc::new(AtomicBool::new(false));
        let stop_c = stop.clone();
        let handle = std::thread::spawn(move || {
            eprintln!("Power source: {}", cpu.debug_source_path());
            let mut last_j = cpu.read_total_j();
            let mut last_t = now();
            while !stop_c.load(Ordering::Relaxed) {
                thread::sleep(Duration::from_millis(period_ms));
                let t1 = now();
                match &cpu {
                    CPUEnergy::Powercap { .. } => {
                        let j1 = cpu.read_total_j();
                        if let (Some(j1), Some(j0)) = (j1, last_j) {
                            let dt = (t1 - last_t).as_secs_f64();
                            let p = if dt > 0.0 { Some((j1 - j0) / dt) } else { None };
                            sink.sample("cpu_pkg", p, Some(j1));
                        }
                        last_j = cpu.read_total_j();
                        last_t = t1;
                    }
                    CPUEnergy::Hwmon { .. } => {
                        let w = cpu.read_watts();
                        sink.sample("cpu_pkg", w, None);
                        last_t = t1;
                    }
                    CPUEnergy::None => { /* nothing; keep file but rows may be empty */ }
                }
            }
        });
        Self {
            stop,
            handle: Some(handle),
        }
    }
    fn stop(mut self) {
        self.stop.store(true, std::sync::atomic::Ordering::Relaxed);
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

// tiny workload to time
fn busy_compute(iters: usize) -> f64 {
    let mut acc = 0.0f64;
    let mut x = 1.000_001_f64;
    for _ in 0..iters {
        x = (x.sin().cos().tan()).abs() + 1e-9;
        acc += x;
    }
    acc
}

#[test]
fn timing_and_power_csv() {
    let tmp = std::env::temp_dir();
    let timing_path = tmp.join(format!("astrum_timing_{}.csv", std::process::id()));
    let power_path = tmp.join(format!("astrum_power_{}.csv", std::process::id()));

    let timing: &'static TimingCsv =
        Box::leak(Box::new(TimingCsv::new(timing_path.to_str().unwrap())));
    let power: &'static PowerCsv = Box::leak(Box::new(PowerCsv::new(power_path.to_str().unwrap())));

    let cpu = CPUEnergy::autodetect();
    let sampler = Sampler::start(cpu.clone(), power, 100);

    {
        let _t = TimeSpan::new("stage_rollout", timing);
        let _ = busy_compute(5_000_000);
    }
    {
        let _t = TimeSpan::new("stage_gradients", timing);
        let _ = busy_compute(6_000_000);
    }
    {
        let _t = TimeSpan::new("stage_line_search", timing);
        let _ = busy_compute(4_000_000);
    }

    thread::sleep(Duration::from_millis(300));
    sampler.stop();

    // Timing CSV must have at least one row with our header format
    let tcsv = fs::read_to_string(&timing_path).expect("timing csv");
    assert!(tcsv.lines().next().unwrap().trim() == "unix_ns,stage,elapsed_ms");
    assert!(tcsv.contains(",stage_rollout,"));

    // Power CSV header always present; rows appear if sensors available
    let pcsv = fs::read_to_string(&power_path).expect("power csv");
    assert!(pcsv.lines().next().unwrap().trim() == "unix_ns,source,power_W,energy_J_total");
    eprintln!("timing -> {}", timing_path.display());
    eprintln!("power  -> {}", power_path.display());
}
