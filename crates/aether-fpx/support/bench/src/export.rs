// aether-fpx-bench/src/export.rs
use std::fs;
use std::path::Path;

use crate::format::opt_f64;
use crate::types::Row;

/// Write benchmark rows to a CSV file.
///
/// Schema:
/// op,mode,precision,run,instructions,cycles,branches,branch_misses,time_elapsed_s,
/// energy_joules,average_power_watts,power_source,iterations
pub fn write_csv<P: AsRef<Path>>(path: P, rows: &[Row]) -> Result<(), String> {
    let mut out = String::new();

    out.push_str(
        "op,mode,precision,run,\
         instructions,cycles,branches,branch_misses,time_elapsed_s,\
         energy_joules,average_power_watts,power_source,iterations\n",
    );

    for r in rows {
        out.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
            r.op,
            r.mode,
            r.precision,
            r.run,
            opt_f64(r.instructions),
            opt_f64(r.cycles),
            opt_f64(r.branches),
            opt_f64(r.branch_misses),
            opt_f64(r.time_elapsed_s),
            opt_f64(r.energy_joules),
            opt_f64(r.average_power_watts),
            r.power_source.as_deref().unwrap_or(""),
            r.iterations,
        ));
    }

    fs::write(path, out).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_csv(name: &str) -> std::path::PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "fpx_export_{name}_{}_{}.csv",
            std::process::id(),
            nanos
        ))
    }

    #[test]
    fn writes_power_columns_when_rows_have_power_readings() {
        let path = temp_csv("power");
        let rows = [Row {
            op: "tan".into(),
            mode: "fpx".into(),
            precision: "f16".into(),
            run: 0,
            iterations: 1_000_000,
            instructions: Some(1.0),
            cycles: Some(2.0),
            branches: Some(3.0),
            branch_misses: Some(4.0),
            time_elapsed_s: Some(0.5),
            energy_joules: Some(1.25),
            average_power_watts: Some(2.5),
            power_source: Some("perf:power/energy-pkg/".into()),
        }];

        write_csv(&path, &rows).unwrap();
        let csv = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(path);

        assert!(csv.contains("energy_joules,average_power_watts,power_source,iterations"));
        assert!(csv.contains(",1.25,2.5,perf:power/energy-pkg/,1000000"));
    }
}
