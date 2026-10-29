// aether-fpx-bench/src/report.rs

use crate::stats::stats;
use crate::types::Row;

pub fn print_summary(op: &str, mode: &str, prec: &str, iters: usize, rows: &[&Row]) {
    let ins = stats(rows.iter().filter_map(|r| r.instructions).collect());
    let cyc = stats(rows.iter().filter_map(|r| r.cycles).collect());
    let br = stats(rows.iter().filter_map(|r| r.branches).collect());
    let bm = stats(rows.iter().filter_map(|r| r.branch_misses).collect());
    let te = stats(rows.iter().filter_map(|r| r.time_elapsed_s).collect());

    let sum_ins: f64 = rows.iter().filter_map(|r| r.instructions).sum();
    let sum_cyc: f64 = rows.iter().filter_map(|r| r.cycles).sum();

    let cpi_agg = if sum_ins > 0.0 {
        Some(sum_cyc / sum_ins)
    } else {
        None
    };
    let ipc_agg = if sum_cyc > 0.0 {
        Some(sum_ins / sum_cyc)
    } else {
        None
    };

    let ins_per_call = ins.map(|(mean, _, _)| mean / iters as f64);
    let cyc_per_call = cyc.map(|(mean, _, _)| mean / iters as f64);
    let ns_per_call = te.map(|(mean, _, _)| 1e9 * mean / iters as f64);

    println!("op={op} mode={mode} prec={prec}");
    print_stat("  instructions", ins);
    print_stat("  cycles      ", cyc);
    print_stat("  branches    ", br);
    print_stat("  br-misses   ", bm);
    print_stat("  elapsed(s)  ", te);

    println!(
        "  CPI (agg)   : {}   IPC (agg): {}",
        fmt(cpi_agg),
        fmt(ipc_agg),
    );

    println!(
        "  per-call    : ins/call={}  cyc/call={}  ns/call={}",
        fmt(ins_per_call),
        fmt(cyc_per_call),
        fmt(ns_per_call),
    );
}

pub fn print_comparison(
    op: &str,
    prec: &str,
    iters: usize,
    lhs_label: &str,
    lhs: &[&Row],
    rhs_label: &str,
    rhs: &[&Row],
) {
    let l = crate::stats::summarize(iters, lhs);
    let r = crate::stats::summarize(iters, rhs);

    let (speedup, speedup_label) = if lhs_label == "nat" {
        (
            match (l.ns_per_call, r.ns_per_call) {
                (Some(nat), Some(other)) if other > 0.0 => Some(nat / other),
                _ => None,
            },
            format!("speedup({lhs_label}/{rhs_label})"),
        )
    } else if rhs_label == "nat" {
        (
            match (l.ns_per_call, r.ns_per_call) {
                (Some(other), Some(nat)) if other > 0.0 => Some(nat / other),
                _ => None,
            },
            format!("speedup({rhs_label}/{lhs_label})"),
        )
    } else {
        (
            match (l.ns_per_call, r.ns_per_call) {
                (Some(a), Some(b)) if a > 0.0 => Some(b / a),
                _ => None,
            },
            format!("speedup({rhs_label}/{lhs_label})"),
        )
    };

    println!("op={op} prec={prec}");
    println!(
        "  ns/call     : lhs({})={}  rhs({})={}  {}={}",
        lhs_label,
        fmt(l.ns_per_call),
        rhs_label,
        fmt(r.ns_per_call),
        speedup_label,
        fmt(speedup),
    );
}

fn print_stat(label: &str, st: Option<(f64, f64, usize)>) {
    if let Some((mean, sd, n)) = st {
        println!("{label}: {mean:.6} (+/-) {sd:.6}  (n={n})");
    } else {
        println!("{label}: (no data)");
    }
}

fn fmt(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{x:.6}"),
        None => "NA".into(),
    }
}

#[macro_export]
macro_rules! report_perf {
    (
        title: $title:expr,

        axes {
            op: [$($op:expr),+ $(,)?],
            mode: [$($mode:expr),+ $(,)?],
            precision: [$($prec:expr),+ $(,)?],
        }

        iters: $iters:expr,
        rows: $rows:expr,

        summary { $($summary:ident),* $(,)? }

        $(, compare { mode: $lhs:expr => $rhs:expr })?
    ) => {{
        println!("\n=== {} ===", $title);

        // bind axes once
        let ops = [$($op),+];
        let modes = [$($mode),+];
        let precs = [$($prec),+];

        for op in ops.iter() {
            for mode in modes.iter() {
                for prec in precs.iter() {
                    let rows: Vec<_> = $rows.iter()
                        .filter(|r|
                            r.op == *op &&
                            r.mode == *mode &&
                            r.precision == *prec
                        )
                        .collect();

                    if !rows.is_empty() {
                        $crate::report::print_summary(op, mode, prec, $iters, &rows);
                    }
                }
            }
        }

        $(
            println!("\nDirect comparison ({} vs {}):", $lhs, $rhs);

            for op in ops.iter() {
                for prec in precs.iter() {
                    let lhs_rows: Vec<_> = $rows.iter()
                        .filter(|r|
                            r.op == *op &&
                            r.mode == $lhs &&
                            r.precision == *prec
                        )
                        .collect();

                    let rhs_rows: Vec<_> = $rows.iter()
                        .filter(|r|
                            r.op == *op &&
                            r.mode == $rhs &&
                            r.precision == *prec
                        )
                        .collect();

                    if !lhs_rows.is_empty() && !rhs_rows.is_empty() {
                        print_comparison(op, prec, $iters, $lhs, &lhs_rows, $rhs, &rhs_rows);
                    }
                }
            }
        )?
    }};
}
