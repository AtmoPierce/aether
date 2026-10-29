// aether-fpx-bench/src/stats.rs
use crate::types::{Row, Summary};

pub fn summarize(iters: usize, rows: &[&Row]) -> Summary {
    let mut s = Summary {
        sum_ins: 0.0,
        sum_cyc: 0.0,
        sum_br: 0.0,
        sum_bm: 0.0,
        sum_time_s: 0.0,
        n: 0,
        cpi_agg: None,
        ipc_agg: None,
        ns_per_call: None,
        ins_per_call: None,
        cyc_per_call: None,
        bm_rate: None,
    };

    for r in rows {
        if let Some(v) = r.instructions {
            s.sum_ins += v;
        }
        if let Some(v) = r.cycles {
            s.sum_cyc += v;
        }
        if let Some(v) = r.branches {
            s.sum_br += v;
        }
        if let Some(v) = r.branch_misses {
            s.sum_bm += v;
        }
        if let Some(v) = r.time_elapsed_s {
            s.sum_time_s += v;
        }
        s.n += 1;
    }

    if s.sum_ins > 0.0 {
        s.cpi_agg = Some(s.sum_cyc / s.sum_ins);
        s.ipc_agg = Some(s.sum_ins / s.sum_cyc);
    }

    let calls = (iters as f64) * (s.n as f64);
    if calls > 0.0 {
        s.ns_per_call = Some(1e9 * s.sum_time_s / calls);
        s.ins_per_call = Some(s.sum_ins / calls);
        s.cyc_per_call = Some(s.sum_cyc / calls);
    }

    if s.sum_br > 0.0 {
        s.bm_rate = Some(s.sum_bm / s.sum_br);
    }

    s
}

/// Compute mean (+/-) sample standard deviation
pub fn stats(xs: Vec<f64>) -> Option<(f64, f64, usize)> {
    let n = xs.len();
    if n == 0 {
        return None;
    }

    let mean = xs.iter().sum::<f64>() / n as f64;

    if n == 1 {
        return Some((mean, 0.0, 1));
    }

    let var = xs
        .iter()
        .map(|x| {
            let d = x - mean;
            d * d
        })
        .sum::<f64>()
        / (n as f64 - 1.0);

    Some((mean, var.sqrt(), n))
}
