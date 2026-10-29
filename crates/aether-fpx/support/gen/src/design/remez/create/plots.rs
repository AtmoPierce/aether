// aether-fpx-gen/src/design/remez/create/plots.rs

use core::f128;
use std::fs;
use std::path::PathBuf;

use super::coefficients::{PolyDomain, TrigKernel};
use super::config::out_dir_for_use;
use crate::remez::RemezApprox;

use aether_core::real::Real;
use aether_viz::plot_error_on_grid;
use aether_fpx_bench::ulp::ulp_dist_f64;

#[inline(always)]
fn tan_runtime_plot_true(x: f64) -> f64 {
    x.tan()
}

#[inline(always)]
fn tan_solve_plot_true(z: f64) -> f64 {
    let x = z.sqrt();

    if x == 0.0 {
        1.0 / 3.0
    } else {
        (x.tan() - x) / (x * x * x)
    }
}

/// Write per-segment SVG plots using aether_viz (plotters).
///
/// Produces:
///  - `{file_prefix}__err_x.svg`      : error in runtime x-domain
///  - `{file_prefix}__err_solve.svg`  : error in solver domain (u)
///  - `{file_prefix}__ulp_x.svg`      : ULP error in runtime x-domain
///  - `{file_prefix}__ulp_solve.svg`  : ULP error in solver domain (u)
pub fn write_segment_plots<const N: usize, const K: usize>(
    _func_name: &str,
    file_prefix: &str,
    dom: PolyDomain,
    kernel: Option<TrigKernel>,
    rz: &RemezApprox<f128, N, K>,
    interval_x: (f128, f128),
    interval_solve: (f128, f128),
    samples: usize,
    f_runtime: fn(f128) -> f128,
    f_solve: fn(f128) -> f128,
) {
    let out_dir: PathBuf = out_dir_for_use("plots", "target/plots");
    fs::create_dir_all(&out_dir).unwrap();

    let path_x = out_dir.join(format!("{file_prefix}__err_x.svg"));
    let path_s = out_dir.join(format!("{file_prefix}__err_solve.svg"));
    let path_ulp_x = out_dir.join(format!("{file_prefix}__ulp_x.svg"));
    let path_ulp_s = out_dir.join(format!("{file_prefix}__ulp_solve.svg"));

    // ---------------------------------------------------------------------------------------------
    // Shared runtime evaluator
    // ---------------------------------------------------------------------------------------------
    let eval_runtime = |x: f128| -> f128 {
        if _func_name == "log" {
            let u = x - f128::ONE;
            let s = u / (x + f128::ONE);
            let z = s * s;
            return 2.0_f128 * s * (f128::ONE + z * rz.eval_poly(z));
        }

        match dom {
            PolyDomain::X => rz.eval_poly(x),
            PolyDomain::Z => {
                let k = kernel.expect("Z-domain plots missing kernel");
                let z = x * x;
                match k {
                    TrigKernel::Cos => rz.eval_poly(z),
                    TrigKernel::Sin => {
                        if x == f128::ZERO {
                            f128::ZERO
                        } else {
                            x * rz.eval_poly(z)
                        }
                    }
                    TrigKernel::Tan => x + x * z * rz.eval_poly(z),
                }
            }
        }
    };

    let runtime_true_f64 = |x: f64| -> f64 {
        if _func_name == "tan" {
            tan_runtime_plot_true(x)
        } else {
            f_runtime(x as f128) as f64
        }
    };

    let solve_true_f64 = |u: f64| -> f64 {
        if _func_name == "tan" {
            tan_solve_plot_true(u)
        } else {
            f_solve(u as f128) as f64
        }
    };

    // ---------------------------------------------------------------------------------------------
    // Runtime-domain error over x
    // ---------------------------------------------------------------------------------------------
    {
        let (a, b) = (interval_x.0 as f64, interval_x.1 as f64);

        let f_approx = |x: f64| -> f64 { eval_runtime(x as f128) as f64 };

        plot_error_on_grid(
            f_approx,
            runtime_true_f64,
            a,
            b,
            samples,
            Some(path_x.to_str().unwrap()),
        )
        .unwrap();
    }

    // ---------------------------------------------------------------------------------------------
    // Runtime-domain ULP error over x
    // ---------------------------------------------------------------------------------------------
    {
        let (a, b) = (interval_x.0 as f64, interval_x.1 as f64);

        let f_ulp = |x: f64| -> f64 {
            let approx = eval_runtime(x as f128) as f64;
            let true_v = runtime_true_f64(x);
            ulp_dist_f64(approx, true_v) as f64
        };

        plot_error_on_grid(
            f_ulp,
            |_| 0.0,
            a,
            b,
            samples,
            Some(path_ulp_x.to_str().unwrap()),
        )
        .unwrap();
    }

    // ---------------------------------------------------------------------------------------------
    // Solver-domain error over u
    // ---------------------------------------------------------------------------------------------
    {
        let (a, b) = (interval_solve.0 as f64, interval_solve.1 as f64);

        let f_approx = |u: f64| -> f64 { rz.eval_poly(u as f128) as f64 };

        plot_error_on_grid(
            f_approx,
            solve_true_f64,
            a,
            b,
            samples,
            Some(path_s.to_str().unwrap()),
        )
        .unwrap();
    }

    // ---------------------------------------------------------------------------------------------
    // Solver-domain ULP error over u
    // ---------------------------------------------------------------------------------------------
    {
        let (a, b) = (interval_solve.0 as f64, interval_solve.1 as f64);

        let f_ulp = |u: f64| -> f64 {
            let approx = rz.eval_poly(u as f128) as f64;
            let true_v = solve_true_f64(u);
            ulp_dist_f64(approx, true_v) as f64
        };

        plot_error_on_grid(
            f_ulp,
            |_| 0.0,
            a,
            b,
            samples,
            Some(path_ulp_s.to_str().unwrap()),
        )
        .unwrap();
    }
}

/// Segment plot data for combining into single HTML
pub struct SegmentPlotData {
    pub segment_idx: usize,
    pub interval_x: (f64, f64),
    pub interval_solve: (f64, f64),
    pub err_x_svg: String,
    pub err_solve_svg: String,
    pub ulp_x_svg: String,
    pub ulp_solve_svg: String,
}
