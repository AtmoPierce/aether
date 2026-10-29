use crate::remez_create;
use crate::remez_piecewise;
use crate::design::remez::RemezApprox;
use crate::design::functions::{sqrt_f, sqrt_df}; // the truth f64 wrappers above

// Single-interval generator (probably just for testing at this point)
remez_create!{
    name        = remez_sqrt_create,
    N           = 23,
    K           = 24,
    func        = sqrt_f,
    deriv       = sqrt_df,
    interval    = [0.0, 1000.0],
    max_error   = 7.0,
    samples     = 100,
    title       = "Remez sqrt(x) on [0, 1000]",
    file_prefix = "remez_sqrt",
    boundary_tol = 1.0,
    enforce_endpoints = false
}

// Piecewise generators
remez_piecewise! {
    name  = remez_sqrt_piecewise,
    func  = sqrt_f,
    deriv = sqrt_df,
    pieces = [
        {
            seg_name        = remez_sqrt_lo,
            N               = 18,
            K               = 19,
            interval        = [0.00, 0.5],
            max_error       = 5e-4,
            samples         = 8192,
            file_prefix     = "sqrt_lo",
            boundary_tol    = 5e-4,
            enforce_endpoints = true
        },
        {
            seg_name        = remez_sqrt_mid,
            N               = 18,
            K               = 19,
            interval        = [0.5, 2.0],
            max_error       = 1e-4,
            samples         = 8192,
            file_prefix     = "sqrt_mid",
            boundary_tol    = 1e-4,
            enforce_endpoints = true
        },
        {
            seg_name        = remez_sqrt_hi,
            N               = 18,
            K               = 19,
            interval        = [2.0, 4.0],
            max_error       = 1e-4,
            samples         = 8192,
            file_prefix     = "sqrt_hi",
            boundary_tol    = 1e-4,
            enforce_endpoints = true
        }
    ]
}

type RzLo = RemezApprox<f64, 18, 19>;
type RzMid = RemezApprox<f64, 18, 19>;
type RzHi = RemezApprox<f64, 18, 19>;

use aether_viz::plot_error_on_grid;

pub fn remez_sqrt_piecewise_merged_plot() {
    // [0.0, 0.5]
    let mut rz_lo = RzLo::new(0.0, 0.5);
    rz_lo.solve(sqrt_f, sqrt_df);

    let mut rz_mid = RzMid::new(0.5, 2.0);
    rz_mid.solve(sqrt_f, sqrt_df);

    let mut rz_hi = RzHi::new(2.0, 4.0);
    rz_hi.solve(sqrt_f, sqrt_df);

    let approx = |x: f64| -> f64 {
        if x < 0.5 {
            return rz_lo.eval_poly(x)
        } 
        if x < 2.0 && x >= 0.5 {
            return rz_mid.eval_poly(x)
        } else {
            return rz_hi.eval_poly(x)
        }
    };

    std::fs::create_dir_all("target/plots").ok();
    let err_file = "target/plots/sqrt_piecewise_error.svg";

    plot_error_on_grid(
        approx,
        sqrt_f,
        0.0,
        4.0,
        8192,
        Some(err_file),
    )
    .unwrap();
}

pub fn remez_sqrt_piecewise_fx_plot() {
    use aether_viz::{XYSeries, PlotStyle, plot_series};
    use std::fs;

    // Solve segments

    let mut rz_lo  = RzLo::new(0.0, 0.5); rz_lo.solve(sqrt_f, sqrt_df);
    let mut rz_mid = RzMid::new(0.5, 2.0); rz_mid.solve(sqrt_f, sqrt_df);
    let mut rz_hi  = RzHi::new(2.0, 4.0);  rz_hi.solve(sqrt_f, sqrt_df);

    let approx = |x: f64| -> f64 {
        if x < 0.5 {
            return rz_lo.eval_poly(x)
        } 
        if x < 2.0 && x >= 0.5 {
            return rz_mid.eval_poly(x)
        } else {
            return rz_hi.eval_poly(x)
        }
    };

    fs::create_dir_all("target/plots").ok();
    let out_file = "target/plots/sqrt_piecewise_fx_vs_poly.svg";

    // Sampling grid
    let n = 3000;
    let a = 0.0;
    let b = 4.0;
    let xs: Vec<f64> = (0..=n).map(|i| a + (b - a)*(i as f64 / n as f64)).collect();

    // curves
    let y_true: Vec<f64> = xs.iter().map(|&x| sqrt_f(x)).collect();
    let y_piecewise: Vec<f64> = xs.iter().map(|&x| approx(x)).collect();
    let y_lo:  Vec<f64> = xs.iter().map(|&x| if x<0.5 { rz_lo.eval_poly(x) } else { f64::NAN }).collect();
    let y_mid: Vec<f64> = xs.iter().map(|&x| if (0.5..2.0).contains(&x) { rz_mid.eval_poly(x) } else { f64::NAN }).collect();
    let y_hi:  Vec<f64> = xs.iter().map(|&x| if x>=2.0 { rz_hi.eval_poly(x) } else { f64::NAN }).collect();

    let series = [
        XYSeries { xs: &xs, ys: &y_true,       label: Some("sqrt_f (true)"),     style: PlotStyle::Line },
        XYSeries { xs: &xs, ys: &y_piecewise,  label: Some("piecewise P(x)"),    style: PlotStyle::Line },
        XYSeries { xs: &xs, ys: &y_lo,         label: Some("P_lo(x)"),           style: PlotStyle::Points },
        XYSeries { xs: &xs, ys: &y_mid,        label: Some("P_mid(x)"),          style: PlotStyle::Points },
        XYSeries { xs: &xs, ys: &y_hi,         label: Some("P_hi(x)"),           style: PlotStyle::Points },
    ];

    plot_series(
        &series,
        "Piecewise Remez sqrt(x) Approximation",
        Some("x"),
        Some("y"),
        Some(out_file),
    ).unwrap();
}
