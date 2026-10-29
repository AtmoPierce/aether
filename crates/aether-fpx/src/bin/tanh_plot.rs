#[cfg(feature = "ocp")]
use std::fs;

#[cfg(feature = "ocp")]
use aether_viz::{plot_series, PlotStyle, XYSeries};
#[cfg(feature = "ocp")]
use aether_fpx::ocp::fp8::{Fp8E4M3, Fp8E5M2};
#[cfg(feature = "ocp")]
use aether_fpx::ocp::trig::{tanh_e4m3, tanh_e5m2};
#[cfg(feature = "ocp")]
use aether_fpx::trig::tanhf;

#[cfg(feature = "ocp")]
#[derive(Clone, Copy)]
struct Sample {
    x: f32,
    core: f32,
    x_q_e4m3: f32,
    x_q_e5m2: f32,
    core_at_x_q_e4m3: f32,
    core_at_x_q_e5m2: f32,
    fpx: f32,
    ocp_e4m3: f32,
    ocp_e5m2: f32,
}

#[cfg(feature = "ocp")]
fn main() {
    let x_min = -4.0f32;
    let x_max = 4.0f32;
    let n = 1601usize;
    let point_stride = 8usize;

    let mut samples = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f32 / (n - 1) as f32;
        let x = x_min + (x_max - x_min) * t;

        let core = x.tanh();
        let x_q_e4m3 = f32::from(Fp8E4M3::from(x));
        let x_q_e5m2 = f32::from(Fp8E5M2::from(x));
        let core_at_x_q_e4m3 = x_q_e4m3.tanh();
        let core_at_x_q_e5m2 = x_q_e5m2.tanh();
        let fpx = tanhf(x);
        let ocp_e4m3 = f32::from(tanh_e4m3(Fp8E4M3::from(x)));
        let ocp_e5m2 = f32::from(tanh_e5m2(Fp8E5M2::from(x)));

        samples.push(Sample {
            x,
            core,
            x_q_e4m3,
            x_q_e5m2,
            core_at_x_q_e4m3,
            core_at_x_q_e5m2,
            fpx,
            ocp_e4m3,
            ocp_e5m2,
        });
    }

    let mut csv = String::from("x,core_tanh_f32,x_q_e4m3,x_q_e5m2,core_tanh_at_x_q_e4m3,core_tanh_at_x_q_e5m2,fpx_tanhf,ocp_tanh_e4m3,ocp_tanh_e5m2\n");
    let mut csv_err = String::from("x,abs_err_fpx,abs_err_ocp_e4m3,abs_err_ocp_e5m2\n");
    let mut xs = Vec::with_capacity(samples.len());
    let mut ys_core = Vec::with_capacity(samples.len());
    let mut ys_core_q_e4 = Vec::with_capacity(samples.len());
    let mut ys_core_q_e5 = Vec::with_capacity(samples.len());
    let mut ys_fpx = Vec::with_capacity(samples.len());
    let mut ys_e4 = Vec::with_capacity(samples.len());
    let mut ys_e5 = Vec::with_capacity(samples.len());
    let mut ys_err_fpx = Vec::with_capacity(samples.len());
    let mut ys_err_e4 = Vec::with_capacity(samples.len());
    let mut ys_err_e5 = Vec::with_capacity(samples.len());

    for s in &samples {
        csv.push_str(&format!(
            "{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            s.x,
            s.core,
            s.x_q_e4m3,
            s.x_q_e5m2,
            s.core_at_x_q_e4m3,
            s.core_at_x_q_e5m2,
            s.fpx,
            s.ocp_e4m3,
            s.ocp_e5m2
        ));

        let err_fpx = (s.fpx - s.core).abs();
        let err_e4 = (s.ocp_e4m3 - s.core).abs();
        let err_e5 = (s.ocp_e5m2 - s.core).abs();

        csv_err.push_str(&format!(
            "{:.9},{:.9},{:.9},{:.9}\n",
            s.x, err_fpx, err_e4, err_e5
        ));

        xs.push(s.x as f64);
        ys_core.push(s.core as f64);
        ys_core_q_e4.push(s.core_at_x_q_e4m3 as f64);
        ys_core_q_e5.push(s.core_at_x_q_e5m2 as f64);
        ys_fpx.push(s.fpx as f64);
        ys_e4.push(s.ocp_e4m3 as f64);
        ys_e5.push(s.ocp_e5m2 as f64);
        ys_err_fpx.push(err_fpx as f64);
        ys_err_e4.push(err_e4 as f64);
        ys_err_e5.push(err_e5 as f64);
    }

    let out_dir = "target/aether-fpx/plots";
    fs::create_dir_all(out_dir).expect("failed to create output directory");

    let csv_path = format!("{}/tanh_compare.csv", out_dir);
    fs::write(&csv_path, csv).expect("failed to write csv");
    let csv_err_path = format!("{}/tanh_compare_abs_error.csv", out_dir);
    fs::write(&csv_err_path, csv_err).expect("failed to write error csv");

    let svg_path = format!("{}/tanh_compare.svg", out_dir);

    let mut xs_points = Vec::with_capacity((xs.len() / point_stride) + 1);
    let mut ys_core_q_e5_points = Vec::with_capacity((ys_core_q_e5.len() / point_stride) + 1);
    let mut ys_core_q_e4_points = Vec::with_capacity((ys_core_q_e4.len() / point_stride) + 1);
    let mut ys_e5_points = Vec::with_capacity((ys_e5.len() / point_stride) + 1);
    let mut ys_e4_points = Vec::with_capacity((ys_e4.len() / point_stride) + 1);

    for i in (0..xs.len()).step_by(point_stride) {
        xs_points.push(xs[i]);
        ys_core_q_e5_points.push(ys_core_q_e5[i]);
        ys_core_q_e4_points.push(ys_core_q_e4[i]);
        ys_e5_points.push(ys_e5[i]);
        ys_e4_points.push(ys_e4[i]);
    }

    let series = [
        XYSeries::new(&xs, &ys_core)
            .with_label("core f32::tanh")
            .with_style(PlotStyle::Line),
        XYSeries::new(&xs_points, &ys_core_q_e5_points)
            .with_label("core tanh(fp8 e5m2-quantized x)")
            .with_style(PlotStyle::Points),
        XYSeries::new(&xs_points, &ys_core_q_e4_points)
            .with_label("core tanh(fp8 e4m3-quantized x)")
            .with_style(PlotStyle::Points),
        XYSeries::new(&xs, &ys_fpx)
            .with_label("aether-fpx::trig::tanhf")
            .with_style(PlotStyle::Line),
        XYSeries::new(&xs_points, &ys_e5_points)
            .with_label("ocp::trig::tanh_e5m2")
            .with_style(PlotStyle::Points),
        XYSeries::new(&xs_points, &ys_e4_points)
            .with_label("ocp::trig::tanh_e4m3")
            .with_style(PlotStyle::Points),
    ];

    plot_series(
        &series,
        "tanh: continuous core + quantized-input core vs fpx and OCP",
        Some("x"),
        Some("tanh(x)"),
        Some(&svg_path),
    )
    .expect("failed to render svg via aether_viz");

    let svg_err_path = format!("{}/tanh_compare_abs_error.svg", out_dir);
    let err_series = [
        XYSeries::new(&xs, &ys_err_fpx)
            .with_label("|aether-fpx::trig::tanhf - core f32::tanh|")
            .with_style(PlotStyle::Line),
        XYSeries::new(&xs, &ys_err_e5)
            .with_label("|ocp::trig::tanh_e5m2 - core f32::tanh|")
            .with_style(PlotStyle::Line),
        XYSeries::new(&xs, &ys_err_e4)
            .with_label("|ocp::trig::tanh_e4m3 - core f32::tanh|")
            .with_style(PlotStyle::Line),
    ];

    plot_series(
        &err_series,
        "absolute error vs core f32::tanh",
        Some("x"),
        Some("abs error"),
        Some(&svg_err_path),
    )
    .expect("failed to render error svg via aether_viz");

    println!("Wrote {}", csv_path);
    println!("Wrote {}", svg_path);
    println!("Wrote {}", csv_err_path);
    println!("Wrote {}", svg_err_path);
}

#[cfg(not(feature = "ocp"))]
fn main() {
    eprintln!(
        "Run with: cargo run --manifest-path crates/aether-fpx/Cargo.toml --features ocp --bin tanh_plot"
    );
}
