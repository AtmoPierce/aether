use aether_core::math::Vector;
use aether_ml::linear::LinearRegression;
use aether_viz::{
    CorrelationSeries, PlotConfig, PlotStyle, XYSeries, plot_correlation, plot_series_with_config,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let xs = [
        0.25, 0.75, 1.10, 1.50, 1.90, 2.20, 2.60, 2.95, 3.35, 3.70, 4.10, 4.45,
    ];
    let ys = [
        1.25, 1.95, 2.70, 3.05, 3.90, 4.35, 5.10, 5.55, 6.15, 6.80, 7.45, 7.95,
    ];

    let features: [Vector<f64, 1>; 12] = core::array::from_fn(|index| Vector::new([xs[index]]));
    let mut model = LinearRegression::<f64, 1> {
        weights: Vector::new([0.0]),
        bias: 0.0,
        lr: 0.02,
    };

    for _ in 0..4_000 {
        model.train_epoch(&features, &ys);
    }

    let correlation = pearson_correlation(&xs, &ys);
    let output_dir = format!("{}/outputs", env!("CARGO_MANIFEST_DIR"));
    std::fs::create_dir_all(&output_dir)?;

    let scatter_path = format!("{output_dir}/correlation_scatter.svg");
    let scatter_title = format!("Feature / target correlation (r = {:.3})", correlation);
    plot_correlation(
        &[CorrelationSeries::new(&xs, &ys).with_label("samples")],
        scatter_title.as_str(),
        Some("feature x"),
        Some("target y"),
        Some(scatter_path.as_str()),
    )?;

    let x_min = xs.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let line_xs = (0..128)
        .map(|index| {
            let t = index as f64 / 127.0;
            x_min + (x_max - x_min) * t
        })
        .collect::<Vec<_>>();
    let line_ys = line_xs
        .iter()
        .map(|&x| model.predict(&Vector::new([x])))
        .collect::<Vec<_>>();

    let fit_path = format!("{output_dir}/correlation_with_fit.svg");
    let fit_title = format!(
        "Linear fit overlay (y = {:.3} x + {:.3})",
        model.weights[0],
        model.bias
    );
    let config = PlotConfig::new(fit_title.as_str())
        .with_x_label("feature x")
        .with_y_label("target y")
        .with_save_path(fit_path.as_str());
    let series = [
        XYSeries::new(&xs, &ys)
            .with_label("samples")
            .with_style(PlotStyle::Points),
        XYSeries::new(&line_xs, &line_ys)
            .with_label("linear regression fit")
            .with_style(PlotStyle::Line),
    ];
    plot_series_with_config(&series, &config)?;

    println!("pearson correlation: {:.6}", correlation);
    println!("fitted slope: {:.6}", model.weights[0]);
    println!("fitted intercept: {:.6}", model.bias);
    println!("wrote {}", scatter_path);
    println!("wrote {}", fit_path);

    Ok(())
}

fn pearson_correlation(xs: &[f64], ys: &[f64]) -> f64 {
    assert_eq!(xs.len(), ys.len(), "correlation inputs must have equal length");
    assert!(!xs.is_empty(), "correlation inputs must not be empty");

    let count = xs.len() as f64;
    let mean_x = xs.iter().sum::<f64>() / count;
    let mean_y = ys.iter().sum::<f64>() / count;

    let mut covariance = 0.0;
    let mut variance_x = 0.0;
    let mut variance_y = 0.0;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let dx = x - mean_x;
        let dy = y - mean_y;
        covariance += dx * dy;
        variance_x += dx * dx;
        variance_y += dy * dy;
    }

    covariance / (variance_x.sqrt() * variance_y.sqrt())
}