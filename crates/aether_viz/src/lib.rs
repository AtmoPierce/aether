/// Plot type selector for a series.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlotStyle {
    /// Continuous line
    Line,
    /// Discrete points
    Points,
}

/// Axis scaling for plots.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxisScale {
    /// Plot values directly.
    Linear,
    /// Plot base-10 logarithm of strictly positive values.
    Log10,
}

/// Configuration for XY-series plots.
#[derive(Debug, Clone, Copy)]
pub struct PlotConfig<'a> {
    pub title: &'a str,
    pub x_label: Option<&'a str>,
    pub y_label: Option<&'a str>,
    pub save_path: Option<&'a str>,
    pub x_scale: AxisScale,
    pub y_scale: AxisScale,
}

impl<'a> PlotConfig<'a> {
    pub fn new(title: &'a str) -> Self {
        Self {
            title,
            x_label: None,
            y_label: None,
            save_path: None,
            x_scale: AxisScale::Linear,
            y_scale: AxisScale::Linear,
        }
    }

    pub fn with_x_label(mut self, label: &'a str) -> Self {
        self.x_label = Some(label);
        self
    }

    pub fn with_y_label(mut self, label: &'a str) -> Self {
        self.y_label = Some(label);
        self
    }

    pub fn with_save_path(mut self, save_path: &'a str) -> Self {
        self.save_path = Some(save_path);
        self
    }

    pub fn with_x_scale(mut self, scale: AxisScale) -> Self {
        self.x_scale = scale;
        self
    }

    pub fn with_y_scale(mut self, scale: AxisScale) -> Self {
        self.y_scale = scale;
        self
    }
}

pub mod realtime;
pub use realtime::{CsvRealtimeConfig, CsvRealtimePlotter};

/// One (x,y) series to draw.
#[derive(Debug, Clone, Copy)]
pub struct XYSeries<'a> {
    pub xs: &'a [f64],
    pub ys: &'a [f64],
    pub label: Option<&'a str>,
    pub style: PlotStyle,
}

impl<'a> XYSeries<'a> {
    pub fn new(xs: &'a [f64], ys: &'a [f64]) -> Self {
        Self { xs, ys, label: None, style: PlotStyle::Line }
    }
    pub fn with_label(mut self, label: &'a str) -> Self { self.label = Some(label); self }
    pub fn with_style(mut self, style: PlotStyle) -> Self { self.style = style; self }
}

/// One histogram series to draw.
#[derive(Debug, Clone, Copy)]
pub struct HistogramSeries<'a> {
    pub values: &'a [f64],
    pub label: Option<&'a str>,
    pub bins: usize,
}

impl<'a> HistogramSeries<'a> {
    pub fn new(values: &'a [f64]) -> Self {
        Self { values, label: None, bins: 40 }
    }
    pub fn with_label(mut self, label: &'a str) -> Self { self.label = Some(label); self }
    pub fn with_bins(mut self, bins: usize) -> Self { self.bins = bins.max(1); self }
}

/// One paired-data series for a correlation plot.
#[derive(Debug, Clone, Copy)]
pub struct CorrelationSeries<'a> {
    pub xs: &'a [f64],
    pub ys: &'a [f64],
    pub label: Option<&'a str>,
}

impl<'a> CorrelationSeries<'a> {
    pub fn new(xs: &'a [f64], ys: &'a [f64]) -> Self {
        Self { xs, ys, label: None }
    }

    pub fn with_label(mut self, label: &'a str) -> Self {
        self.label = Some(label);
        self
    }
}

#[cfg(feature = "svg")]
use plotters::prelude::SVGBackend;
#[cfg(feature = "bitmap")]
use plotters::prelude::BitMapBackend;

use plotters::coord::Shift;
use plotters::prelude::*;
use plotters::series::LineSeries;

#[cfg(feature = "svg")]
fn with_svg_root<F>(
    save_path: Option<&str>,
    size: (u32, u32),
    draw: F,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: for<'a> FnOnce(DrawingArea<SVGBackend<'a>, Shift>) -> Result<(), Box<dyn std::error::Error>>,
{
    if let Some(path) = save_path {
        let area = SVGBackend::new(path, size).into_drawing_area();
        draw(area)?;
    } else {
        let mut sink = String::new();
        let area = SVGBackend::with_string(&mut sink, size).into_drawing_area();
        draw(area)?;
    }
    Ok(())
}

#[cfg(feature = "bitmap")]
fn with_bitmap_root<F>(
    save_path: Option<&str>,
    size: (u32, u32),
    draw: F,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: FnOnce(DrawingArea<BitMapBackend<'_>, Shift>) -> Result<(), Box<dyn std::error::Error>>,
{
    if let Some(path) = save_path {
        let area = BitMapBackend::new(path, size).into_drawing_area();
        draw(area)?;
    } else {
        let mut buf = vec![0u8; (size.0 as usize) * (size.1 as usize) * 3];
        let area = BitMapBackend::with_buffer(&mut buf, size).into_drawing_area();
        draw(area)?;
    }
    Ok(())
}

/// Plot e(x) = f_approx(x) - f_true(x) on [a, b].
pub fn plot_error_on_grid<F, G>(
    f_approx: F,
    f_true: G,
    a: f64,
    b: f64,
    samples: usize,
    save_path: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: Fn(f64) -> f64,
    G: Fn(f64) -> f64,
{
    let size = (1200, 400);

    let mut min_e = f64::INFINITY;
    let mut max_e = f64::NEG_INFINITY;
    let mut pts = Vec::with_capacity(samples + 1);
    for i in 0..=samples {
        let t = i as f64 / samples as f64;
        let x = a + (b - a) * t;
        let e = f_approx(x) - f_true(x);
        min_e = min_e.min(e);
        max_e = max_e.max(e);
        pts.push((x, e));
    }
    if min_e == max_e {
        let pad = if min_e == 0.0 {
            1.0
        } else {
            min_e.abs() * 0.1
        };
        min_e -= pad;
        max_e += pad;
    }

    #[cfg(feature = "svg")]
    {
        return with_svg_root(save_path, size, |root| {
            root.fill(&WHITE)?;
            let mut chart = ChartBuilder::on(&root)
                .margin(10)
                .caption("Error: f_approx(x) - f_true(x)", ("monospace", 20))
                .x_label_area_size(35)
                .y_label_area_size(60)
                .build_cartesian_2d(a..b, min_e..max_e)?;
            let y_label_fmt = |v: &f64| format_axis_value(*v);
            chart.configure_mesh()
                .x_desc("x")
                .y_desc("error")
                .y_label_formatter(&y_label_fmt)
                .label_style(("monospace", 16))
                .axis_desc_style(("monospace", 16))
                .draw()?;
            chart.draw_series(LineSeries::new(pts, &BLACK))?;
            if save_path.is_some() { root.present()?; }
            Ok(())
        });
    }

    #[cfg(all(not(feature = "svg"), feature = "bitmap"))]
    {
        return with_bitmap_root(save_path, size, |root| {
            root.fill(&WHITE)?;
            let mut chart = ChartBuilder::on(&root)
                .margin(10)
                .caption("Error: f_approx(x) - f_true(x)", ("monospace", 20))
                .x_label_area_size(35)
                .y_label_area_size(60)
                .build_cartesian_2d(a..b, min_e..max_e)?;
            let y_label_fmt = |v: &f64| format_axis_value(*v);
            chart.configure_mesh()
                .x_desc("x")
                .y_desc("error")
                .y_label_formatter(&y_label_fmt)
                .label_style(("monospace", 16))
                .axis_desc_style(("monospace", 16))
                .draw()?;
            chart.draw_series(LineSeries::new(pts, &BLACK))?;
            if save_path.is_some() { root.present()?; }
            Ok(())
        });
    }

    #[allow(unreachable_code)]
    Err("no backend enabled; enable `svg` or `bitmap` feature".into())
}

/// Plot multiple series on one chart with colored legend markers.
pub fn plot_series(
    series: &[XYSeries<'_>],
    title: &str,
    x_label: Option<&str>,
    y_label: Option<&str>,
    save_path: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut config = PlotConfig::new(title);
    if let Some(label) = x_label {
        config = config.with_x_label(label);
    }
    if let Some(label) = y_label {
        config = config.with_y_label(label);
    }
    if let Some(path) = save_path {
        config = config.with_save_path(path);
    }
    plot_series_with_config(series, &config)
}

/// Plot multiple series with optional per-axis scaling.
pub fn plot_series_scaled(
    series: &[XYSeries<'_>],
    title: &str,
    x_label: Option<&str>,
    y_label: Option<&str>,
    save_path: Option<&str>,
    x_scale: AxisScale,
    y_scale: AxisScale,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut config = PlotConfig::new(title)
        .with_x_scale(x_scale)
        .with_y_scale(y_scale);
    if let Some(label) = x_label {
        config = config.with_x_label(label);
    }
    if let Some(label) = y_label {
        config = config.with_y_label(label);
    }
    if let Some(path) = save_path {
        config = config.with_save_path(path);
    }
    plot_series_with_config(series, &config)
}

/// Plot multiple series using a reusable plot configuration.
pub fn plot_series_with_config(
    series: &[XYSeries<'_>],
    config: &PlotConfig<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    assert!(!series.is_empty(), "plot_series: need at least one series");

    struct PreparedSeries<'a> {
        points: Vec<(f64, f64)>,
        label: Option<&'a str>,
        style: PlotStyle,
    }

    let prepared = series
        .iter()
        .map(|s| {
            assert_eq!(s.xs.len(), s.ys.len(), "x/y length mismatch in a series");

            let points = s
                .xs
                .iter()
                .zip(s.ys.iter())
                .filter_map(|(&x, &y)| {
                    transform_point(x, y, config.x_scale, config.y_scale)
                })
                .collect::<Vec<_>>();

            PreparedSeries {
                points,
                label: s.label,
                style: s.style,
            }
        })
        .filter(|s| !s.points.is_empty())
        .collect::<Vec<_>>();

    if prepared.is_empty() {
        return Err("plot_series: no finite points to plot".into());
    }

    // global bounds
    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;
    for s in &prepared {
        for &(x, y) in &s.points {
            x_min = x_min.min(x);
            x_max = x_max.max(x);
            y_min = y_min.min(y);
            y_max = y_max.max(y);
        }
    }

    if x_min == x_max { x_max = x_min + 1.0; }
    if y_min == y_max { y_max = y_min + 1.0; }

    let x_pad = (x_max - x_min).abs() * 0.03;
    let y_pad = (y_max - y_min).abs() * 0.03;
    let x_range = (x_min - x_pad)..(x_max + x_pad);
    let y_range = (y_min - y_pad)..(y_max + y_pad);
    let x_desc = axis_label(config.x_label, config.x_scale);
    let y_desc = axis_label(config.y_label, config.y_scale);

    let size = (1920, 1080);

    #[cfg(feature = "svg")]
    {
        return with_svg_root(config.save_path, size, |root| {
            root.fill(&WHITE)?;

            let mut chart = ChartBuilder::on(&root)
                .margin(20)
            .caption(config.title, ("monospace", 24))
                .x_label_area_size(40)
                .y_label_area_size(80)
                .build_cartesian_2d(x_range.clone(), y_range.clone())?;

            let y_label_fmt = |v: &f64| format_axis_value(*v);
            let mut mesh = chart.configure_mesh();
            if let Some(ref lbl) = x_desc { mesh.x_desc(lbl); }
            if let Some(ref lbl) = y_desc { mesh.y_desc(lbl); }
            mesh.y_label_formatter(&y_label_fmt)
                .label_style(("monospace", 18))
                .axis_desc_style(("monospace", 18))
                .draw()?;

            for (idx, s) in prepared.iter().enumerate() {
                let color = Palette99::pick(idx).mix(1.0);

                match s.style {
                    PlotStyle::Line => {
                        let ds = chart.draw_series(LineSeries::new(
                            s.points.iter().copied(),
                            &color,
                        ))?;
                        if let Some(lbl) = s.label {
                            let _ = ds
                                .label(lbl)
                                .legend(move |(x, y)| {
                                    PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))
                                });
                        }
                    }
                    PlotStyle::Points => {
                        let ds = chart.draw_series(
                            s.points.iter().map(|&(x, y)| Circle::new((x, y), 3, color.filled())),
                        )?;
                        if let Some(lbl) = s.label {
                            let _ = ds
                                .label(lbl)
                                .legend(move |(x, y)| {
                                    Circle::new((x + 10, y), 4, color.filled())
                                });
                        }
                    }
                }
            }

            if prepared.iter().any(|s| s.label.is_some()) {
                chart.configure_series_labels()
                    .border_style(&BLACK)
                    .background_style(WHITE.mix(0.8))
                    .label_font(("monospace", 16))
                    .draw()?;
            }

            if config.save_path.is_some() { root.present()?; }
            Ok(())
        });
    }

    #[cfg(all(not(feature = "svg"), feature = "bitmap"))]
    {
        return with_bitmap_root(config.save_path, size, |root| {
            root.fill(&WHITE)?;

            let mut chart = ChartBuilder::on(&root)
                .margin(20)
                .caption(config.title, ("monospace", 24))
                .x_label_area_size(40)
                .y_label_area_size(80)
                .build_cartesian_2d(x_range.clone(), y_range.clone())?;

            let y_label_fmt = |v: &f64| format_axis_value(*v);
            let mut mesh = chart.configure_mesh();
            if let Some(ref lbl) = x_desc { mesh.x_desc(lbl); }
            if let Some(ref lbl) = y_desc { mesh.y_desc(lbl); }
            mesh.y_label_formatter(&y_label_fmt)
                .label_style(("monospace", 18))
                .axis_desc_style(("monospace", 18))
                .draw()?;

            for (idx, s) in prepared.iter().enumerate() {
                let color = Palette99::pick(idx).mix(1.0);

                match s.style {
                    PlotStyle::Line => {
                        let ds = chart.draw_series(LineSeries::new(
                            s.points.iter().copied(),
                            &color,
                        ))?;
                        if let Some(lbl) = s.label {
                            let _ = ds
                                .label(lbl)
                                .legend(move |(x, y)| {
                                    PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))
                                });
                        }
                    }
                    PlotStyle::Points => {
                        let ds = chart.draw_series(
                            s.points.iter().map(|&(x, y)| Circle::new((x, y), 3, color.filled())),
                        )?;
                        if let Some(lbl) = s.label {
                            let _ = ds
                                .label(lbl)
                                .legend(move |(x, y)| {
                                    Circle::new((x + 10, y), 4, color.filled())
                                });
                        }
                    }
                }
            }

            if prepared.iter().any(|s| s.label.is_some()) {
                chart.configure_series_labels()
                    .border_style(&BLACK)
                    .background_style(WHITE.mix(0.8))
                    .label_font(("monospace", 16))
                    .draw()?;
            }

            if config.save_path.is_some() { root.present()?; }
            Ok(())
        });
    }

    #[allow(unreachable_code)]
    Err("no backend enabled; enable `svg` or `bitmap` feature".into())
}

/// Plot one or more histograms using shared x-axis bounds.
pub fn plot_histograms(
    series: &[HistogramSeries<'_>],
    title: &str,
    x_label: Option<&str>,
    y_label: Option<&str>,
    save_path: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    assert!(!series.is_empty(), "plot_histograms: need at least one series");

    struct PreparedHistogram<'a> {
        counts: Vec<f64>,
        label: Option<&'a str>,
    }

    let finite_values = series
        .iter()
        .flat_map(|s| s.values.iter().copied())
        .filter(|value| value.is_finite())
        .collect::<Vec<_>>();

    if finite_values.is_empty() {
        return Err("plot_histograms: no finite values to plot".into());
    }

    let mut x_min = finite_values.iter().copied().fold(f64::INFINITY, f64::min);
    let mut x_max = finite_values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if x_min == x_max {
        x_min -= 0.5;
        x_max += 0.5;
    }

    let bin_count = series
        .iter()
        .map(|s| s.bins.max(1))
        .max()
        .unwrap_or(1);
    let bin_width = (x_max - x_min) / bin_count as f64;

    let prepared = series
        .iter()
        .map(|s| {
            let mut counts = vec![0.0_f64; bin_count];
            for value in s.values.iter().copied().filter(|value| value.is_finite()) {
                let mut index = ((value - x_min) / bin_width).floor() as usize;
                if index >= bin_count {
                    index = bin_count - 1;
                }
                counts[index] += 1.0;
            }

            PreparedHistogram {
                counts,
                label: s.label,
            }
        })
        .collect::<Vec<_>>();

    let y_max = prepared
        .iter()
        .flat_map(|s| s.counts.iter().copied())
        .fold(0.0_f64, f64::max)
        .max(1.0);

    let size = (1920, 1080);

    #[cfg(feature = "svg")]
    {
        return with_svg_root(save_path, size, |root| {
            root.fill(&WHITE)?;

            let mut chart = ChartBuilder::on(&root)
                .margin(20)
                .caption(title, ("monospace", 24))
                .x_label_area_size(40)
                .y_label_area_size(80)
                .build_cartesian_2d(x_min..x_max, 0.0_f64..(y_max * 1.1))?;

            let y_label_fmt = |v: &f64| format_axis_value(*v);
            let mut mesh = chart.configure_mesh();
            if let Some(lbl) = x_label { mesh.x_desc(lbl); }
            if let Some(lbl) = y_label { mesh.y_desc(lbl); }
            mesh.y_label_formatter(&y_label_fmt)
                .label_style(("monospace", 18))
                .axis_desc_style(("monospace", 18))
                .draw()?;

            for (idx, histogram) in prepared.iter().enumerate() {
                let fill = Palette99::pick(idx).mix(0.35).filled();
                let stroke = Palette99::pick(idx).mix(1.0).stroke_width(1);
                let ds = chart.draw_series((0..bin_count).map(|bin_index| {
                    let left = x_min + bin_index as f64 * bin_width;
                    let right = left + bin_width;
                    Rectangle::new([(left, 0.0_f64), (right, histogram.counts[bin_index])], fill)
                }))?;

                if let Some(lbl) = histogram.label {
                    let legend_fill = Palette99::pick(idx).mix(0.35).filled();
                    let _ = ds.label(lbl).legend(move |(x, y)| {
                        Rectangle::new([(x, y - 4), (x + 14, y + 4)], legend_fill)
                    });
                }

                chart.draw_series((0..bin_count).map(|bin_index| {
                    let left = x_min + bin_index as f64 * bin_width;
                    let right = left + bin_width;
                    Rectangle::new([(left, 0.0_f64), (right, histogram.counts[bin_index])], stroke)
                }))?;
            }

            if prepared.iter().any(|s| s.label.is_some()) {
                chart.configure_series_labels()
                    .border_style(&BLACK)
                    .background_style(WHITE.mix(0.8))
                    .label_font(("monospace", 16))
                    .draw()?;
            }

            if save_path.is_some() { root.present()?; }
            Ok(())
        });
    }

    #[cfg(all(not(feature = "svg"), feature = "bitmap"))]
    {
        return with_bitmap_root(save_path, size, |root| {
            root.fill(&WHITE)?;

            let mut chart = ChartBuilder::on(&root)
                .margin(20)
                .caption(title, ("monospace", 24))
                .x_label_area_size(40)
                .y_label_area_size(80)
                .build_cartesian_2d(x_min..x_max, 0.0_f64..(y_max * 1.1))?;

            let y_label_fmt = |v: &f64| format_axis_value(*v);
            let mut mesh = chart.configure_mesh();
            if let Some(lbl) = x_label { mesh.x_desc(lbl); }
            if let Some(lbl) = y_label { mesh.y_desc(lbl); }
            mesh.y_label_formatter(&y_label_fmt)
                .label_style(("monospace", 18))
                .axis_desc_style(("monospace", 18))
                .draw()?;

            for (idx, histogram) in prepared.iter().enumerate() {
                let fill = Palette99::pick(idx).mix(0.35).filled();
                let stroke = Palette99::pick(idx).mix(1.0).stroke_width(1);
                let mut ds = chart.draw_series((0..bin_count).map(|bin_index| {
                    let left = x_min + bin_index as f64 * bin_width;
                    let right = left + bin_width;
                    Rectangle::new([(left, 0.0_f64), (right, histogram.counts[bin_index])], fill)
                }))?;

                if let Some(lbl) = histogram.label {
                    let legend_fill = Palette99::pick(idx).mix(0.35).filled();
                    let _ = ds.label(lbl).legend(move |(x, y)| {
                        Rectangle::new([(x, y - 4), (x + 14, y + 4)], legend_fill)
                    });
                }

                chart.draw_series((0..bin_count).map(|bin_index| {
                    let left = x_min + bin_index as f64 * bin_width;
                    let right = left + bin_width;
                    Rectangle::new([(left, 0.0_f64), (right, histogram.counts[bin_index])], stroke)
                }))?;
            }

            if prepared.iter().any(|s| s.label.is_some()) {
                chart.configure_series_labels()
                    .border_style(&BLACK)
                    .background_style(WHITE.mix(0.8))
                    .label_font(("monospace", 16))
                    .draw()?;
            }

            if save_path.is_some() { root.present()?; }
            Ok(())
        });
    }

    #[allow(unreachable_code)]
    Err("no backend enabled; enable `svg` or `bitmap` feature".into())
}

/// Plot one or more correlation datasets as scatter plots.
pub fn plot_correlation(
    series: &[CorrelationSeries<'_>],
    title: &str,
    x_label: Option<&str>,
    y_label: Option<&str>,
    save_path: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut config = PlotConfig::new(title);
    if let Some(label) = x_label {
        config = config.with_x_label(label);
    }
    if let Some(label) = y_label {
        config = config.with_y_label(label);
    }
    if let Some(path) = save_path {
        config = config.with_save_path(path);
    }
    plot_correlation_with_config(series, &config)
}

/// Plot one or more correlation datasets using a reusable plot configuration.
pub fn plot_correlation_with_config(
    series: &[CorrelationSeries<'_>],
    config: &PlotConfig<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    assert!(!series.is_empty(), "plot_correlation: need at least one series");

    let converted = series
        .iter()
        .map(|s| {
            let mut xy = XYSeries::new(s.xs, s.ys).with_style(PlotStyle::Points);
            if let Some(label) = s.label {
                xy = xy.with_label(label);
            }
            xy
        })
        .collect::<Vec<_>>();

    plot_series_with_config(&converted, config)
}

pub fn plot(
    xs: &[f64],
    ys: &[f64],
    title: &str,
    x_label: Option<&str>,
    y_label: Option<&str>,
    save_path: Option<&str>,
    style: PlotStyle,
) -> Result<(), Box<dyn std::error::Error>> {
    let s = XYSeries { xs, ys, label: None, style };
    plot_series(&[s], title, x_label, y_label, save_path)
}

pub fn plot_psd(
    series: &[XYSeries<'_>],
    title: &str,
    save_path: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut config = PlotConfig::new(title)
        .with_x_label("frequency [Hz]")
        .with_y_label("power spectral density")
        .with_y_scale(AxisScale::Log10);
    if let Some(path) = save_path {
        config = config.with_save_path(path);
    }
    plot_series_with_config(series, &config)
}

fn format_axis_value(v: f64) -> String {
    if !v.is_finite() {
        return format!("{v}");
    }

    if v == 0.0 {
        return "0".to_string();
    }

    let av = v.abs();

    if av < 1.0e-3 || av >= 1.0e4 {
        format!("{:.3e}", v)
    } else {
        let s = format!("{:.6}", v);
        let s = s.trim_end_matches('0').trim_end_matches('.');
        if s == "-0" { "0".to_string() } else { s.to_string() }
    }
}

fn transform_point(
    x: f64,
    y: f64,
    x_scale: AxisScale,
    y_scale: AxisScale,
) -> Option<(f64, f64)> {
    let x = transform_value(x, x_scale)?;
    let y = transform_value(y, y_scale)?;
    if x.is_finite() && y.is_finite() {
        Some((x, y))
    } else {
        None
    }
}

fn transform_value(value: f64, scale: AxisScale) -> Option<f64> {
    match scale {
        AxisScale::Linear => value.is_finite().then_some(value),
        AxisScale::Log10 => (value.is_finite() && value > 0.0).then(|| value.log10()),
    }
}

fn axis_label(label: Option<&str>, scale: AxisScale) -> Option<String> {
    label.map(|label| match scale {
        AxisScale::Linear => label.to_string(),
        AxisScale::Log10 => format!("log10({label})"),
    })
}