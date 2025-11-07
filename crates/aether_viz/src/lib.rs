/// Plot type selector for a series.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlotStyle {
    /// Continuous line
    Line,
    /// Discrete points
    Points,
}

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
    if min_e == max_e { max_e = min_e + 1.0; }

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
            chart.configure_mesh()
                .x_desc("x").y_desc("error")
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
            chart.configure_mesh()
                .x_desc("x").y_desc("error")
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
    assert!(!series.is_empty(), "plot_series: need at least one series");

    // global bounds
    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;
    for s in series {
        assert_eq!(s.xs.len(), s.ys.len(), "x/y length mismatch in a series");
        for &x in s.xs { x_min = x_min.min(x); x_max = x_max.max(x); }
        for &y in s.ys { y_min = y_min.min(y); y_max = y_max.max(y); }
    }
    if x_min == x_max { x_max = x_min + 1.0; }
    if y_min == y_max { y_max = y_min + 1.0; }

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
                .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

            let mut mesh = chart.configure_mesh();
            if let Some(lbl) = x_label { mesh.x_desc(lbl); }
            if let Some(lbl) = y_label { mesh.y_desc(lbl); }
            mesh.label_style(("monospace", 18))
                .axis_desc_style(("monospace", 18))
                .draw()?;

            for (idx, s) in series.iter().enumerate() {
                let color = Palette99::pick(idx).mix(1.0);

                match s.style {
                    PlotStyle::Line => {
                        let mut ds = chart.draw_series(LineSeries::new(
                            s.xs.iter().zip(s.ys).map(|(&x, &y)| (x, y)),
                            &color,
                        ))?;
                        if let Some(lbl) = s.label {
                            ds = ds
                                .label(lbl)
                                .legend(move |(x, y)| {
                                    PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))
                                });
                        }
                    }
                    PlotStyle::Points => {
                        let mut ds = chart.draw_series(
                            s.xs.iter().zip(s.ys).map(|(&x, &y)| Circle::new((x, y), 3, color.filled())),
                        )?;
                        if let Some(lbl) = s.label {
                            ds = ds
                                .label(lbl)
                                .legend(move |(x, y)| {
                                    Circle::new((x + 10, y), 4, color.filled())
                                });
                        }
                    }
                }
            }

            if series.iter().any(|s| s.label.is_some()) {
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
                .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

            let mut mesh = chart.configure_mesh();
            if let Some(lbl) = x_label { mesh.x_desc(lbl); }
            if let Some(lbl) = y_label { mesh.y_desc(lbl); }
            mesh.label_style(("monospace", 18))
                .axis_desc_style(("monospace", 18))
                .draw()?;

            for (idx, s) in series.iter().enumerate() {
                let color = Palette99::pick(idx).mix(1.0);

                match s.style {
                    PlotStyle::Line => {
                        let mut ds = chart.draw_series(LineSeries::new(
                            s.xs.iter().zip(s.ys).map(|(&x, &y)| (x, y)),
                            &color,
                        ))?;
                        if let Some(lbl) = s.label {
                            ds = ds
                                .label(lbl)
                                .legend(move |(x, y)| {
                                    PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))
                                });
                        }
                    }
                    PlotStyle::Points => {
                        let mut ds = chart.draw_series(
                            s.xs.iter().zip(s.ys).map(|(&x, &y)| Circle::new((x, y), 3, color.filled())),
                        )?;
                        if let Some(lbl) = s.label {
                            ds = ds
                                .label(lbl)
                                .legend(move |(x, y)| {
                                    Circle::new((x + 10, y), 4, color.filled())
                                });
                        }
                    }
                }
            }

            if series.iter().any(|s| s.label.is_some()) {
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
