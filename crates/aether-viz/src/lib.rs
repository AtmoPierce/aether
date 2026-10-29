/// Plot type selector for a series.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlotStyle {
    /// Continuous line
    Line,
    /// Continuous line using an explicit `Palette99` color index.
    LineWithPalette(usize),
    /// Dashed line using an explicit `Palette99` color index.
    DashedLineWithPalette(usize),
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
    pub size: (u32, u32),
    pub font_scale: f64,
    pub legend_in_header: bool,
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
            size: (1920, 1080),
            font_scale: 1.0,
            legend_in_header: false,
        }
    }

    /// Use thesis-scale typography and place the legend above the data.
    pub fn publication(title: &'a str) -> Self {
        Self::new(title)
            .with_size(1200, 675)
            .with_font_scale(1.4)
            .with_legend_in_header(true)
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

    pub fn with_size(mut self, width: u32, height: u32) -> Self {
        assert!(width > 0 && height > 0, "plot dimensions must be nonzero");
        self.size = (width, height);
        self
    }

    pub fn with_font_scale(mut self, scale: f64) -> Self {
        assert!(
            scale.is_finite() && scale > 0.0,
            "plot font scale must be positive and finite"
        );
        self.font_scale = scale;
        self
    }

    pub fn with_legend_in_header(mut self, enabled: bool) -> Self {
        self.legend_in_header = enabled;
        self
    }

    fn assert_valid(&self) {
        assert!(
            self.size.0 > 0 && self.size.1 > 0,
            "plot dimensions must be nonzero"
        );
        assert!(
            self.font_scale.is_finite() && self.font_scale > 0.0,
            "plot font scale must be positive and finite"
        );
    }
}

pub mod realtime;
pub use realtime::{CsvRealtimeConfig, CsvRealtimePlotter};

pub mod rtos;

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
        Self {
            xs,
            ys,
            label: None,
            style: PlotStyle::Line,
        }
    }
    pub fn with_label(mut self, label: &'a str) -> Self {
        self.label = Some(label);
        self
    }
    pub fn with_style(mut self, style: PlotStyle) -> Self {
        self.style = style;
        self
    }
}

/// One histogram series to draw.
#[derive(Debug, Clone, Copy)]
pub struct HistogramSeries<'a> {
    pub values: &'a [f64],
    pub label: Option<&'a str>,
    pub bins: usize,
    /// Explicit print-safe palette color index. Series order is used when omitted.
    pub palette_index: Option<usize>,
}

impl<'a> HistogramSeries<'a> {
    pub fn new(values: &'a [f64]) -> Self {
        Self {
            values,
            label: None,
            bins: 40,
            palette_index: None,
        }
    }
    pub fn with_label(mut self, label: &'a str) -> Self {
        self.label = Some(label);
        self
    }
    pub fn with_bins(mut self, bins: usize) -> Self {
        self.bins = bins.max(1);
        self
    }
    pub fn with_palette(mut self, palette_index: usize) -> Self {
        self.palette_index = Some(palette_index);
        self
    }
}

/// Layout and typography configuration for histogram plots.
///
/// [`Self::publication`] is intended for figures that will be reduced to a
/// single-column or thesis text width. The default retains the historical
/// 1920-by-1080 rendering and typography used by [`plot_histograms`].
#[derive(Debug, Clone, Copy)]
pub struct HistogramConfig<'a> {
    pub title: &'a str,
    pub x_label: Option<&'a str>,
    pub y_label: Option<&'a str>,
    pub save_path: Option<&'a str>,
    pub size: (u32, u32),
    pub font_scale: f64,
    /// Fraction of the chart width reserved at the right for the legend.
    pub legend_width_fraction: f64,
    /// Draw the title and legend above the chart instead of over the data.
    pub legend_in_header: bool,
    /// Draw a vertical guide at this x-value.
    pub reference_line: Option<f64>,
    /// Draw histogram interval labels horizontally below the x-axis.
    pub horizontal_bin_labels: bool,
}

impl<'a> HistogramConfig<'a> {
    pub fn new(title: &'a str) -> Self {
        Self {
            title,
            x_label: None,
            y_label: None,
            save_path: None,
            size: (1920, 1080),
            font_scale: 1.0,
            legend_width_fraction: 0.0,
            legend_in_header: false,
            reference_line: None,
            horizontal_bin_labels: false,
        }
    }

    /// Use publication-scale text and reserve space so the legend does not
    /// obscure histogram bars.
    pub fn publication(title: &'a str) -> Self {
        Self::new(title)
            .with_size(1200, 675)
            .with_font_scale(1.4)
            .with_legend_in_header(true)
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

    pub fn with_size(mut self, width: u32, height: u32) -> Self {
        assert!(
            width > 0 && height > 0,
            "histogram dimensions must be nonzero"
        );
        self.size = (width, height);
        self
    }

    pub fn with_font_scale(mut self, scale: f64) -> Self {
        assert!(
            scale.is_finite() && scale > 0.0,
            "histogram font scale must be positive and finite"
        );
        self.font_scale = scale;
        self
    }

    pub fn with_legend_width_fraction(mut self, fraction: f64) -> Self {
        assert!(
            fraction.is_finite() && (0.0..0.5).contains(&fraction),
            "histogram legend width fraction must be in [0, 0.5)"
        );
        self.legend_width_fraction = fraction;
        self
    }

    pub fn with_legend_in_header(mut self, enabled: bool) -> Self {
        self.legend_in_header = enabled;
        self
    }

    pub fn with_reference_line(mut self, value: f64) -> Self {
        assert!(value.is_finite(), "histogram reference line must be finite");
        self.reference_line = Some(value);
        self
    }

    pub fn with_horizontal_bin_labels(mut self, enabled: bool) -> Self {
        self.horizontal_bin_labels = enabled;
        self
    }

    fn assert_valid(&self) {
        assert!(
            self.size.0 > 0 && self.size.1 > 0,
            "histogram dimensions must be nonzero"
        );
        assert!(
            self.font_scale.is_finite() && self.font_scale > 0.0,
            "histogram font scale must be positive and finite"
        );
        assert!(
            self.legend_width_fraction.is_finite()
                && (0.0..0.5).contains(&self.legend_width_fraction),
            "histogram legend width fraction must be in [0, 0.5)"
        );
    }
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
        Self {
            xs,
            ys,
            label: None,
        }
    }

    pub fn with_label(mut self, label: &'a str) -> Self {
        self.label = Some(label);
        self
    }
}

#[cfg(feature = "bitmap")]
use plotters::prelude::BitMapBackend;
#[cfg(feature = "svg")]
use plotters::prelude::SVGBackend;

use plotters::coord::Shift;
use plotters::coord::types::RangedCoordf64;
use plotters::element::DashedPathElement;
use plotters::prelude::*;
use plotters::series::{DashedLineSeries, LineSeries};
use plotters::style::FontTransform;
use plotters::style::text_anchor::{HPos, Pos, VPos};

const FONT_FAMILY: &str = "monospace";
const SMALL_TITLE_FONT_SIZE: u32 = 28;
const TITLE_FONT_SIZE: u32 = 36;
const AXIS_TICK_FONT_SIZE: u32 = 24;
const AXIS_DESC_FONT_SIZE: u32 = 26;
const LEGEND_FONT_SIZE: u32 = 24;
const HISTOGRAM_BIN_LABEL_FONT_SIZE: u32 = 18;
const SMALL_X_LABEL_AREA_SIZE: u32 = 55;
const SMALL_Y_LABEL_AREA_SIZE: u32 = 80;
const X_LABEL_AREA_SIZE: u32 = 90;
const Y_LABEL_AREA_SIZE: u32 = 120;
const HISTOGRAM_X_LABEL_AREA_SIZE: u32 = 190;
const HISTOGRAM_BIN_LABEL_OFFSET: i32 = 12;
const LEGEND_MARKER_WIDTH: i32 = 28;

fn scaled_u32(value: u32, scale: f64) -> u32 {
    ((f64::from(value) * scale).round() as u32).max(1)
}

fn scaled_i32(value: i32, scale: f64) -> i32 {
    ((f64::from(value) * scale).round() as i32).max(1)
}

const PRINT_SAFE_COLORS: [RGBColor; 10] = [
    RGBColor(0, 114, 178),
    RGBColor(213, 94, 0),
    RGBColor(0, 158, 115),
    RGBColor(204, 121, 167),
    RGBColor(86, 180, 233),
    RGBColor(86, 86, 86),
    RGBColor(148, 103, 189),
    RGBColor(214, 39, 40),
    RGBColor(44, 160, 44),
    RGBColor(0, 0, 0),
];

#[derive(Clone, Copy)]
enum HeaderLegendGlyph {
    Fill(usize),
    SolidLine(usize),
    DashedLine(usize),
    Point(usize),
}

#[derive(Clone, Copy)]
struct HeaderLegendItem<'a> {
    label: &'a str,
    glyph: HeaderLegendGlyph,
}

fn draw_publication_header<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    title: &str,
    items: &[HeaderLegendItem<'_>],
    title_font_size: u32,
    legend_font_size: u32,
    marker_width: i32,
    marker_half_height: i32,
) -> Result<DrawingArea<DB, Shift>, DrawingAreaErrorKind<DB::ErrorType>> {
    let (width, _) = root.dim_in_pixel();
    let width = width as i32;
    let outer_margin = (legend_font_size as i32 / 2).max(8);
    let marker_text_gap = (legend_font_size as i32 / 3).max(8);
    // Keep the next glyph visibly separate from the preceding label.  A gap
    // comparable to marker_text_gap makes a line glyph between two labels
    // ambiguous; one em preserves the conventional "glyph, then label"
    // grouping even after the publication figure is reduced in the thesis.
    let item_gap = (legend_font_size as i32).max(18);
    let available_width = (width - outer_margin * 2).max(1);
    let legend_style = TextStyle::from((FONT_FAMILY, legend_font_size).into_font());

    let mut rows: Vec<Vec<(usize, i32)>> = Vec::new();
    let mut row = Vec::new();
    let mut row_width = 0_i32;
    for (index, item) in items.iter().enumerate() {
        let (text_width, _) = root.estimate_text_size(item.label, &legend_style)?;
        let item_width = marker_width + marker_text_gap + text_width as i32 + item_gap;
        if !row.is_empty() && row_width + item_width - item_gap > available_width {
            rows.push(row);
            row = Vec::new();
            row_width = 0;
        }
        row.push((index, item_width));
        row_width += item_width;
    }
    if !row.is_empty() {
        rows.push(row);
    }

    let title_height = title_font_size as i32 + 12;
    let row_height = (legend_font_size as i32 + 12).max(marker_half_height * 2 + 10);
    let header_height = title_height + row_height * rows.len() as i32 + 8;
    let (header, chart_area) = root.split_vertically(header_height as u32);

    let title_style = TextStyle::from((FONT_FAMILY, title_font_size).into_font())
        .pos(Pos::new(HPos::Center, VPos::Center));
    header.draw(&Text::new(
        title,
        (width / 2, title_height / 2),
        title_style,
    ))?;

    for (row_index, row) in rows.iter().enumerate() {
        let total_width = row.iter().map(|(_, item_width)| *item_width).sum::<i32>() - item_gap;
        let mut x = (width - total_width) / 2;
        let y = title_height + row_height * row_index as i32 + row_height / 2;

        for (item_index, item_width) in row {
            let item = items[*item_index];
            let color = match item.glyph {
                HeaderLegendGlyph::Fill(index)
                | HeaderLegendGlyph::SolidLine(index)
                | HeaderLegendGlyph::DashedLine(index)
                | HeaderLegendGlyph::Point(index) => print_safe_color(index),
            };
            match item.glyph {
                HeaderLegendGlyph::Fill(_) => header.draw(&Rectangle::new(
                    [
                        (x, y - marker_half_height),
                        (x + marker_width, y + marker_half_height),
                    ],
                    color.filled(),
                ))?,
                HeaderLegendGlyph::SolidLine(_) => header.draw(&PathElement::new(
                    vec![(x, y), (x + marker_width, y)],
                    color.stroke_width(3),
                ))?,
                HeaderLegendGlyph::DashedLine(_) => header.draw(&DashedPathElement::new(
                    vec![(x, y), (x + marker_width, y)],
                    5,
                    3,
                    color.stroke_width(3),
                ))?,
                HeaderLegendGlyph::Point(_) => header.draw(&Circle::new(
                    (x + marker_width / 2, y),
                    marker_half_height,
                    color.filled(),
                ))?,
            }

            let text_style = legend_style.clone().pos(Pos::new(HPos::Left, VPos::Center));
            header.draw(&Text::new(
                item.label,
                (x + marker_width + marker_text_gap, y),
                text_style,
            ))?;
            x += *item_width;
        }
    }

    Ok(chart_area)
}

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

fn print_safe_color(index: usize) -> RGBAColor {
    PRINT_SAFE_COLORS[index % PRINT_SAFE_COLORS.len()].to_rgba()
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
        let pad = if min_e == 0.0 { 1.0 } else { min_e.abs() * 0.1 };
        min_e -= pad;
        max_e += pad;
    }

    #[cfg(feature = "svg")]
    {
        return with_svg_root(save_path, size, |root| {
            root.fill(&WHITE)?;
            let mut chart = ChartBuilder::on(&root)
                .margin(10)
                .caption(
                    "Error: f_approx(x) - f_true(x)",
                    (FONT_FAMILY, SMALL_TITLE_FONT_SIZE),
                )
                .x_label_area_size(SMALL_X_LABEL_AREA_SIZE)
                .y_label_area_size(SMALL_Y_LABEL_AREA_SIZE)
                .build_cartesian_2d(a..b, min_e..max_e)?;
            let y_label_fmt = |v: &f64| format_axis_value(*v);
            chart
                .configure_mesh()
                .x_desc("x")
                .y_desc("error")
                .y_label_formatter(&y_label_fmt)
                .label_style((FONT_FAMILY, AXIS_TICK_FONT_SIZE))
                .axis_desc_style((FONT_FAMILY, AXIS_DESC_FONT_SIZE))
                .draw()?;
            chart.draw_series(LineSeries::new(pts, &BLACK))?;
            if save_path.is_some() {
                root.present()?;
            }
            Ok(())
        });
    }

    #[cfg(all(not(feature = "svg"), feature = "bitmap"))]
    {
        return with_bitmap_root(save_path, size, |root| {
            root.fill(&WHITE)?;
            let mut chart = ChartBuilder::on(&root)
                .margin(10)
                .caption(
                    "Error: f_approx(x) - f_true(x)",
                    (FONT_FAMILY, SMALL_TITLE_FONT_SIZE),
                )
                .x_label_area_size(SMALL_X_LABEL_AREA_SIZE)
                .y_label_area_size(SMALL_Y_LABEL_AREA_SIZE)
                .build_cartesian_2d(a..b, min_e..max_e)?;
            let y_label_fmt = |v: &f64| format_axis_value(*v);
            chart
                .configure_mesh()
                .x_desc("x")
                .y_desc("error")
                .y_label_formatter(&y_label_fmt)
                .label_style((FONT_FAMILY, AXIS_TICK_FONT_SIZE))
                .axis_desc_style((FONT_FAMILY, AXIS_DESC_FONT_SIZE))
                .draw()?;
            chart.draw_series(LineSeries::new(pts, &BLACK))?;
            if save_path.is_some() {
                root.present()?;
            }
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
    config.assert_valid();
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

            let points =
                s.xs.iter()
                    .zip(s.ys.iter())
                    .filter_map(|(&x, &y)| transform_point(x, y, config.x_scale, config.y_scale))
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

    if x_min == x_max {
        x_max = x_min + 1.0;
    }
    if y_min == y_max {
        y_max = y_min + 1.0;
    }

    let x_pad = (x_max - x_min).abs() * 0.03;
    let y_pad = (y_max - y_min).abs() * 0.03;
    let x_range = (x_min - x_pad)..(x_max + x_pad);
    let y_range = (y_min - y_pad)..(y_max + y_pad);
    let y_zero_threshold = (y_max - y_min).abs() * 1.0e-12;
    let x_desc = axis_label(config.x_label, config.x_scale);
    let y_desc = axis_label(config.y_label, config.y_scale);
    let size = config.size;
    let margin = scaled_u32(20, config.font_scale);
    let title_font_size = scaled_u32(TITLE_FONT_SIZE, config.font_scale);
    let axis_tick_font_size = scaled_u32(AXIS_TICK_FONT_SIZE, config.font_scale);
    let axis_desc_font_size = scaled_u32(AXIS_DESC_FONT_SIZE, config.font_scale);
    let legend_font_size = scaled_u32(LEGEND_FONT_SIZE, config.font_scale);
    let x_label_area_size = scaled_u32(X_LABEL_AREA_SIZE, config.font_scale);
    let y_label_area_size = scaled_u32(Y_LABEL_AREA_SIZE, config.font_scale);
    let legend_marker_width = scaled_i32(LEGEND_MARKER_WIDTH, config.font_scale);
    let header_marker_width = scaled_i32(17, config.font_scale);
    let legend_marker_half_height = scaled_i32(6, config.font_scale);
    let header_items = prepared
        .iter()
        .enumerate()
        .filter_map(|(index, series)| {
            let label = series.label?;
            let palette_index = match series.style {
                PlotStyle::LineWithPalette(index) | PlotStyle::DashedLineWithPalette(index) => {
                    index
                }
                _ => index,
            };
            let glyph = match series.style {
                PlotStyle::Line | PlotStyle::LineWithPalette(_) => {
                    HeaderLegendGlyph::SolidLine(palette_index)
                }
                PlotStyle::DashedLineWithPalette(_) => HeaderLegendGlyph::DashedLine(palette_index),
                PlotStyle::Points => HeaderLegendGlyph::Point(palette_index),
            };
            Some(HeaderLegendItem { label, glyph })
        })
        .collect::<Vec<_>>();
    let use_header_legend = config.legend_in_header && !header_items.is_empty();

    #[cfg(feature = "svg")]
    {
        return with_svg_root(config.save_path, size, |root| {
            root.fill(&WHITE)?;
            let chart_root = if use_header_legend {
                draw_publication_header(
                    &root,
                    config.title,
                    &header_items,
                    title_font_size,
                    legend_font_size,
                    header_marker_width,
                    legend_marker_half_height,
                )?
            } else {
                root.clone()
            };
            let mut chart_builder = ChartBuilder::on(&chart_root);
            chart_builder
                .margin(margin)
                .x_label_area_size(x_label_area_size)
                .y_label_area_size(y_label_area_size);
            if !use_header_legend {
                chart_builder.caption(config.title, (FONT_FAMILY, title_font_size));
            }
            let mut chart = chart_builder.build_cartesian_2d(x_range.clone(), y_range.clone())?;

            let y_label_fmt = |v: &f64| {
                if v.abs() <= y_zero_threshold {
                    "0".to_string()
                } else {
                    format_axis_value(*v)
                }
            };
            let mut mesh = chart.configure_mesh();
            if let Some(ref lbl) = x_desc {
                mesh.x_desc(lbl);
            }
            if let Some(ref lbl) = y_desc {
                mesh.y_desc(lbl);
            }
            mesh.y_label_formatter(&y_label_fmt)
                .label_style((FONT_FAMILY, axis_tick_font_size))
                .axis_desc_style((FONT_FAMILY, axis_desc_font_size))
                .draw()?;

            for (idx, s) in prepared.iter().enumerate() {
                let palette_index = match s.style {
                    PlotStyle::LineWithPalette(index) | PlotStyle::DashedLineWithPalette(index) => {
                        index
                    }
                    _ => idx,
                };
                let color = print_safe_color(palette_index);

                match s.style {
                    PlotStyle::Line | PlotStyle::LineWithPalette(_) => {
                        let ds =
                            chart.draw_series(LineSeries::new(s.points.iter().copied(), &color))?;
                        if let Some(lbl) = s.label {
                            let _ = ds.label(lbl).legend(move |(x, y)| {
                                PathElement::new(
                                    vec![(x, y), (x + legend_marker_width, y)],
                                    color.stroke_width(3),
                                )
                            });
                        }
                    }
                    PlotStyle::DashedLineWithPalette(_) => {
                        let ds = chart.draw_series(DashedLineSeries::new(
                            s.points.iter().copied(),
                            10,
                            6,
                            color.stroke_width(2),
                        ))?;
                        if let Some(lbl) = s.label {
                            let _ = ds.label(lbl).legend(move |(x, y)| {
                                DashedPathElement::new(
                                    vec![(x, y), (x + legend_marker_width, y)],
                                    5,
                                    3,
                                    color.stroke_width(3),
                                )
                            });
                        }
                    }
                    PlotStyle::Points => {
                        let ds = chart.draw_series(
                            s.points
                                .iter()
                                .map(|&(x, y)| Circle::new((x, y), 3, color.filled())),
                        )?;
                        if let Some(lbl) = s.label {
                            let _ = ds.label(lbl).legend(move |(x, y)| {
                                Circle::new(
                                    (x + legend_marker_width / 2, y),
                                    legend_marker_half_height,
                                    color.filled(),
                                )
                            });
                        }
                    }
                }
            }

            if !use_header_legend && prepared.iter().any(|s| s.label.is_some()) {
                chart
                    .configure_series_labels()
                    .border_style(&BLACK)
                    .background_style(WHITE.mix(0.8))
                    .legend_area_size(legend_marker_width + scaled_i32(12, config.font_scale))
                    .label_font((FONT_FAMILY, legend_font_size))
                    .draw()?;
            }

            if config.save_path.is_some() {
                root.present()?;
            }
            Ok(())
        });
    }

    #[cfg(all(not(feature = "svg"), feature = "bitmap"))]
    {
        return with_bitmap_root(config.save_path, size, |root| {
            root.fill(&WHITE)?;
            let chart_root = if use_header_legend {
                draw_publication_header(
                    &root,
                    config.title,
                    &header_items,
                    title_font_size,
                    legend_font_size,
                    header_marker_width,
                    legend_marker_half_height,
                )?
            } else {
                root.clone()
            };
            let mut chart_builder = ChartBuilder::on(&chart_root);
            chart_builder
                .margin(margin)
                .x_label_area_size(x_label_area_size)
                .y_label_area_size(y_label_area_size);
            if !use_header_legend {
                chart_builder.caption(config.title, (FONT_FAMILY, title_font_size));
            }
            let mut chart = chart_builder.build_cartesian_2d(x_range.clone(), y_range.clone())?;

            let y_label_fmt = |v: &f64| {
                if v.abs() <= y_zero_threshold {
                    "0".to_string()
                } else {
                    format_axis_value(*v)
                }
            };
            let mut mesh = chart.configure_mesh();
            if let Some(ref lbl) = x_desc {
                mesh.x_desc(lbl);
            }
            if let Some(ref lbl) = y_desc {
                mesh.y_desc(lbl);
            }
            mesh.y_label_formatter(&y_label_fmt)
                .label_style((FONT_FAMILY, axis_tick_font_size))
                .axis_desc_style((FONT_FAMILY, axis_desc_font_size))
                .draw()?;

            for (idx, s) in prepared.iter().enumerate() {
                let palette_index = match s.style {
                    PlotStyle::LineWithPalette(index) | PlotStyle::DashedLineWithPalette(index) => {
                        index
                    }
                    _ => idx,
                };
                let color = print_safe_color(palette_index);

                match s.style {
                    PlotStyle::Line | PlotStyle::LineWithPalette(_) => {
                        let ds =
                            chart.draw_series(LineSeries::new(s.points.iter().copied(), &color))?;
                        if let Some(lbl) = s.label {
                            let _ = ds.label(lbl).legend(move |(x, y)| {
                                PathElement::new(
                                    vec![(x, y), (x + legend_marker_width, y)],
                                    color.stroke_width(3),
                                )
                            });
                        }
                    }
                    PlotStyle::DashedLineWithPalette(_) => {
                        let ds = chart.draw_series(DashedLineSeries::new(
                            s.points.iter().copied(),
                            10,
                            6,
                            color.stroke_width(2),
                        ))?;
                        if let Some(lbl) = s.label {
                            let _ = ds.label(lbl).legend(move |(x, y)| {
                                DashedPathElement::new(
                                    vec![(x, y), (x + legend_marker_width, y)],
                                    5,
                                    3,
                                    color.stroke_width(3),
                                )
                            });
                        }
                    }
                    PlotStyle::Points => {
                        let ds = chart.draw_series(
                            s.points
                                .iter()
                                .map(|&(x, y)| Circle::new((x, y), 3, color.filled())),
                        )?;
                        if let Some(lbl) = s.label {
                            let _ = ds.label(lbl).legend(move |(x, y)| {
                                Circle::new(
                                    (x + legend_marker_width / 2, y),
                                    legend_marker_half_height,
                                    color.filled(),
                                )
                            });
                        }
                    }
                }
            }

            if !use_header_legend && prepared.iter().any(|s| s.label.is_some()) {
                chart
                    .configure_series_labels()
                    .border_style(&BLACK)
                    .background_style(WHITE.mix(0.8))
                    .legend_area_size(legend_marker_width + scaled_i32(12, config.font_scale))
                    .label_font((FONT_FAMILY, legend_font_size))
                    .draw()?;
            }

            if config.save_path.is_some() {
                root.present()?;
            }
            Ok(())
        });
    }

    #[allow(unreachable_code)]
    Err("no backend enabled; enable `svg` or `bitmap` feature".into())
}

/// Plot one or more histograms using shared x-axis bounds.
///
/// When multiple series are provided, each bin is divided into adjacent bars
/// so every series remains visible without translucent overlays.
pub fn plot_histograms(
    series: &[HistogramSeries<'_>],
    title: &str,
    x_label: Option<&str>,
    y_label: Option<&str>,
    save_path: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut config = HistogramConfig::new(title);
    if let Some(label) = x_label {
        config = config.with_x_label(label);
    }
    if let Some(label) = y_label {
        config = config.with_y_label(label);
    }
    if let Some(path) = save_path {
        config = config.with_save_path(path);
    }
    plot_histograms_with_config(series, &config)
}

/// Plot one or more histograms using a reusable layout and typography
/// configuration.
pub fn plot_histograms_with_config(
    series: &[HistogramSeries<'_>],
    config: &HistogramConfig<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    config.assert_valid();
    assert!(
        !series.is_empty(),
        "plot_histograms: need at least one series"
    );

    struct PreparedHistogram<'a> {
        counts: Vec<f64>,
        label: Option<&'a str>,
        palette_index: usize,
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
    let mut x_max = finite_values
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    if x_min == x_max {
        x_min -= 0.5;
        x_max += 0.5;
    }

    let bin_count = series.iter().map(|s| s.bins.max(1)).max().unwrap_or(1);
    let bin_width = (x_max - x_min) / bin_count as f64;

    let prepared = series
        .iter()
        .enumerate()
        .map(|(series_index, s)| {
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
                palette_index: s.palette_index.unwrap_or(series_index),
            }
        })
        .collect::<Vec<_>>();

    let y_max = prepared
        .iter()
        .flat_map(|s| s.counts.iter().copied())
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let has_legend = prepared.iter().any(|s| s.label.is_some());
    let chart_x_max =
        if has_legend && !config.legend_in_header && config.legend_width_fraction > 0.0 {
            let data_width = x_max - x_min;
            x_max + data_width * config.legend_width_fraction / (1.0 - config.legend_width_fraction)
        } else {
            x_max
        };
    let size = config.size;
    let margin = scaled_u32(20, config.font_scale);
    let title_font_size = scaled_u32(TITLE_FONT_SIZE, config.font_scale);
    let axis_tick_font_size = scaled_u32(AXIS_TICK_FONT_SIZE, config.font_scale);
    let axis_desc_font_size = scaled_u32(AXIS_DESC_FONT_SIZE, config.font_scale);
    let legend_font_size = scaled_u32(LEGEND_FONT_SIZE, config.font_scale);
    let histogram_bin_label_font_size =
        scaled_u32(HISTOGRAM_BIN_LABEL_FONT_SIZE, config.font_scale);
    let histogram_x_label_area_size = scaled_u32(
        if config.horizontal_bin_labels {
            HISTOGRAM_X_LABEL_AREA_SIZE
        } else {
            X_LABEL_AREA_SIZE
        },
        config.font_scale,
    );
    let y_label_area_size = scaled_u32(Y_LABEL_AREA_SIZE, config.font_scale);
    let histogram_bin_label_offset = scaled_i32(HISTOGRAM_BIN_LABEL_OFFSET, config.font_scale);
    let legend_marker_width = scaled_i32(LEGEND_MARKER_WIDTH, config.font_scale);
    let header_marker_width = scaled_i32(17, config.font_scale);
    let legend_marker_half_height = scaled_i32(7, config.font_scale);
    let header_items = prepared
        .iter()
        .filter_map(|histogram| {
            Some(HeaderLegendItem {
                label: histogram.label?,
                glyph: HeaderLegendGlyph::Fill(histogram.palette_index),
            })
        })
        .collect::<Vec<_>>();
    let use_header_legend = config.legend_in_header && !header_items.is_empty();

    #[cfg(feature = "svg")]
    {
        return with_svg_root(config.save_path, size, |root| {
            root.fill(&WHITE)?;
            let chart_root = if use_header_legend {
                draw_publication_header(
                    &root,
                    config.title,
                    &header_items,
                    title_font_size,
                    legend_font_size,
                    header_marker_width,
                    legend_marker_half_height,
                )?
            } else {
                root.clone()
            };
            let mut chart_builder = ChartBuilder::on(&chart_root);
            chart_builder
                .margin(margin)
                .x_label_area_size(histogram_x_label_area_size)
                .y_label_area_size(y_label_area_size);
            if !use_header_legend {
                chart_builder.caption(config.title, (FONT_FAMILY, title_font_size));
            }
            let mut chart =
                chart_builder.build_cartesian_2d(x_min..chart_x_max, 0.0_f64..(y_max * 1.1))?;

            let y_label_fmt = |v: &f64| format_axis_value(*v);
            let mut mesh = chart.configure_mesh();
            if let Some(lbl) = config.x_label {
                mesh.x_desc(lbl);
            }
            if let Some(lbl) = config.y_label {
                mesh.y_desc(lbl);
            }
            mesh.x_labels(8)
                .y_labels(8)
                .disable_x_mesh()
                .disable_y_mesh()
                .y_label_formatter(&y_label_fmt)
                .label_style((FONT_FAMILY, axis_tick_font_size))
                .axis_desc_style((FONT_FAMILY, axis_desc_font_size))
                .draw()?;

            if config.horizontal_bin_labels {
                draw_histogram_bin_labels(
                    &root,
                    &chart,
                    x_min,
                    bin_width,
                    bin_count,
                    histogram_bin_label_font_size,
                    histogram_bin_label_offset,
                    true,
                )?;
            }

            for (idx, histogram) in prepared.iter().enumerate() {
                let color = print_safe_color(histogram.palette_index);
                let fill = color.filled();
                let stroke = color.stroke_width(2);
                let ds = chart.draw_series((0..bin_count).map(|bin_index| {
                    let (left, right) =
                        histogram_bar_bounds(x_min, bin_width, bin_index, idx, prepared.len());
                    Rectangle::new(
                        [(left, 0.0_f64), (right, histogram.counts[bin_index])],
                        fill,
                    )
                }))?;

                if let Some(lbl) = histogram.label {
                    let legend_fill = color.filled();
                    let _ = ds.label(lbl).legend(move |(x, y)| {
                        Rectangle::new(
                            [
                                (x, y - legend_marker_half_height),
                                (x + legend_marker_width, y + legend_marker_half_height),
                            ],
                            legend_fill,
                        )
                    });
                }

                chart.draw_series((0..bin_count).map(|bin_index| {
                    let (left, right) =
                        histogram_bar_bounds(x_min, bin_width, bin_index, idx, prepared.len());
                    Rectangle::new(
                        [(left, 0.0_f64), (right, histogram.counts[bin_index])],
                        stroke,
                    )
                }))?;
            }

            if let Some(value) = config.reference_line {
                chart.draw_series(std::iter::once(DashedPathElement::new(
                    vec![(value, 0.0_f64), (value, y_max * 1.1)],
                    6,
                    4,
                    BLACK.stroke_width(2),
                )))?;
            }

            if has_legend && !use_header_legend {
                chart
                    .configure_series_labels()
                    .border_style(&BLACK)
                    .background_style(WHITE.mix(0.8))
                    .legend_area_size(legend_marker_width + scaled_i32(12, config.font_scale))
                    .label_font((FONT_FAMILY, legend_font_size))
                    .draw()?;
            }

            if config.save_path.is_some() {
                root.present()?;
            }
            Ok(())
        });
    }

    #[cfg(all(not(feature = "svg"), feature = "bitmap"))]
    {
        return with_bitmap_root(config.save_path, size, |root| {
            root.fill(&WHITE)?;
            let chart_root = if use_header_legend {
                draw_publication_header(
                    &root,
                    config.title,
                    &header_items,
                    title_font_size,
                    legend_font_size,
                    header_marker_width,
                    legend_marker_half_height,
                )?
            } else {
                root.clone()
            };
            let mut chart_builder = ChartBuilder::on(&chart_root);
            chart_builder
                .margin(margin)
                .x_label_area_size(histogram_x_label_area_size)
                .y_label_area_size(y_label_area_size);
            if !use_header_legend {
                chart_builder.caption(config.title, (FONT_FAMILY, title_font_size));
            }
            let mut chart =
                chart_builder.build_cartesian_2d(x_min..chart_x_max, 0.0_f64..(y_max * 1.1))?;

            let y_label_fmt = |v: &f64| format_axis_value(*v);
            let mut mesh = chart.configure_mesh();
            if let Some(lbl) = config.x_label {
                mesh.x_desc(lbl);
            }
            if let Some(lbl) = config.y_label {
                mesh.y_desc(lbl);
            }
            mesh.x_labels(8)
                .y_labels(8)
                .disable_x_mesh()
                .disable_y_mesh()
                .y_label_formatter(&y_label_fmt)
                .label_style((FONT_FAMILY, axis_tick_font_size))
                .axis_desc_style((FONT_FAMILY, axis_desc_font_size))
                .draw()?;

            if config.horizontal_bin_labels {
                draw_histogram_bin_labels(
                    &root,
                    &chart,
                    x_min,
                    bin_width,
                    bin_count,
                    histogram_bin_label_font_size,
                    histogram_bin_label_offset,
                    true,
                )?;
            }

            for (idx, histogram) in prepared.iter().enumerate() {
                let color = print_safe_color(histogram.palette_index);
                let fill = color.filled();
                let stroke = color.stroke_width(2);
                let ds = chart.draw_series((0..bin_count).map(|bin_index| {
                    let (left, right) =
                        histogram_bar_bounds(x_min, bin_width, bin_index, idx, prepared.len());
                    Rectangle::new(
                        [(left, 0.0_f64), (right, histogram.counts[bin_index])],
                        fill,
                    )
                }))?;

                if let Some(lbl) = histogram.label {
                    let legend_fill = color.filled();
                    let _ = ds.label(lbl).legend(move |(x, y)| {
                        Rectangle::new(
                            [
                                (x, y - legend_marker_half_height),
                                (x + legend_marker_width, y + legend_marker_half_height),
                            ],
                            legend_fill,
                        )
                    });
                }

                chart.draw_series((0..bin_count).map(|bin_index| {
                    let (left, right) =
                        histogram_bar_bounds(x_min, bin_width, bin_index, idx, prepared.len());
                    Rectangle::new(
                        [(left, 0.0_f64), (right, histogram.counts[bin_index])],
                        stroke,
                    )
                }))?;
            }

            if let Some(value) = config.reference_line {
                chart.draw_series(std::iter::once(DashedPathElement::new(
                    vec![(value, 0.0_f64), (value, y_max * 1.1)],
                    6,
                    4,
                    BLACK.stroke_width(2),
                )))?;
            }

            if has_legend && !use_header_legend {
                chart
                    .configure_series_labels()
                    .border_style(&BLACK)
                    .background_style(WHITE.mix(0.8))
                    .legend_area_size(legend_marker_width + scaled_i32(12, config.font_scale))
                    .label_font((FONT_FAMILY, legend_font_size))
                    .draw()?;
            }

            if config.save_path.is_some() {
                root.present()?;
            }
            Ok(())
        });
    }

    #[allow(unreachable_code)]
    Err("no backend enabled; enable `svg` or `bitmap` feature".into())
}

fn histogram_bar_bounds(
    x_min: f64,
    bin_width: f64,
    bin_index: usize,
    series_index: usize,
    series_count: usize,
) -> (f64, f64) {
    debug_assert!(series_count > 0);
    debug_assert!(series_index < series_count);

    let bar_width = bin_width / series_count as f64;
    let left = x_min + bin_index as f64 * bin_width + series_index as f64 * bar_width;
    (left, left + bar_width)
}

fn histogram_bin_label(x_min: f64, bin_width: f64, bin_index: usize, bin_count: usize) -> String {
    let left = x_min + bin_index as f64 * bin_width;
    let right = left + bin_width;
    let left = format_histogram_edge(left, bin_width);
    let right = format_histogram_edge(right, bin_width);

    if bin_index + 1 == bin_count {
        format!("[{left}, {right}]")
    } else {
        format!("[{left}, {right})")
    }
}

fn format_histogram_edge(value: f64, bin_width: f64) -> String {
    let width = bin_width.abs();
    if !value.is_finite() || !width.is_finite() || width == 0.0 {
        return format_axis_value(value);
    }

    let magnitude = width.log10().floor() as i32;
    let decimals = (1 - magnitude).clamp(0, 17) as usize;
    if value != 0.0 && (value.abs() >= 1.0e6 || value.abs() < 1.0e-4) {
        let value_magnitude = value.abs().log10().floor() as i32;
        let scientific_decimals = (1 + value_magnitude - magnitude).clamp(0, 16) as usize;
        return format!("{value:.scientific_decimals$e}");
    }

    let mut formatted = format!("{value:.decimals$}");
    if decimals > 0 {
        while formatted.ends_with('0') {
            formatted.pop();
        }
        if formatted.ends_with('.') {
            formatted.pop();
        }
    }
    if formatted == "-0" {
        "0".to_string()
    } else {
        formatted
    }
}

fn draw_histogram_bin_labels<DB>(
    root: &DrawingArea<DB, Shift>,
    chart: &ChartContext<'_, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    x_min: f64,
    bin_width: f64,
    bin_count: usize,
    font_size: u32,
    label_offset: i32,
    horizontal_labels: bool,
) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>>
where
    DB: DrawingBackend,
{
    let style = (FONT_FAMILY, font_size).into_text_style(root);
    let style = if horizontal_labels {
        style.pos(Pos::new(HPos::Center, VPos::Top))
    } else {
        style
            .transform(FontTransform::Rotate90)
            .pos(Pos::new(HPos::Right, VPos::Center))
    };
    let (_, axis_y) = chart.backend_coord(&(x_min, 0.0_f64));
    let max_label_chars = (0..bin_count)
        .map(|bin_index| histogram_bin_label(x_min, bin_width, bin_index, bin_count).len())
        .max()
        .unwrap_or(0);
    let label_depth = (max_label_chars as f64 * f64::from(font_size) * 0.6).ceil() as i32;
    let label_y = axis_y + label_offset + label_depth;

    for bin_index in 0..bin_count {
        let left = x_min + bin_index as f64 * bin_width;
        let center = left + bin_width * 0.5;
        let (label_x, _) = chart.backend_coord(&(center, 0.0_f64));
        let label = histogram_bin_label(x_min, bin_width, bin_index, bin_count);
        root.draw(&Text::new(label, (label_x, label_y), style.clone()))?;
    }

    Ok(())
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
    assert!(
        !series.is_empty(),
        "plot_correlation: need at least one series"
    );

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
    let s = XYSeries {
        xs,
        ys,
        label: None,
        style,
    };
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
        if s == "-0" {
            "0".to_string()
        } else {
            s.to_string()
        }
    }
}

fn transform_point(x: f64, y: f64, x_scale: AxisScale, y_scale: AxisScale) -> Option<(f64, f64)> {
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

#[cfg(test)]
mod tests {
    use super::{
        HistogramConfig, HistogramSeries, PlotConfig, PlotStyle, XYSeries, histogram_bar_bounds,
        histogram_bin_label, print_safe_color,
    };

    #[test]
    fn histogram_can_use_an_explicit_palette_color() {
        let values = [1.0, 2.0];
        let series = HistogramSeries::new(&values).with_palette(7);

        assert_eq!(series.palette_index, Some(7));
    }

    #[test]
    fn single_histogram_uses_the_entire_bin() {
        assert_eq!(histogram_bar_bounds(-1.0, 2.0, 3, 0, 1), (5.0, 7.0));
    }

    #[test]
    fn multiple_histograms_tile_each_bin_without_overlap() {
        let bars = (0..3)
            .map(|series_index| histogram_bar_bounds(2.0, 6.0, 1, series_index, 3))
            .collect::<Vec<_>>();

        assert_eq!(bars, vec![(8.0, 10.0), (10.0, 12.0), (12.0, 14.0)]);
    }

    #[test]
    fn histogram_bin_labels_show_half_open_ranges_until_the_last_bin() {
        assert_eq!(histogram_bin_label(0.0, 0.25, 0, 4), "[0, 0.25)");
        assert_eq!(histogram_bin_label(0.0, 0.25, 3, 4), "[0.75, 1]");
    }

    #[test]
    fn histogram_bin_labels_follow_the_bin_resolution() {
        assert_eq!(
            histogram_bin_label(49.178_589, 5.334_091, 0, 10),
            "[49.2, 54.5)"
        );
    }

    #[test]
    fn histogram_bin_labels_keep_narrow_ranges_distinct() {
        assert_eq!(histogram_bin_label(1.0, 1.0e-8, 0, 10), "[1, 1.00000001)");
        assert_eq!(
            histogram_bin_label(1_000_000.0, 1.0, 0, 10),
            "[1.0000000e6, 1.0000010e6)"
        );
    }

    #[test]
    fn publication_histograms_use_readable_layout_defaults() {
        let config = HistogramConfig::publication("title");

        assert_eq!(config.size, (1200, 675));
        assert_eq!(config.font_scale, 1.4);
        assert_eq!(config.legend_width_fraction, 0.0);
        assert!(config.legend_in_header);
    }

    #[test]
    #[should_panic(expected = "histogram legend width fraction must be in [0, 0.5)")]
    fn histogram_config_is_revalidated_before_rendering() {
        let values = [1.0, 2.0];
        let series = [HistogramSeries::new(&values)];
        let mut config = HistogramConfig::new("title");
        config.legend_width_fraction = 1.0;

        let _ = super::plot_histograms_with_config(&series, &config);
    }

    #[cfg(feature = "svg")]
    #[test]
    fn publication_histogram_renders_scaled_typography() {
        let values = [1.0, 1.5, 2.0];
        let series = [HistogramSeries::new(&values)
            .with_label("series")
            .with_bins(8)];
        let path = std::env::temp_dir().join(format!(
            "aether-viz-publication-histogram-{}.svg",
            std::process::id()
        ));
        let path_string = path.to_string_lossy().into_owned();
        let config = HistogramConfig::publication("title")
            .with_x_label("x")
            .with_y_label("Samples")
            .with_reference_line(1.0)
            .with_horizontal_bin_labels(true)
            .with_save_path(&path_string);

        super::plot_histograms_with_config(&series, &config).unwrap();
        let svg = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(path).unwrap();

        assert!(svg.contains("<svg width=\"1200\" height=\"675\""));
        assert!(svg.contains("font-size=\"40.32258064516129\""));
        assert!(svg.contains("font-size=\"27.419354838709676\""));
        assert!(svg.contains("Samples"));
        assert!(svg.contains("rotate(90,") || svg.contains("rotate(270,"));
        assert!(svg.contains("[1, 1.12)"));
        assert!(svg.matches("<text").count() >= 24);
        assert!(svg.matches("<polyline").count() > 20);
        assert!(svg.contains("series"));
    }

    #[cfg(feature = "svg")]
    #[test]
    fn publication_series_renders_scaled_typography_and_header_legend() {
        let xs = [0.0, 1.0, 2.0];
        let ys = [0.0, 1.0, 0.0];
        let series = [XYSeries::new(&xs, &ys)
            .with_label("series")
            .with_style(PlotStyle::Line)];
        let path = std::env::temp_dir().join(format!(
            "aether-viz-publication-series-{}.svg",
            std::process::id()
        ));
        let path_string = path.to_string_lossy().into_owned();
        let config = PlotConfig::publication("title")
            .with_x_label("x")
            .with_y_label("y")
            .with_save_path(&path_string);

        super::plot_series_with_config(&series, &config).unwrap();
        let svg = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(path).unwrap();

        assert!(svg.contains("<svg width=\"1200\" height=\"675\""));
        assert!(svg.contains("font-size=\"40.32258064516129\""));
        assert!(svg.contains("font-size=\"27.419354838709676\""));
        assert!(svg.contains("series"));
    }

    #[test]
    fn print_safe_palette_repeats_without_light_yellow() {
        assert_eq!(print_safe_color(0), print_safe_color(10));

        for index in 0..20 {
            let color = print_safe_color(index);
            assert!(
                !(color.0 > 180 && color.1 > 180 && color.2 < 120),
                "color {index} is too yellow for print: {color:?}"
            );
        }
    }
}
