use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use aether_viz::{CsvRealtimeConfig, CsvRealtimePlotter};
use eframe::egui;
use egui_plot::{Legend, Line, Plot, PlotPoints, Points};
use glob::glob;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <csv_file_or_pattern> <x_column> [y_columns_csv] [hz=5] [style=line|points]",
            args.first().map(String::as_str).unwrap_or("csv_tail_egui")
        );
        eprintln!(
            "Example: csv_tail_egui \"target/missions/orbit_raise_100_to_200km/runs/run_*/states/earth.crab\" time altitude 5"
        );
        return Ok(());
    }

    let csv_arg = args[1].clone();
    let x_column = args[2].clone();

    let mut y_columns = Vec::new();
    let mut hz = 5.0_f64;
    let mut point_mode = false;

    for token in args.iter().skip(3) {
        if let Ok(parsed_hz) = token.parse::<f64>() {
            hz = parsed_hz;
            continue;
        }

        if token.eq_ignore_ascii_case("points") {
            point_mode = true;
            continue;
        }
        if token.eq_ignore_ascii_case("line") {
            point_mode = false;
            continue;
        }

        if y_columns.is_empty() {
            y_columns = token
                .split(',')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect::<Vec<_>>();
        }
    }

    let polling_enabled = hz > 0.0;
    let clamped_hz = if polling_enabled { hz.max(0.1) } else { 5.0 };

    let has_pattern = csv_arg.contains('*') || csv_arg.contains('?') || csv_arg.contains('[');
    let csv_paths = resolve_csv_path_pattern_all(&csv_arg)?;
    if csv_paths.is_empty() {
        return Err(format!("No files matched: {}", csv_arg).into());
    }

    let mut sources = Vec::new();
    for path in csv_paths {
        let cfg = CsvRealtimeConfig {
            csv_path: path.clone(),
            output_path: "".into(),
            title: format!("{} vs {}", y_columns.join(", "), x_column),
            x_column: x_column.clone(),
            y_columns: y_columns.clone(),
            x_label: Some(x_column.clone()),
            y_label: Some("value".to_string()),
            style: aether_viz::PlotStyle::Line,
            has_header: true,
            poll_hz: clamped_hz,
            max_points: None,
        };

        let mut plotter = CsvRealtimePlotter::new(cfg);
        plotter.initialize()?;

        let name = derive_run_label(&path);
        sources.push(PlotSource {
            path,
            name,
            plotter,
            enabled: true,
        });
    }

    let app = CsvTailApp {
        sources,
        x_label: x_column,
        y_labels: y_columns.clone(),
        selected_x_map: y_columns
            .iter()
            .filter(|c| classify_component(c) == Some(ComponentKind::X))
            .cloned()
            .map(|c| (c, true))
            .collect(),
        selected_y_map: y_columns
            .iter()
            .filter(|c| classify_component(c) == Some(ComponentKind::Y))
            .cloned()
            .map(|c| (c, true))
            .collect(),
        initialized_selection: !y_columns.is_empty(),
        poll_interval: Duration::from_secs_f64(1.0 / clamped_hz),
        source_refresh_interval: Duration::from_secs_f64((1.0 / clamped_hz).max(1.0)),
        last_poll: Instant::now(),
        last_source_refresh: Instant::now(),
        last_error: None,
        source_pattern: if has_pattern { Some(csv_arg) } else { None },
        polling_enabled,
        point_mode,
    };

    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1280.0, 720.0]),
        ..Default::default()
    };

    let run = eframe::run_native(
        "aether_viz :: csv_tail_egui",
        native_options,
        Box::new(move |_cc| Ok(Box::new(app))),
    );

    if let Err(err) = run {
        return Err(format!("Failed to run egui viewer: {}", err).into());
    }

    Ok(())
}

#[derive(Debug)]
struct PlotSource {
    path: PathBuf,
    name: String,
    plotter: CsvRealtimePlotter,
    enabled: bool,
}

struct CsvTailApp {
    sources: Vec<PlotSource>,
    x_label: String,
    y_labels: Vec<String>,
    selected_x_map: BTreeMap<String, bool>,
    selected_y_map: BTreeMap<String, bool>,
    initialized_selection: bool,
    poll_interval: Duration,
    source_refresh_interval: Duration,
    last_poll: Instant,
    last_source_refresh: Instant,
    last_error: Option<String>,
    source_pattern: Option<String>,
    polling_enabled: bool,
    point_mode: bool,
}

impl eframe::App for CsvTailApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if self.polling_enabled && self.last_poll.elapsed() >= self.poll_interval {
            if self.last_source_refresh.elapsed() >= self.source_refresh_interval {
                self.refresh_sources_from_pattern();
                self.last_source_refresh = Instant::now();
            }
            let mut err_msg = None;
            for source in &mut self.sources {
                if let Err(err) = source.plotter.update_once() {
                    err_msg = Some(err.to_string());
                }
            }
            self.last_error = err_msg;
            self.last_poll = Instant::now();
        }

        egui::TopBottomPanel::top("toolbar")
            .resizable(true)
            .default_height(220.0)
            .min_height(140.0)
            .show(ctx, |ui| {
            ui.horizontal(|ui| {
                if self.polling_enabled {
                    ui.label(format!("Poll: {:.2} Hz", 1.0 / self.poll_interval.as_secs_f64()));
                } else {
                    ui.label("Poll: off (static)");
                }
                ui.separator();
                let points: usize = self
                    .sources
                    .iter()
                    .filter(|s| s.enabled)
                    .map(|s| s.plotter.x_values().len())
                    .sum();
                ui.label(format!("Points(total): {}", points));
                ui.separator();
                let enabled_sources = self.sources.iter().filter(|s| s.enabled).count();
                ui.label(format!("Sources: {}/{}", enabled_sources, self.sources.len()));
                if let Some(first) = self.sources.first() {
                    ui.separator();
                    ui.label(format!("x: {}", first.plotter.x_column()));
                }
                if let Some(err) = &self.last_error {
                    ui.separator();
                    ui.colored_label(egui::Color32::RED, format!("Error: {}", err));
                }
            });

            ui.separator();
            ui.collapsing("Files", |ui| {
                ui.horizontal(|ui| {
                    if ui.button("All").clicked() {
                        for source in &mut self.sources {
                            source.enabled = true;
                        }
                    }
                    if ui.button("None").clicked() {
                        for source in &mut self.sources {
                            source.enabled = false;
                        }
                    }
                });

                egui::ScrollArea::vertical().max_height(90.0).show(ui, |ui| {
                    for source in &mut self.sources {
                        let label = format!("{} ({})", source.name, source.path.display());
                        ui.checkbox(&mut source.enabled, label);
                    }
                });
            });

            let columns = self
                .sources
                .first()
                .map(|s| s.plotter.available_columns().to_vec())
                .unwrap_or_default();

            if !columns.is_empty() {
                let x_candidates = columns
                    .iter()
                    .filter(|col| classify_component(col) == Some(ComponentKind::X))
                    .cloned()
                    .collect::<Vec<_>>();
                let y_candidates = columns
                    .iter()
                    .filter(|col| *col != &self.x_label)
                    .cloned()
                    .collect::<Vec<_>>();

                if !self.initialized_selection {
                    if let Some(first_x) = x_candidates.first() {
                        self.selected_x_map.entry(first_x.clone()).or_insert(true);
                    }
                    if let Some(first_y) = y_candidates.first() {
                        self.selected_y_map.entry(first_y.clone()).or_insert(true);
                    }
                    self.apply_selection();
                    self.initialized_selection = true;
                }

                ui.separator();
                let mut changed = false;
                egui::ScrollArea::both()
                    .max_height(140.0)
                    .show(ui, |ui| {
                    ui.horizontal(|ui| {
                        ui.label("Component columns:");
                        ui.columns(2, |cols| {
                            cols[0].heading("X components");
                            for col in &x_candidates {
                                if col == &self.x_label {
                                    continue;
                                }
                                let selected = self.selected_x_map.entry(col.clone()).or_insert(false);
                                if cols[0].checkbox(selected, col).changed() {
                                    changed = true;
                                }
                            }

                            cols[1].heading("Y / scalar columns");
                            for col in &y_candidates {
                                if col == &self.x_label {
                                    continue;
                                }
                                let selected = self.selected_y_map.entry(col.clone()).or_insert(false);
                                if cols[1].checkbox(selected, col).changed() {
                                    changed = true;
                                }
                            }
                        });
                    });
                });
                if changed {
                    self.apply_selection();
                }
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            if self.sources.is_empty() {
                ui.heading("No data sources loaded");
                return;
            }

            let (selected_x_cols, selected_y_cols) = selected_component_columns(
                &self.y_labels,
                &self.selected_x_map,
                &self.selected_y_map,
            );

            let x_unit = common_unit(&selected_x_cols);
            let y_unit = common_unit(&selected_y_cols);
            let use_dual_axis = !selected_x_cols.is_empty()
                && !selected_y_cols.is_empty()
                && x_unit != y_unit;

            let (x_group_min, x_group_max) =
                min_max_for_columns_multi(&self.sources, &self.y_labels, &selected_x_cols);
            let (y_group_min, y_group_max) =
                min_max_for_columns_multi(&self.sources, &self.y_labels, &selected_y_cols);

            Plot::new("csv_realtime_plot")
                .legend(Legend::default())
                .x_axis_label(self.x_label.clone())
                .y_axis_label("value")
                .show(ui, |plot_ui| {
                    for source in &self.sources {
                        if !source.enabled {
                            continue;
                        }
                        let xs = source.plotter.x_values();
                        let ys_all = source.plotter.y_values();
                        if xs.is_empty() || ys_all.is_empty() {
                            continue;
                        }

                        for (idx, ys) in ys_all.iter().enumerate() {
                            if ys.len() != xs.len() {
                                continue;
                            }

                            let col_name = self
                                .y_labels
                                .get(idx)
                                .cloned()
                                .unwrap_or_else(|| format!("series_{}", idx));

                            let is_y_component = self
                                .selected_y_map
                                .get(&col_name)
                                .copied()
                                .unwrap_or(false);

                            let mut points = xs
                                .iter()
                                .zip(ys.iter())
                                .filter(|(x, y)| x.is_finite() && y.is_finite())
                                .map(|(x, y)| [*x, *y])
                                .collect::<Vec<_>>();

                            if use_dual_axis && is_y_component {
                                for p in &mut points {
                                    p[1] = map_range(
                                        p[1],
                                        y_group_min,
                                        y_group_max,
                                        x_group_min,
                                        x_group_max,
                                    );
                                }
                            }

                            let label = format!("{}::{}", source.name, col_name);
                            if self.point_mode {
                                plot_ui.points(
                                    Points::new(PlotPoints::from(points))
                                        .name(label)
                                        .radius(2.0),
                                );
                            } else {
                                plot_ui.line(Line::new(PlotPoints::from(points)).name(label));
                            }
                        }
                    }
                });

            if use_dual_axis {
                ui.separator();
                ui.label(format!(
                    "Dual-axis mode: X unit '{}' vs Y unit '{}'. Y-series are scaled onto X-axis for plotting.",
                    x_unit.clone().unwrap_or_else(|| "(none)".to_string()),
                    y_unit.clone().unwrap_or_else(|| "(none)".to_string())
                ));
                ui.label(format!(
                    "X range [{:.6}, {:.6}] | Y range [{:.6}, {:.6}]",
                    x_group_min, x_group_max, y_group_min, y_group_max
                ));
            }
        });

        ctx.request_repaint_after(Duration::from_millis(16));
    }
}

impl CsvTailApp {
    fn create_source(
        &self,
        path: PathBuf,
        x_column: &str,
        y_columns: &[String],
    ) -> Result<PlotSource, Box<dyn Error>> {
        let hz = if self.polling_enabled {
            (1.0 / self.poll_interval.as_secs_f64()).max(0.1)
        } else {
            5.0
        };
        let cfg = CsvRealtimeConfig {
            csv_path: path.clone(),
            output_path: "".into(),
            title: format!("{} vs {}", y_columns.join(", "), x_column),
            x_column: x_column.to_string(),
            y_columns: y_columns.to_vec(),
            x_label: Some(x_column.to_string()),
            y_label: Some("value".to_string()),
            style: aether_viz::PlotStyle::Line,
            has_header: true,
            poll_hz: hz,
            max_points: None,
        };

        let mut plotter = CsvRealtimePlotter::new(cfg);
        plotter.initialize()?;

        Ok(PlotSource {
            path: path.clone(),
            name: derive_run_label(&path),
            plotter,
            enabled: true,
        })
    }

    fn refresh_sources_from_pattern(&mut self) {
        let Some(pattern) = self.source_pattern.as_ref() else {
            return;
        };

        let Ok(paths) = resolve_csv_path_pattern_all(pattern) else {
            return;
        };

        for path in paths {
            if self.sources.iter().any(|source| source.path == path) {
                continue;
            }

            match self.create_source(path.clone(), &self.x_label, &self.y_labels) {
                Ok(mut source) => {
                    if !self.y_labels.is_empty() {
                        let _ = source.plotter.set_y_columns(self.y_labels.clone());
                    }
                    self.sources.push(source);
                }
                Err(err) => {
                    self.last_error = Some(format!(
                        "Failed to load new source {}: {}",
                        path.display(),
                        err
                    ));
                }
            }
        }
    }

    fn apply_selection(&mut self) {
        let mut set = BTreeSet::new();
        for (name, on) in &self.selected_x_map {
            if *on {
                set.insert(name.clone());
            }
        }
        for (name, on) in &self.selected_y_map {
            if *on {
                set.insert(name.clone());
            }
        }

        let new_cols = set.into_iter().collect::<Vec<_>>();
        let mut err_msg = None;
        for source in &mut self.sources {
            if let Err(err) = source.plotter.set_y_columns(new_cols.clone()) {
                err_msg = Some(err.to_string());
            }
        }

        if err_msg.is_none() {
            self.y_labels = new_cols;
        }
        self.last_error = err_msg;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ComponentKind {
    X,
    Y,
}

fn classify_component(name: &str) -> Option<ComponentKind> {
    let lower = name.to_ascii_lowercase();
    if lower.contains("_x[") || lower.ends_with("_x") || lower.contains("_x_") {
        return Some(ComponentKind::X);
    }
    if lower.contains("_y[") || lower.ends_with("_y") || lower.contains("_y_") {
        return Some(ComponentKind::Y);
    }
    None
}

fn extract_unit(name: &str) -> Option<String> {
    let start = name.rfind('[')?;
    let end = name.rfind(']')?;
    if end <= start + 1 {
        return None;
    }
    Some(name[start + 1..end].trim().to_string())
}

fn common_unit(columns: &[String]) -> Option<String> {
    let mut unit: Option<String> = None;
    for col in columns {
        let u = extract_unit(col);
        match (&unit, u) {
            (None, Some(v)) => unit = Some(v),
            (None, None) => unit = Some(String::new()),
            (Some(existing), Some(v)) if existing == &v => {}
            (Some(existing), None) if existing.is_empty() => {}
            _ => return None,
        }
    }
    unit.filter(|u| !u.is_empty())
}

fn selected_component_columns(
    labels: &[String],
    selected_x_map: &BTreeMap<String, bool>,
    selected_y_map: &BTreeMap<String, bool>,
) -> (Vec<String>, Vec<String>) {
    let x_cols = labels
        .iter()
        .filter(|name| selected_x_map.get(*name).copied().unwrap_or(false))
        .cloned()
        .collect::<Vec<_>>();
    let y_cols = labels
        .iter()
        .filter(|name| selected_y_map.get(*name).copied().unwrap_or(false))
        .cloned()
        .collect::<Vec<_>>();
    (x_cols, y_cols)
}

fn min_max_for_columns_multi(
    sources: &[PlotSource],
    labels: &[String],
    selected_cols: &[String],
) -> (f64, f64) {
    let mut min_v = f64::INFINITY;
    let mut max_v = f64::NEG_INFINITY;

    for source in sources {
        let xs = source.plotter.x_values();
        let ys_all = source.plotter.y_values();
        for (idx, ys) in ys_all.iter().enumerate() {
            let Some(label) = labels.get(idx) else {
                continue;
            };
            if !selected_cols.contains(label) || ys.len() != xs.len() {
                continue;
            }
            for &y in ys {
                if y.is_finite() {
                    min_v = min_v.min(y);
                    max_v = max_v.max(y);
                }
            }
        }
    }

    if !min_v.is_finite() || !max_v.is_finite() {
        (0.0, 1.0)
    } else if (max_v - min_v).abs() < f64::EPSILON {
        (min_v, min_v + 1.0)
    } else {
        (min_v, max_v)
    }
}

fn map_range(value: f64, src_min: f64, src_max: f64, dst_min: f64, dst_max: f64) -> f64 {
    if (src_max - src_min).abs() < f64::EPSILON {
        return dst_min;
    }
    let t = (value - src_min) / (src_max - src_min);
    dst_min + t * (dst_max - dst_min)
}

fn resolve_csv_path_pattern_all(input: &str) -> Result<Vec<PathBuf>, Box<dyn Error>> {
    let has_pattern = input.contains('*') || input.contains('?') || input.contains('[');
    if !has_pattern {
        return Ok(vec![PathBuf::from(input)]);
    }

    let mut matches: Vec<PathBuf> = glob(input)?
        .filter_map(Result::ok)
        .filter(|p| p.is_file())
        .collect();

    if matches.is_empty() {
        return Err(format!("No files matched pattern: {}", input).into());
    }

    matches.sort();
    Ok(matches)
}

fn derive_run_label(path: &Path) -> String {
    let parts = path
        .components()
        .map(|comp| comp.as_os_str().to_string_lossy().into_owned())
        .collect::<Vec<_>>();

    for s in parts.iter().rev() {
        if s.starts_with("run_") {
            return s.clone();
        }
    }

    path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("source")
        .to_string()
}
