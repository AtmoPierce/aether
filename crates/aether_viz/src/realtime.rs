use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::thread;
use std::time::{Duration, SystemTime};

use crate::{plot_series, PlotStyle, XYSeries};

#[derive(Debug, Clone)]
pub struct CsvRealtimeConfig {
    pub csv_path: PathBuf,
    pub output_path: PathBuf,
    pub title: String,
    pub x_column: String,
    pub y_columns: Vec<String>,
    pub x_label: Option<String>,
    pub y_label: Option<String>,
    pub style: PlotStyle,
    pub has_header: bool,
    pub poll_hz: f64,
    pub max_points: Option<usize>,
}

impl Default for CsvRealtimeConfig {
    fn default() -> Self {
        Self {
            csv_path: PathBuf::from("flight.crab"),
            output_path: PathBuf::from("plot.svg"),
            title: "Realtime CSV Plot".to_string(),
            x_column: "time".to_string(),
            y_columns: vec!["value".to_string()],
            x_label: Some("x".to_string()),
            y_label: Some("y".to_string()),
            style: PlotStyle::Line,
            has_header: true,
            poll_hz: 5.0,
            max_points: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct FileIdentity {
    created: Option<SystemTime>,
    len: u64,
    #[cfg(unix)]
    dev: u64,
    #[cfg(unix)]
    ino: u64,
}

impl FileIdentity {
    fn from_path(path: &Path) -> std::io::Result<Self> {
        let md = std::fs::metadata(path)?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::MetadataExt;
            Ok(Self {
                created: md.created().ok(),
                len: md.len(),
                dev: md.dev(),
                ino: md.ino(),
            })
        }
        #[cfg(not(unix))]
        {
            Ok(Self {
                created: md.created().ok(),
                len: md.len(),
            })
        }
    }
}

#[derive(Debug)]
pub struct CsvRealtimePlotter {
    cfg: CsvRealtimeConfig,
    headers: Vec<String>,
    x_index: Option<usize>,
    y_indices: Vec<usize>,
    x_values: Vec<f64>,
    y_values: Vec<Vec<f64>>,
    partial_line: String,
    file_pos: u64,
    last_identity: Option<FileIdentity>,
}

impl CsvRealtimePlotter {
    pub fn new(cfg: CsvRealtimeConfig) -> Self {
        let y_count = cfg.y_columns.len();
        Self {
            cfg,
            headers: Vec::new(),
            x_index: None,
            y_indices: Vec::new(),
            x_values: Vec::new(),
            y_values: vec![Vec::new(); y_count],
            partial_line: String::new(),
            file_pos: 0,
            last_identity: None,
        }
    }

    pub fn initialize(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        self.reload_from_start()
    }

    pub fn run_forever(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        self.initialize()?;
        self.render()?;

        let hz = if self.cfg.poll_hz > 0.0 { self.cfg.poll_hz } else { 5.0 };
        let dt = Duration::from_secs_f64(1.0 / hz);

        loop {
            let changed = self.update_once()?;
            if changed {
                self.render()?;
            }
            thread::sleep(dt);
        }
    }

    pub fn update_once(&mut self) -> Result<bool, Box<dyn std::error::Error>> {
        if !self.cfg.csv_path.exists() {
            return Ok(false);
        }

        let current_id = FileIdentity::from_path(&self.cfg.csv_path)?;

        if self.should_reset_for_new_file(current_id) {
            self.reset_series();
            self.reload_from_start()?;
            return Ok(true);
        }

        let mut changed = false;
        if current_id.len > self.file_pos {
            changed |= self.read_new_bytes(self.file_pos)?;
        }

        self.last_identity = Some(current_id);
        Ok(changed)
    }

    pub fn render(&self) -> Result<(), Box<dyn std::error::Error>> {
        if self.x_values.is_empty() {
            return Ok(());
        }

        let output_path = self
            .cfg
            .output_path
            .to_str()
            .ok_or("output path is not valid UTF-8")?;

        let mut series = Vec::with_capacity(self.y_values.len());
        for (idx, ys) in self.y_values.iter().enumerate() {
            if ys.len() != self.x_values.len() {
                continue;
            }
            let label = self.cfg.y_columns.get(idx).map(String::as_str);
            let mut s = XYSeries::new(&self.x_values, ys).with_style(self.cfg.style);
            if let Some(lbl) = label {
                s = s.with_label(lbl);
            }
            series.push(s);
        }

        if series.is_empty() {
            return Ok(());
        }

        plot_series(
            &series,
            &self.cfg.title,
            self.cfg.x_label.as_deref(),
            self.cfg.y_label.as_deref(),
            Some(output_path),
        )
    }

    fn should_reset_for_new_file(&self, current: FileIdentity) -> bool {
        let Some(previous) = self.last_identity else {
            return false;
        };

        if current.len < self.file_pos {
            return true;
        }

        #[cfg(unix)]
        {
            if previous.dev != current.dev || previous.ino != current.ino {
                return true;
            }
        }

        if previous.created != current.created {
            return true;
        }

        false
    }

    fn reload_from_start(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        if !self.cfg.csv_path.exists() {
            return Ok(());
        }

        let mut file = File::open(&self.cfg.csv_path)?;
        let mut text = String::new();
        file.read_to_string(&mut text)?;
        self.file_pos = text.as_bytes().len() as u64;
        self.last_identity = Some(FileIdentity::from_path(&self.cfg.csv_path)?);

        let _ = self.parse_text_chunk(&text)?;
        Ok(())
    }

    fn read_new_bytes(&mut self, start: u64) -> Result<bool, Box<dyn std::error::Error>> {
        let mut file = File::open(&self.cfg.csv_path)?;
        file.seek(SeekFrom::Start(start))?;
        let mut text = String::new();
        file.read_to_string(&mut text)?;
        self.file_pos = file.stream_position()?;
        if text.is_empty() {
            return Ok(false);
        }
        self.parse_text_chunk(&text)
    }

    fn parse_text_chunk(&mut self, text: &str) -> Result<bool, Box<dyn std::error::Error>> {
        let mut changed = false;
        let mut buf = String::new();
        if !self.partial_line.is_empty() {
            buf.push_str(&self.partial_line);
            self.partial_line.clear();
        }
        buf.push_str(text);

        let has_trailing_newline = buf.ends_with('\n');
        let mut lines: Vec<&str> = buf.lines().collect();
        if !has_trailing_newline {
            if let Some(last) = lines.pop() {
                self.partial_line.push_str(last);
            }
        }

        for raw in lines {
            let line = raw.trim();
            if line.is_empty() {
                continue;
            }

            if self.cfg.has_header && self.headers.is_empty() {
                self.headers = split_csv_line(line)
                    .into_iter()
                    .map(|s| s.to_string())
                    .collect();
                self.resolve_column_indices()?;
                continue;
            }

            if self.x_index.is_none() || self.y_indices.is_empty() {
                self.resolve_column_indices()?;
            }

            let values = split_csv_line(line);
            if let Some(xi) = self.x_index {
                if xi >= values.len() {
                    continue;
                }
                let x = match values[xi].parse::<f64>() {
                    Ok(v) if v.is_finite() => v,
                    _ => continue,
                };

                let mut y_row = Vec::with_capacity(self.y_indices.len());
                let mut all_ok = true;
                for &yi in &self.y_indices {
                    if yi >= values.len() {
                        all_ok = false;
                        break;
                    }
                    match values[yi].parse::<f64>() {
                        Ok(v) if v.is_finite() => y_row.push(v),
                        _ => {
                            all_ok = false;
                            break;
                        }
                    }
                }
                if !all_ok {
                    continue;
                }

                self.x_values.push(x);
                for (idx, val) in y_row.into_iter().enumerate() {
                    if let Some(col) = self.y_values.get_mut(idx) {
                        col.push(val);
                    }
                }

                self.trim_if_needed();
                changed = true;
            }
        }

        Ok(changed)
    }

    fn resolve_column_indices(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        if self.cfg.has_header {
            if self.headers.is_empty() {
                return Ok(());
            }

            self.x_index = Some(
                self
                .headers
                .iter()
                .position(|h| h == &self.cfg.x_column)
                .ok_or_else(|| format!("x column '{}' not found", self.cfg.x_column))?,
            );

            self.y_indices.clear();
            for y in &self.cfg.y_columns {
                let yi = self
                    .headers
                    .iter()
                    .position(|h| h == y)
                    .ok_or_else(|| format!("y column '{}' not found", y))?;
                self.y_indices.push(yi);
            }
            return Ok(());
        }

        self.x_index = Some(self.cfg.x_column.parse::<usize>()?);
        self.y_indices = self
            .cfg
            .y_columns
            .iter()
            .map(|s| s.parse::<usize>())
            .collect::<Result<Vec<_>, _>>()?;
        Ok(())
    }

    fn trim_if_needed(&mut self) {
        let Some(max_points) = self.cfg.max_points else {
            return;
        };

        if self.x_values.len() <= max_points {
            return;
        }

        let drop_n = self.x_values.len() - max_points;
        self.x_values.drain(0..drop_n);
        for ys in &mut self.y_values {
            if ys.len() >= drop_n {
                ys.drain(0..drop_n);
            } else {
                ys.clear();
            }
        }
    }

    fn reset_series(&mut self) {
        self.headers.clear();
        self.x_index = None;
        self.y_indices.clear();
        self.x_values.clear();
        for ys in &mut self.y_values {
            ys.clear();
        }
        self.partial_line.clear();
        self.file_pos = 0;
        self.last_identity = None;
    }

    pub fn x_values(&self) -> &[f64] {
        &self.x_values
    }

    pub fn y_values(&self) -> &[Vec<f64>] {
        &self.y_values
    }

    pub fn available_columns(&self) -> &[String] {
        &self.headers
    }

    pub fn x_column(&self) -> &str {
        &self.cfg.x_column
    }

    pub fn selected_y_columns(&self) -> &[String] {
        &self.cfg.y_columns
    }

    pub fn set_y_columns(
        &mut self,
        y_columns: Vec<String>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        self.cfg.y_columns = y_columns;
        self.y_values = vec![Vec::new(); self.cfg.y_columns.len()];
        self.reset_series();
        self.reload_from_start()
    }
}

fn split_csv_line(line: &str) -> Vec<&str> {
    line.split(',').map(|s| s.trim()).collect()
}
