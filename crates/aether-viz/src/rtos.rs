//! Import, summarize, compare, and plot Pulsar RTOS trace bundles.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::error::Error;
use std::fs;
use std::io::{Error as IoError, ErrorKind};
use std::path::{Path, PathBuf};

use crate::{plot_histograms, plot_series, HistogramSeries, PlotStyle, XYSeries};

pub const RTOS_SCHEMA_VERSION: &str = "pulsar-rtos-v1";
pub const METADATA_FILE: &str = "metadata.csv";
pub const TRACE_FILE: &str = "trace.csv";

fn invalid_data(message: impl Into<String>) -> Box<dyn Error> {
    IoError::new(ErrorKind::InvalidData, message.into()).into()
}

#[derive(Debug, Clone, Default)]
pub struct RtosMetadata {
    values: BTreeMap<String, String>,
    units: BTreeMap<String, String>,
}

impl RtosMetadata {
    pub fn get(&self, key: &str) -> Option<&str> {
        self.values.get(key).map(String::as_str)
    }

    pub fn unit(&self, key: &str) -> Option<&str> {
        self.units
            .get(key)
            .map(String::as_str)
            .filter(|unit| !unit.is_empty())
    }

    pub fn entries(&self) -> impl Iterator<Item = (&str, &str, Option<&str>)> {
        self.values.iter().map(|(key, value)| {
            (
                key.as_str(),
                value.as_str(),
                self.units
                    .get(key)
                    .map(String::as_str)
                    .filter(|unit| !unit.is_empty()),
            )
        })
    }

    pub fn parse_u64(&self, key: &str) -> Result<Option<u64>, Box<dyn Error>> {
        self.get(key)
            .map(|value| {
                value.parse::<u64>().map_err(|error| {
                    invalid_data(format!(
                        "metadata {key}={value:?} is not an integer: {error}"
                    ))
                })
            })
            .transpose()
    }

    pub fn parse_bool(&self, key: &str) -> Result<Option<bool>, Box<dyn Error>> {
        self.get(key)
            .map(|value| match value {
                "true" => Ok(true),
                "false" => Ok(false),
                _ => Err(invalid_data(format!(
                    "metadata {key}={value:?} is not true or false"
                ))),
            })
            .transpose()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RtosEventKind {
    Release,
    Ready,
    Wake,
    RunStart,
    Resume,
    Preempt,
    Block,
    Complete,
    IsrEnter,
    IsrExit,
    Calibration,
    QueueDepth,
    BufferUse,
    Fault,
    Other(String),
}

impl RtosEventKind {
    fn parse(value: &str) -> Self {
        match value {
            "release" => Self::Release,
            "ready" => Self::Ready,
            "wake" => Self::Wake,
            "run_start" => Self::RunStart,
            "resume" => Self::Resume,
            "preempt" => Self::Preempt,
            "block" => Self::Block,
            "complete" => Self::Complete,
            "isr_enter" => Self::IsrEnter,
            "isr_exit" => Self::IsrExit,
            "calibration" => Self::Calibration,
            "queue_depth" => Self::QueueDepth,
            "buffer_use" => Self::BufferUse,
            "fault" => Self::Fault,
            other => Self::Other(other.to_string()),
        }
    }
}

#[derive(Debug, Clone)]
pub struct RtosTraceRecord {
    pub sequence: u64,
    pub timestamp_ticks: u32,
    pub unwrapped_ticks: u64,
    pub timestamp_ns: Option<f64>,
    pub event: RtosEventKind,
    pub core_id: u16,
    pub subject_id: Option<u16>,
    pub correlation_id: Option<u32>,
    pub payload: Option<u32>,
}

#[derive(Debug, Clone)]
pub struct RtosRun {
    pub directory: PathBuf,
    pub run_id: String,
    pub metadata: RtosMetadata,
    pub records: Vec<RtosTraceRecord>,
    pub warnings: Vec<String>,
}

impl RtosRun {
    pub fn is_measured(&self) -> bool {
        self.metadata.get("data_origin") == Some("measured")
    }

    pub fn clock_hz(&self) -> Result<Option<u64>, Box<dyn Error>> {
        self.metadata.parse_u64("clock_hz")
    }

    pub fn capture_duration_ticks(&self) -> Result<Option<u64>, Box<dyn Error>> {
        if self.records.is_empty() {
            return Ok(None);
        }
        let cores = self
            .records
            .iter()
            .map(|record| record.core_id)
            .collect::<BTreeSet<_>>();
        if cores.len() > 1 && self.metadata.parse_bool("clock_domains_synchronized")? != Some(true)
        {
            return Ok(None);
        }
        let min = self
            .records
            .iter()
            .map(|record| record.unwrapped_ticks)
            .min();
        let max = self
            .records
            .iter()
            .map(|record| record.unwrapped_ticks)
            .max();
        Ok(min.zip(max).map(|(min, max)| max.saturating_sub(min)))
    }

    pub fn capture_duration_ns(&self) -> Result<Option<f64>, Box<dyn Error>> {
        Ok(self
            .capture_duration_ticks()?
            .zip(self.clock_hz()?)
            .map(|(ticks, hz)| ticks as f64 * 1.0e9 / hz as f64))
    }

    pub fn wake_latencies(&self) -> Vec<WakeLatencySample> {
        let mut pending = HashMap::<(u16, u32), &RtosTraceRecord>::new();
        let mut samples = Vec::new();

        for record in &self.records {
            let key = record.subject_id.zip(record.correlation_id);
            match (&record.event, key) {
                (RtosEventKind::Wake, Some(key)) => {
                    pending.insert(key, record);
                }
                (RtosEventKind::RunStart | RtosEventKind::Resume, Some(key)) => {
                    let Some(wake) = pending.remove(&key) else {
                        continue;
                    };
                    let latency_ticks =
                        record.timestamp_ticks.wrapping_sub(wake.timestamp_ticks) as u64;
                    let latency_ns = self
                        .metadata
                        .get("clock_hz")
                        .and_then(|value| value.parse::<u64>().ok())
                        .map(|hz| latency_ticks as f64 * 1.0e9 / hz as f64);
                    samples.push(WakeLatencySample {
                        subject_id: key.0,
                        correlation_id: key.1,
                        wake_core_id: wake.core_id,
                        run_core_id: record.core_id,
                        wake_sequence: wake.sequence,
                        run_sequence: record.sequence,
                        latency_ticks,
                        latency_ns,
                    });
                }
                _ => {}
            }
        }
        samples
    }

    pub fn wake_latency_summary(&self) -> Result<Option<WakeLatencySummary>, Box<dyn Error>> {
        let samples = self.wake_latencies();
        if samples.is_empty() {
            return Ok(None);
        }
        let ticks = samples
            .iter()
            .map(|sample| sample.latency_ticks as f64)
            .collect::<Vec<_>>();
        let nanoseconds = samples
            .iter()
            .map(|sample| sample.latency_ns)
            .collect::<Option<Vec<_>>>();
        Ok(Some(WakeLatencySummary {
            ticks: DistributionStats::from_values(&ticks).expect("wake samples are non-empty"),
            nanoseconds: nanoseconds
                .as_deref()
                .and_then(DistributionStats::from_values),
            capture_duration_ticks: self.capture_duration_ticks()?,
            capture_duration_ns: self.capture_duration_ns()?,
        }))
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WakeLatencySample {
    pub subject_id: u16,
    pub correlation_id: u32,
    pub wake_core_id: u16,
    pub run_core_id: u16,
    pub wake_sequence: u64,
    pub run_sequence: u64,
    pub latency_ticks: u64,
    pub latency_ns: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistributionStats {
    pub count: usize,
    pub min: f64,
    pub median: f64,
    pub p99: f64,
    pub max: f64,
}

impl DistributionStats {
    pub fn from_values(values: &[f64]) -> Option<Self> {
        let mut sorted = values
            .iter()
            .copied()
            .filter(|value| value.is_finite())
            .collect::<Vec<_>>();
        if sorted.is_empty() {
            return None;
        }
        sorted.sort_by(f64::total_cmp);
        Some(Self {
            count: sorted.len(),
            min: sorted[0],
            median: nearest_rank(&sorted, 0.50),
            p99: nearest_rank(&sorted, 0.99),
            max: sorted[sorted.len() - 1],
        })
    }
}

fn nearest_rank(sorted: &[f64], percentile: f64) -> f64 {
    let rank = (percentile * sorted.len() as f64).ceil().max(1.0) as usize;
    sorted[rank.min(sorted.len()) - 1]
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WakeLatencySummary {
    pub ticks: DistributionStats,
    pub nanoseconds: Option<DistributionStats>,
    pub capture_duration_ticks: Option<u64>,
    pub capture_duration_ns: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MetricComparison {
    pub candidate: f64,
    pub baseline: f64,
    pub absolute_change: f64,
    pub relative_change_percent: Option<f64>,
}

pub fn compare_metric(candidate: f64, baseline: f64) -> MetricComparison {
    let absolute_change = candidate - baseline;
    MetricComparison {
        candidate,
        baseline,
        absolute_change,
        relative_change_percent: (baseline != 0.0).then_some(absolute_change / baseline * 100.0),
    }
}

pub fn load_rtos_run(directory: impl AsRef<Path>) -> Result<RtosRun, Box<dyn Error>> {
    let directory = directory.as_ref();
    let metadata = load_metadata(&directory.join(METADATA_FILE))?;
    let schema = metadata
        .get("schema_version")
        .ok_or_else(|| invalid_data("metadata is missing schema_version"))?;
    if schema != RTOS_SCHEMA_VERSION {
        return Err(invalid_data(format!(
            "unsupported RTOS schema {schema:?}; expected {RTOS_SCHEMA_VERSION:?}"
        )));
    }
    let run_id = metadata
        .get("run_id")
        .ok_or_else(|| invalid_data("metadata is missing run_id"))?
        .to_string();
    let clock_hz = metadata.parse_u64("clock_hz")?;
    if clock_hz == Some(0) {
        return Err(invalid_data("metadata clock_hz must be greater than zero"));
    }
    let records = load_trace(&directory.join(TRACE_FILE), &run_id, clock_hz)?;
    let warnings = validation_warnings(&metadata, &records)?;

    Ok(RtosRun {
        directory: directory.to_path_buf(),
        run_id,
        metadata,
        records,
        warnings,
    })
}

fn load_metadata(path: &Path) -> Result<RtosMetadata, Box<dyn Error>> {
    let mut reader = csv::Reader::from_path(path)?;
    let headers = reader.headers()?.clone();
    let key_index = column(&headers, "key")?;
    let value_index = column(&headers, "value")?;
    let unit_index = column(&headers, "unit")?;
    let mut metadata = RtosMetadata::default();

    for row in reader.records() {
        let row = row?;
        let key = required_field(&row, key_index, "key")?;
        let value = required_field(&row, value_index, "value")?;
        let unit = row.get(unit_index).unwrap_or_default().trim();
        if metadata
            .values
            .insert(key.to_string(), value.to_string())
            .is_some()
        {
            return Err(invalid_data(format!("duplicate metadata key {key:?}")));
        }
        metadata.units.insert(key.to_string(), unit.to_string());
    }
    Ok(metadata)
}

fn load_trace(
    path: &Path,
    expected_run_id: &str,
    clock_hz: Option<u64>,
) -> Result<Vec<RtosTraceRecord>, Box<dyn Error>> {
    let mut reader = csv::Reader::from_path(path)?;
    let headers = reader.headers()?.clone();
    let schema_index = column(&headers, "schema_version")?;
    let run_index = column(&headers, "run_id")?;
    let sequence_index = column(&headers, "sequence")?;
    let timestamp_index = column(&headers, "timestamp_ticks")?;
    let event_index = column(&headers, "event")?;
    let core_index = column(&headers, "core_id")?;
    let subject_index = column(&headers, "subject_id")?;
    let correlation_index = column(&headers, "correlation_id")?;
    let payload_index = column(&headers, "payload")?;
    let mut wrap_state = HashMap::<u16, (u32, u64)>::new();
    let mut records = Vec::new();
    let mut previous_sequence = None;

    for row in reader.records() {
        let row = row?;
        let schema = required_field(&row, schema_index, "schema_version")?;
        if schema != RTOS_SCHEMA_VERSION {
            return Err(invalid_data(format!(
                "trace row uses unsupported schema {schema:?}"
            )));
        }
        let run_id = required_field(&row, run_index, "run_id")?;
        if run_id != expected_run_id {
            return Err(invalid_data(format!(
                "trace run_id {run_id:?} does not match metadata {expected_run_id:?}"
            )));
        }
        let sequence = parse_required::<u64>(&row, sequence_index, "sequence")?;
        if previous_sequence.is_some_and(|previous| sequence <= previous) {
            return Err(invalid_data("trace sequence must increase strictly"));
        }
        previous_sequence = Some(sequence);
        let timestamp_ticks = parse_required::<u32>(&row, timestamp_index, "timestamp_ticks")?;
        let core_id = parse_required::<u16>(&row, core_index, "core_id")?;
        let state = wrap_state.entry(core_id).or_insert((timestamp_ticks, 0));
        let difference = state.0.abs_diff(timestamp_ticks);
        let record_epoch = if timestamp_ticks < state.0 && difference > (1_u32 << 31) {
            // A large backwards jump is a forward counter wrap.
            state.1 += 1_u64 << 32;
            state.0 = timestamp_ticks;
            state.1
        } else if timestamp_ticks > state.0
            && difference > (1_u32 << 31)
            && state.1 >= (1_u64 << 32)
        {
            // A nested writer reserved first but completed after a wrap. This
            // record belongs to the preceding epoch and must not move state.
            state.1 - (1_u64 << 32)
        } else {
            // Small backwards movements are nested-writer ordering, not wrap.
            if timestamp_ticks >= state.0 {
                state.0 = timestamp_ticks;
            }
            state.1
        };
        let unwrapped_ticks = record_epoch + timestamp_ticks as u64;
        let timestamp_ns = clock_hz.map(|hz| unwrapped_ticks as f64 * 1.0e9 / hz as f64);

        records.push(RtosTraceRecord {
            sequence,
            timestamp_ticks,
            unwrapped_ticks,
            timestamp_ns,
            event: RtosEventKind::parse(required_field(&row, event_index, "event")?),
            core_id,
            subject_id: parse_optional(&row, subject_index, "subject_id")?,
            correlation_id: parse_optional(&row, correlation_index, "correlation_id")?,
            payload: parse_optional(&row, payload_index, "payload")?,
        });
    }
    Ok(records)
}

fn validation_warnings(
    metadata: &RtosMetadata,
    records: &[RtosTraceRecord],
) -> Result<Vec<String>, Box<dyn Error>> {
    let mut warnings = Vec::new();
    if metadata.get("data_origin") != Some("measured") {
        warnings.push(format!(
            "data_origin is {:?}; do not present this run as measured hardware data",
            metadata.get("data_origin").unwrap_or("missing")
        ));
    }
    if metadata.parse_bool("capture_complete")? != Some(true) {
        warnings.push("capture_complete is not true".to_string());
    }
    if let Some(dropped) = metadata.parse_u64("trace_overflow_count")? {
        if dropped > 0 {
            warnings.push(format!("trace overflow dropped {dropped} records"));
        }
    }
    if metadata.get("clock_hz").is_none() {
        warnings.push("clock_hz is missing; values remain in counter ticks".to_string());
    }
    if metadata
        .get("clock_frequency_basis")
        .is_some_and(|basis| basis.contains("not_externally_verified"))
    {
        warnings.push(
            "clock frequency is configured but not externally verified; converted time inherits that uncertainty"
                .to_string(),
        );
    }
    if metadata.get("firmware_commit") == Some("not_recorded") {
        warnings.push("firmware commit was not recorded".to_string());
    }
    match metadata.get("firmware_tree_dirty") {
        Some("true") => warnings.push(
            "firmware was built from a dirty tree; commit alone does not reproduce this run"
                .to_string(),
        ),
        Some("false") | None => {}
        Some("not_recorded") => warnings.push("firmware tree state was not recorded".to_string()),
        Some(_) => {
            let _ = metadata.parse_bool("firmware_tree_dirty")?;
        }
    }
    if let Some(reported) = metadata.parse_u64("trace_records")? {
        if reported != records.len() as u64 {
            warnings.push(format!(
                "metadata reports {reported} trace records but the file contains {}",
                records.len()
            ));
        }
    }
    let cores = records
        .iter()
        .map(|record| record.core_id)
        .collect::<BTreeSet<_>>();
    if cores.len() > 1 && metadata.parse_bool("clock_domains_synchronized")? != Some(true) {
        warnings.push(
            "multiple cores use unsynchronized clocks; cross-core durations are invalid"
                .to_string(),
        );
    }
    if let (Some(requested), Some(completed)) = (
        metadata.parse_u64("requested_samples")?,
        metadata.parse_u64("completed_samples")?,
    ) {
        if requested != completed {
            warnings.push(format!(
                "requested {requested} samples but completed {completed}"
            ));
        }
    }
    if records.is_empty() {
        warnings.push("trace contains no records".to_string());
    }
    Ok(warnings)
}

fn column(headers: &csv::StringRecord, name: &str) -> Result<usize, Box<dyn Error>> {
    headers
        .iter()
        .position(|header| header.trim() == name)
        .ok_or_else(|| invalid_data(format!("CSV is missing required column {name:?}")))
}

fn required_field<'a>(
    row: &'a csv::StringRecord,
    index: usize,
    name: &str,
) -> Result<&'a str, Box<dyn Error>> {
    let value = row.get(index).unwrap_or_default().trim();
    if value.is_empty() {
        Err(invalid_data(format!("CSV field {name:?} is empty")))
    } else {
        Ok(value)
    }
}

fn parse_required<T>(row: &csv::StringRecord, index: usize, name: &str) -> Result<T, Box<dyn Error>>
where
    T: std::str::FromStr,
    T::Err: std::fmt::Display,
{
    let value = required_field(row, index, name)?;
    value
        .parse::<T>()
        .map_err(|error| invalid_data(format!("CSV field {name}={value:?} is invalid: {error}")))
}

fn parse_optional<T>(
    row: &csv::StringRecord,
    index: usize,
    name: &str,
) -> Result<Option<T>, Box<dyn Error>>
where
    T: std::str::FromStr,
    T::Err: std::fmt::Display,
{
    let value = row.get(index).unwrap_or_default().trim();
    if value.is_empty() {
        Ok(None)
    } else {
        value.parse::<T>().map(Some).map_err(|error| {
            invalid_data(format!("CSV field {name}={value:?} is invalid: {error}"))
        })
    }
}

pub fn import_probe_log(
    log_path: impl AsRef<Path>,
    output_directory: impl AsRef<Path>,
) -> Result<(), Box<dyn Error>> {
    let text = fs::read_to_string(log_path)?;
    let metadata = extract_section(
        &text,
        "PULSAR_RUN_METADATA_BEGIN",
        "PULSAR_RUN_METADATA_END",
    )?;
    let trace = extract_section(&text, "PULSAR_TRACE_BEGIN", "PULSAR_TRACE_END")?;
    let output_directory = output_directory.as_ref();
    fs::create_dir_all(output_directory)?;
    fs::write(output_directory.join(METADATA_FILE), metadata)?;
    fs::write(output_directory.join(TRACE_FILE), trace)?;
    Ok(())
}

fn extract_section(text: &str, begin: &str, end: &str) -> Result<String, Box<dyn Error>> {
    let mut active = false;
    let mut complete = false;
    let mut output = String::new();
    for line in text.lines() {
        let line = line.trim();
        if line.contains(begin) {
            if active {
                return Err(invalid_data(format!("duplicate {begin} marker")));
            }
            active = true;
            continue;
        }
        if line.contains(end) {
            if !active {
                return Err(invalid_data(format!("{end} appeared before {begin}")));
            }
            complete = true;
            break;
        }
        if active {
            output.push_str(line);
            output.push('\n');
        }
    }
    if !active || !complete {
        return Err(invalid_data(format!(
            "incomplete log section {begin}..{end}"
        )));
    }
    Ok(output)
}

pub fn plot_wake_latency_series(
    run: &RtosRun,
    save_path: impl AsRef<Path>,
) -> Result<(), Box<dyn Error>> {
    let samples = run.wake_latencies();
    if samples.is_empty() {
        return Err(invalid_data("run has no correlated wake latency samples"));
    }
    let x = (1..=samples.len())
        .map(|value| value as f64)
        .collect::<Vec<_>>();
    let use_ns = samples.iter().all(|sample| sample.latency_ns.is_some());
    let y = samples
        .iter()
        .map(|sample| {
            if use_ns {
                sample.latency_ns.expect("checked above")
            } else {
                sample.latency_ticks as f64
            }
        })
        .collect::<Vec<_>>();
    let path = save_path
        .as_ref()
        .to_str()
        .ok_or_else(|| invalid_data("plot path is not valid UTF-8"))?;
    plot_series(
        &[XYSeries::new(&x, &y)
            .with_label(&run.run_id)
            .with_style(PlotStyle::Line)],
        "Pulsar wake-to-run latency",
        Some("sample"),
        Some(if use_ns {
            "latency [ns]"
        } else {
            "latency [ticks]"
        }),
        Some(path),
    )
}

pub fn plot_wake_latency_histogram(
    runs: &[(&RtosRun, &str)],
    save_path: impl AsRef<Path>,
) -> Result<(), Box<dyn Error>> {
    if runs.is_empty() {
        return Err(invalid_data("at least one run is required"));
    }
    let sample_sets = runs
        .iter()
        .map(|(run, _)| run.wake_latencies())
        .collect::<Vec<_>>();
    let use_ns = sample_sets
        .iter()
        .flatten()
        .all(|sample| sample.latency_ns.is_some());
    let values = sample_sets
        .iter()
        .map(|samples| {
            samples
                .iter()
                .map(|sample| {
                    if use_ns {
                        sample.latency_ns.expect("checked above")
                    } else {
                        sample.latency_ticks as f64
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    if values.iter().all(Vec::is_empty) {
        return Err(invalid_data("runs have no correlated wake latency samples"));
    }
    let histograms = values
        .iter()
        .zip(runs.iter())
        .filter(|(values, _)| !values.is_empty())
        .map(|(values, (_, label))| HistogramSeries::new(values).with_label(label).with_bins(20))
        .collect::<Vec<_>>();
    let path = save_path
        .as_ref()
        .to_str()
        .ok_or_else(|| invalid_data("plot path is not valid UTF-8"))?;
    plot_histograms(
        &histograms,
        "Pulsar wake-to-run latency distribution",
        Some(if use_ns {
            "latency [ns]"
        } else {
            "latency [ticks]"
        }),
        Some("sample count"),
        Some(path),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEST_ID: AtomicU64 = AtomicU64::new(0);

    fn test_directory() -> PathBuf {
        let id = TEST_ID.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!("aether-rtos-{}-{id}", std::process::id()))
    }

    #[test]
    fn imports_wraps_and_summarizes_without_turning_missing_values_into_zero() {
        let directory = test_directory();
        fs::create_dir_all(&directory).unwrap();
        fs::write(
            directory.join(METADATA_FILE),
            concat!(
                "key,value,unit\n",
                "schema_version,pulsar-rtos-v1,\n",
                "run_id,synthetic-wrap-test,\n",
                "data_origin,synthetic,\n",
                "clock_hz,1000000000,Hz\n",
                "clock_domains_synchronized,true,\n",
                "capture_complete,true,\n",
                "trace_overflow_count,0,records\n",
                "requested_samples,1,samples\n",
                "completed_samples,1,samples\n",
            ),
        )
        .unwrap();
        fs::write(
            directory.join(TRACE_FILE),
            concat!(
                "schema_version,run_id,sequence,timestamp_ticks,event,core_id,subject_id,correlation_id,payload\n",
                "pulsar-rtos-v1,synthetic-wrap-test,0,4294967280,wake,0,1,7,\n",
                "pulsar-rtos-v1,synthetic-wrap-test,1,20,run_start,0,1,7,\n",
            ),
        )
        .unwrap();

        let run = load_rtos_run(&directory).unwrap();
        assert!(!run.is_measured());
        assert!(run
            .warnings
            .iter()
            .any(|warning| warning.contains("synthetic")));
        assert_eq!(run.records[0].payload, None);
        assert_eq!(run.records[1].unwrapped_ticks, (1_u64 << 32) + 20);
        let samples = run.wake_latencies();
        assert_eq!(samples.len(), 1);
        assert_eq!(samples[0].latency_ticks, 36);
        assert_eq!(samples[0].latency_ns, Some(36.0));
        let summary = run.wake_latency_summary().unwrap().unwrap();
        assert_eq!(summary.ticks.count, 1);
        assert_eq!(summary.ticks.p99, 36.0);
        assert_eq!(summary.capture_duration_ticks, Some(36));
    }

    #[test]
    fn extracts_versioned_bundle_from_probe_log() {
        let directory = test_directory();
        fs::create_dir_all(&directory).unwrap();
        let log = directory.join("probe.log");
        fs::write(
            &log,
            concat!(
                "runner noise\n",
                "PULSAR_RUN_METADATA_BEGIN\n",
                "key,value,unit\n",
                "schema_version,pulsar-rtos-v1,\n",
                "run_id,import-test,\n",
                "PULSAR_RUN_METADATA_END\n",
                "PULSAR_TRACE_BEGIN\n",
                "schema_version,run_id,sequence,timestamp_ticks,event,core_id,subject_id,correlation_id,payload\n",
                "PULSAR_TRACE_END\n",
            ),
        )
        .unwrap();
        let output = directory.join("run");
        import_probe_log(&log, &output).unwrap();
        assert!(fs::read_to_string(output.join(METADATA_FILE))
            .unwrap()
            .contains("schema_version,pulsar-rtos-v1"));
        assert!(fs::read_to_string(output.join(TRACE_FILE))
            .unwrap()
            .starts_with("schema_version,run_id"));
    }

    #[test]
    fn comparison_handles_zero_baseline_without_false_percentage() {
        let comparison = compare_metric(5.0, 0.0);
        assert_eq!(comparison.absolute_change, 5.0);
        assert_eq!(comparison.relative_change_percent, None);
    }
}
