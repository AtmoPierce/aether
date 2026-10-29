// aether-fpx-bench/src/types.rs
#[derive(Debug, Clone)]
pub struct BenchConfig {
    pub bin_name: String,
    pub ops: Vec<String>,
    pub modes: Vec<String>,
    pub precisions: Vec<String>,
    pub events: Vec<String>,
    pub runs: usize,
    pub iters: usize,
}

#[derive(Debug, Clone)]
pub struct Row {
    pub op: String,
    pub mode: String,
    pub precision: String,
    pub run: usize,
    /// Function calls made by this sample; retained so reports cannot be
    /// accidentally normalized with a different collection-time value.
    pub iterations: usize,

    pub instructions: Option<f64>,
    pub cycles: Option<f64>,
    pub branches: Option<f64>,
    pub branch_misses: Option<f64>,

    pub time_elapsed_s: Option<f64>,
    /// Energy attributed to the benchmark by an optional perf energy PMU.
    pub energy_joules: Option<f64>,
    /// Mean power over the measured wall-clock interval.
    pub average_power_watts: Option<f64>,
    /// Machine-readable description of the energy/power measurement source.
    pub power_source: Option<String>,
}

#[derive(Debug, Clone)]
pub struct Summary {
    pub sum_ins: f64,
    pub sum_cyc: f64,
    pub sum_br: f64,
    pub sum_bm: f64,
    pub sum_time_s: f64,
    pub n: usize,

    pub cpi_agg: Option<f64>,
    pub ipc_agg: Option<f64>,
    pub ns_per_call: Option<f64>,
    pub ins_per_call: Option<f64>,
    pub cyc_per_call: Option<f64>,
    pub bm_rate: Option<f64>,
}
