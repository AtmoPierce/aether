use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use crate::coefficients::*;
use crate::generators::*;
use serde::Deserialize;

/* ================================================================================================
 * Acceptance
 * ================================================================================================ */

#[derive(Clone, Debug, Default, Deserialize)]
pub struct AcceptanceSpec {
    pub max_ulp_f128: Option<u128>,
    pub max_ulp_f64: Option<u64>,
    pub max_ulp_f32: Option<u32>,
    pub abs_tol_f128: Option<f64>,
    pub abs_tol_f64: Option<f64>,
    pub abs_tol_f32: Option<f32>,
}

#[derive(Clone, Debug)]
pub struct AcceptanceResolved {
    pub max_ulp_f128: u128,
    pub max_ulp_f64: u64,
    pub max_ulp_f32: u32,
    pub abs_tol_f128: f128,
    pub abs_tol_f64: f64,
    pub abs_tol_f32: f32,
}

impl AcceptanceResolved {
    pub fn fallback() -> Self {
        Self {
            max_ulp_f128: 5,
            max_ulp_f64: 5,
            max_ulp_f32: 5,
            abs_tol_f128: 1.0e-18_f128,
            abs_tol_f64: 1.0e-16,
            abs_tol_f32: 2.0e-8,
        }
    }
}

fn apply_acceptance(out: &mut AcceptanceResolved, spec: Option<&AcceptanceSpec>) {
    let Some(s) = spec else { return };
    if let Some(v) = s.max_ulp_f128 {
        out.max_ulp_f128 = v;
    }
    if let Some(v) = s.max_ulp_f64 {
        out.max_ulp_f64 = v;
    }
    if let Some(v) = s.max_ulp_f32 {
        out.max_ulp_f32 = v;
    }
    if let Some(v) = s.abs_tol_f128 {
        out.abs_tol_f128 = v as f128;
    }
    if let Some(v) = s.abs_tol_f64 {
        out.abs_tol_f64 = v;
    }
    if let Some(v) = s.abs_tol_f32 {
        out.abs_tol_f32 = v;
    }
}

pub fn resolve_acceptance(
    defaults: &DefaultsConfig,
    func: &FunctionSpec,
    sweep: Option<&SweepSpec>,
    seg: Option<&SegmentSpec>,
) -> AcceptanceResolved {
    let mut out = AcceptanceResolved::fallback();
    apply_acceptance(&mut out, defaults.acceptance.as_ref());
    apply_acceptance(&mut out, func.acceptance.as_ref());
    if let Some(sw) = sweep {
        apply_acceptance(&mut out, sw.acceptance.as_ref());
    }
    if let Some(sg) = seg {
        apply_acceptance(&mut out, sg.acceptance.as_ref());
    }
    out
}

fn resolve_acceptance_from_spec(
    base: &AcceptanceResolved,
    spec: Option<&AcceptanceSpec>,
) -> AcceptanceResolved {
    let mut out = base.clone();
    apply_acceptance(&mut out, spec);
    out
}

fn meets_acceptance_spec(
    base: &AcceptanceResolved,
    spec: Option<&AcceptanceSpec>,
    e: &CandidateErrors,
) -> bool {
    let acc = resolve_acceptance_from_spec(base, spec);
    meets_acceptance(&acc, e)
}

#[derive(Clone, Debug)]
struct CandidateErrors {
    max_abs128: f128,
    max_ulp128: u128,
    max_abs64: f64,
    max_ulp64: u64,
    max_abs32: f32,
    max_ulp32: u32,
}

fn compute_candidate_errors(
    name: &str,
    dom: PolyDomain,
    interval_x: (f64, f64),
    c128: &CoeffSet128,
    c64: &CoeffSet64,
    c32: &CoeffSet32,
    sweep_samples: usize,
    default_samples: usize,
) -> CandidateErrors {
    let samples = if sweep_samples != 0 {
        sweep_samples
    } else {
        default_samples
    };
    let (max_abs128, max_ulp128) = max_err_ulp_f128(name, dom, interval_x, c128, samples);
    let (max_abs64, max_ulp64) = max_err_ulp_f64(name, dom, interval_x, c64, samples);
    let (max_abs32, max_ulp32) = max_err_ulp_f32(name, dom, interval_x, c32, samples);
    CandidateErrors {
        max_abs128,
        max_ulp128,
        max_abs64,
        max_ulp64,
        max_abs32,
        max_ulp32,
    }
}

fn meets_acceptance(acc: &AcceptanceResolved, e: &CandidateErrors) -> bool {
    if e.max_ulp128 > acc.max_ulp_f128 {
        return false;
    }
    if e.max_ulp64 > acc.max_ulp_f64 {
        return false;
    }
    if e.max_ulp32 > acc.max_ulp_f32 {
        return false;
    }
    if e.max_abs128 > acc.abs_tol_f128 {
        return false;
    }
    if e.max_abs64 > acc.abs_tol_f64 {
        return false;
    }
    if e.max_abs32 > acc.abs_tol_f32 {
        return false;
    }
    true
}

pub fn pick_nk_for_interval(
    name: &str,
    func: &FunctionSpec,
    cfg: &FpxRemezConfig,
    dom: PolyDomain,
    interval_x: (f128, f128),
    acc: &AcceptanceResolved,
) -> Nk {
    let sweep = func
        .sweep
        .as_ref()
        .unwrap_or_else(|| panic!("function `{}` has no sweep spec", name));

    if sweep.candidates.is_empty() {
        panic!("function `{}` has no nk candidates", name);
    }

    let mut best: Option<(Nk, CandidateErrors)> = None;

    for cand in &sweep.candidates {
        let nk = Nk {
            n: cand.n,
            k: cand.k,
        };

        if nk.k != nk.n + 1 {
            continue;
        }

        let (c128, c64, c32, _c16) = solve_segment_coeffs(name, func, cfg, dom, interval_x, nk);

        let errs = compute_candidate_errors(
            name,
            dom,
            (interval_x.0 as f64, interval_x.1 as f64),
            &c128,
            &c64,
            &c32,
            sweep.samples,
            cfg.defaults.samples,
        );

        if meets_acceptance_spec(acc, sweep.acceptance.as_ref(), &errs) {
            best = Some((nk, errs));
            break;
        }
    }

    if let Some((nk, _)) = best {
        return nk;
    }

    panic!("no nk candidate for `{}` met acceptance", name);
}

pub fn candidate_meets_acceptance(
    name: &str,
    dom: PolyDomain,
    interval_x: (f64, f64),
    c128: &CoeffSet128,
    c64: &CoeffSet64,
    c32: &CoeffSet32,
    acc: &AcceptanceResolved,
    sweep_samples: usize,
    default_samples: usize,
) -> bool {
    let errs = compute_candidate_errors(
        name,
        dom,
        interval_x,
        c128,
        c64,
        c32,
        sweep_samples,
        default_samples,
    );
    meets_acceptance(acc, &errs)
}
fn max_err_ulp_f128(
    name: &str,
    dom: PolyDomain,
    interval_x: (f64, f64),
    c128: &CoeffSet128,
    samples: usize,
) -> (f128, u128) {
    let f_ref = ref_f128(name);

    let n = samples.max(2);
    let a = interval_x.0 as f128;
    let b = interval_x.1 as f128;
    let step = (b - a) / ((n - 1) as f128);

    let mut max_abs = 0.0_f128;
    let mut max_ulp = 0u128;

    for i in 0..n {
        let x = a + step * (i as f128);
        let y_ref = f_ref(x);
        let y_poly = if name == "log" {
            eval_coeff128_log_atanh(x, c128)
        } else {
            let shape = match dom {
                PolyDomain::X => SegmentShape::X,
                PolyDomain::Z => SegmentShape::Z(trig_kernel(name)),
            };
            eval_coeff128_shape(x, c128, shape)
        };
        let abs = (y_ref - y_poly).abs();
        let ulp = ulp_diff_f128(y_ref, y_poly);
        if abs > max_abs {
            max_abs = abs;
        }
        if ulp > max_ulp {
            max_ulp = ulp;
        }
    }

    (max_abs, max_ulp)
}

fn max_err_ulp_f64(
    name: &str,
    dom: PolyDomain,
    interval_x: (f64, f64),
    c64: &CoeffSet64,
    samples: usize,
) -> (f64, u64) {
    let f_ref = ref_f64(name);

    let n = samples.max(2);
    let a = interval_x.0;
    let b = interval_x.1;
    let step = (b - a) / ((n - 1) as f64);

    let mut max_abs = 0.0f64;
    let mut max_ulp = 0u64;

    for i in 0..n {
        let x = a + step * (i as f64);
        let y_ref = f_ref(x);
        let y_poly = if name == "log" {
            eval_coeff64_log_atanh(x, c64)
        } else {
            let shape = match dom {
                PolyDomain::X => SegmentShape::X,
                PolyDomain::Z => SegmentShape::Z(trig_kernel(name)),
            };
            eval_coeff64_shape(x, c64, shape)
        };
        let abs = (y_ref - y_poly).abs();
        let ulp = ulp_diff_f64(y_ref, y_poly);
        if abs > max_abs {
            max_abs = abs;
        }
        if ulp > max_ulp {
            max_ulp = ulp;
        }
    }

    (max_abs, max_ulp)
}

fn max_err_ulp_f32(
    name: &str,
    dom: PolyDomain,
    interval_x: (f64, f64),
    c32: &CoeffSet32,
    samples: usize,
) -> (f32, u32) {
    let f_ref = ref_f32(name);

    let n = samples.max(2);
    let a = interval_x.0 as f32;
    let b = interval_x.1 as f32;
    let step = (b - a) / ((n - 1) as f32);

    let mut max_abs = 0.0f32;
    let mut max_ulp = 0u32;

    for i in 0..n {
        let x = a + step * (i as f32);
        let y_ref = f_ref(x);
        let y_poly = if name == "log" {
            eval_coeff32_log_atanh(x, c32)
        } else {
            let shape = match dom {
                PolyDomain::X => SegmentShape::X,
                PolyDomain::Z => SegmentShape::Z(trig_kernel(name)),
            };
            eval_coeff32_shape(x, c32, shape)
        };
        let abs = (y_ref - y_poly).abs();
        let ulp = ulp_diff_f32(y_ref, y_poly);
        if abs > max_abs {
            max_abs = abs;
        }
        if ulp > max_ulp {
            max_ulp = ulp;
        }
    }

    (max_abs, max_ulp)
}

fn ref_f64(name: &str) -> fn(f64) -> f64 {
    match name {
        "sin" | "sin_f" | "sinf" => f64::sin,
        "cos" | "cos_f" | "cosf" => f64::cos,
        "tan" | "tan_f" | "tanf" => f64::tan,
        "sinh" | "sinh_f" | "sinhf" => f64::sinh,
        "cosh" | "cosh_f" | "coshf" => f64::cosh,
        "tanh" | "tanh_f" | "tanhf" => f64::tanh,
        "asin" | "asin_f" | "asinf" => f64::asin,
        "acos" | "acos_f" | "acosf" => f64::acos,
        "atan" | "atan_f" | "atanf" => f64::atan,
        "log" | "log_f" | "logf" => f64::ln,
        "log1p" | "log1p_f" | "log1pf" => f64::ln_1p,
        "exp" | "exp_f" | "expf" => f64::exp,
        "expm1" | "expm1_f" | "expm1f" => f64::exp_m1,
        _ => panic!("no f64 reference mapping for {name}"),
    }
}

fn ref_f32(name: &str) -> fn(f32) -> f32 {
    match name {
        "sin" | "sin_f" | "sinf" => f32::sin,
        "cos" | "cos_f" | "cosf" => f32::cos,
        "tan" | "tan_f" | "tanf" => f32::tan,
        "sinh" | "sinh_f" | "sinhf" => f32::sinh,
        "cosh" | "cosh_f" | "coshf" => f32::cosh,
        "tanh" | "tanh_f" | "tanhf" => f32::tanh,
        "asin" | "asin_f" | "asinf" => f32::asin,
        "acos" | "acos_f" | "acosf" => f32::acos,
        "atan" | "atan_f" | "atanf" => f32::atan,
        "log" | "log_f" | "logf" => f32::ln,
        "log1p" | "log1p_f" | "log1pf" => f32::ln_1p,
        "exp" | "exp_f" | "expf" => f32::exp,
        "expm1" | "expm1_f" | "expm1f" => f32::exp_m1,
        _ => panic!("no f32 reference mapping for {name}"),
    }
}

fn ref_f128(name: &str) -> fn(f128) -> f128 {
    match name {
        "sin" | "sin_f" | "sinf" => f128::sin,
        "cos" | "cos_f" | "cosf" => f128::cos,
        "tan" | "tan_f" | "tanf" => f128::tan,
        "sinh" | "sinh_f" | "sinhf" => f128::sinh,
        "cosh" | "cosh_f" | "coshf" => f128::cosh,
        "tanh" | "tanh_f" | "tanhf" => f128::tanh,
        "asin" | "asin_f" | "asinf" => f128::asin,
        "acos" | "acos_f" | "acosf" => f128::acos,
        "atan" | "atan_f" | "atanf" => f128::atan,
        "log" | "log_f" | "logf" => f128::ln,
        "log1p" | "log1p_f" | "log1pf" => f128::ln_1p,
        "exp" | "exp_f" | "expf" => f128::exp,
        "expm1" | "expm1_f" | "expm1f" => f128::exp_m1,
        _ => panic!("no f128 reference mapping for {name}"),
    }
}

/* ================================================================================================
 * Public Config
 * ================================================================================================ */

#[derive(Debug, Clone, Default)]
pub struct FpxRemezConfig {
    pub output: OutputConfig,
    pub defaults: DefaultsConfig,
    pub functions: BTreeMap<String, FunctionSpec>,
}

#[derive(Debug, Clone)]
pub struct OutputConfig {
    pub rs_out_dir: PathBuf,
    pub json_out_dir: PathBuf,
    pub c_out_dir: PathBuf,
    pub gpu_out_dir: PathBuf,
    pub plots_out_dir: PathBuf,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            rs_out_dir: PathBuf::from("target/rs"),
            json_out_dir: PathBuf::from("target/coeffs"),
            c_out_dir: PathBuf::from("target/c"),
            gpu_out_dir: PathBuf::from("target/gpu"),
            plots_out_dir: PathBuf::from("target/plots"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct DefaultsConfig {
    pub samples: usize,
    pub acceptance: Option<AcceptanceSpec>,
}

impl Default for DefaultsConfig {
    fn default() -> Self {
        Self {
            samples: 10_000,
            acceptance: None,
        }
    }
}

/* ================================================================================================
 * Function / Sweep / Segment
 * ================================================================================================ */

#[derive(Debug, Clone)]
pub struct FunctionSpec {
    pub file_prefix: String,
    pub quadrants: Vec<String>,
    pub mirror: bool,
    pub segments: Vec<SegmentSpec>,
    pub sweep: Option<SweepSpec>,
    pub acceptance: Option<AcceptanceSpec>,
}

#[derive(Debug, Clone)]
pub struct SegmentSpec {
    pub seg: Option<String>,
    pub n: usize,
    pub k: usize,
    pub interval: [TomlNumber; 2],
    pub samples: Option<usize>,
    pub acceptance: Option<AcceptanceSpec>,
}

#[derive(Debug, Clone)]
pub struct SweepSpec {
    pub enabled: bool,
    pub interval: [TomlNumber; 2],
    pub min_width: f64,
    pub max_segments: usize,
    pub split_strategy: SplitStrategy,
    pub merge_pass: bool,
    pub samples: usize,
    pub candidates: Vec<NkCandidate>,
    pub w_madds: f64,
    pub w_segments: f64,
    pub w_selectors: f64,
    pub acceptance: Option<AcceptanceSpec>,
}

#[derive(Debug, Clone)]
pub struct NkCandidate {
    pub n: usize,
    pub k: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitStrategy {
    Mid,
    Worst,
}

/* ================================================================================================
 * Numbers / Expressions
 * ================================================================================================ */

#[derive(Debug, Clone)]
pub enum TomlNumber {
    Num(f64),
    Expr(String),
}

impl TomlNumber {
    pub fn as_f64(&self) -> f64 {
        match self {
            TomlNumber::Num(x) => *x,
            TomlNumber::Expr(s) => parse_math_expr(s),
        }
    }
}

fn parse_math_expr(s: &str) -> f64 {
    let t = s.trim().replace(' ', "");

    if let Ok(v) = t.parse::<f64>() {
        return v;
    }

    let (sign, rest) = if let Some(r) = t.strip_prefix('-') {
        (-1.0_f64, r)
    } else if let Some(r) = t.strip_prefix('+') {
        (1.0_f64, r)
    } else {
        (1.0_f64, t.as_str())
    };

    fn parse_unsigned(rest: &str) -> f64 {
        let pi = core::f64::consts::PI;

        if let Ok(v) = rest.parse::<f64>() {
            return v;
        }

        // Handles: (pi/4)^2, (2*pi/3)^2, (0.5*pi)^2
        if let Some(inner) = rest.strip_prefix('(').and_then(|x| x.strip_suffix(")^2")) {
            let v = parse_math_expr(inner);
            return v * v;
        }

        // Handles: pi^2, pi/4^2, 2*pi/3^2 as "square the whole parsed expression".
        if let Some(base) = rest.strip_suffix("^2") {
            let v = parse_math_expr(base);
            return v * v;
        }

        if rest == "pi" {
            return pi;
        }

        if let Some(d) = rest.strip_prefix("pi/") {
            return pi / parse_math_expr(d);
        }

        // Must come before "*pi", otherwise "2*pi/3" gets caught too early.
        if let Some((k, d)) = rest.split_once("*pi/") {
            return parse_math_expr(k) * pi / parse_math_expr(d);
        }

        if let Some((k, _)) = rest.split_once("*pi") {
            return parse_math_expr(k) * pi;
        }

        if let Some(k) = rest.strip_suffix("*pi") {
            return parse_math_expr(k) * pi;
        }

        panic!("unsupported interval expr: {rest}");
    }

    sign * parse_unsigned(rest)
}

/* ================================================================================================
 * TOML DECODING
 * ================================================================================================ */

mod toml_decode {
    use super::*;

    /* -------------------------------- Manifest -------------------------------- */

    #[derive(Debug, Clone, Deserialize, Default)]
    pub struct TomlManifest {
        #[serde(default)]
        pub output: TomlOutput,
        #[serde(default)]
        pub defaults: TomlDefaults,
        #[serde(default)]
        pub inputs: TomlInputs,
        #[serde(default)]
        pub enable: BTreeMap<String, bool>,
        #[serde(default, rename = "override")]
        pub override_: BTreeMap<String, TomlFunctionOverride>,
    }

    #[derive(Debug, Clone, Deserialize, Default)]
    pub struct TomlInputs {
        #[serde(default)]
        pub files: Vec<String>,
    }

    #[derive(Debug, Clone, Deserialize, Default)]
    pub struct TomlFunctionOverride {
        #[serde(default)]
        pub acceptance: Option<AcceptanceSpec>,
        #[serde(default)]
        pub sweep: Option<TomlSweep>,
    }

    /* -------------------------------- Library -------------------------------- */

    #[derive(Debug, Clone, Deserialize, Default)]
    pub struct TomlLibrary {
        #[serde(default)]
        pub functions: BTreeMap<String, TomlFunction>,
    }

    /* -------------------------------- Shared -------------------------------- */

    #[derive(Debug, Clone, Deserialize, Default)]
    pub struct TomlOutput {
        pub rs_out_dir: Option<String>,
        pub json_out_dir: Option<String>,
        pub c_out_dir: Option<String>,
        pub gpu_out_dir: Option<String>,
        pub plots_out_dir: Option<String>,
    }

    #[derive(Debug, Clone, Deserialize, Default)]
    pub struct TomlDefaults {
        pub samples: Option<usize>,
        #[serde(default)]
        pub acceptance: Option<AcceptanceSpec>,
    }

    #[derive(Debug, Clone, Deserialize)]
    pub struct TomlFunction {
        pub file_prefix: String,
        #[serde(default)]
        pub quadrants: Vec<String>,
        #[serde(default)]
        pub mirror: bool,
        #[serde(default)]
        pub segments: Vec<TomlSegment>,
        #[serde(default)]
        pub sweep: Option<TomlSweep>,
        #[serde(default)]
        pub acceptance: Option<AcceptanceSpec>,
    }

    #[derive(Debug, Clone, Deserialize)]
    #[serde(untagged)]
    pub enum TomlNum {
        Num(f64),
        Expr(String),
    }

    #[derive(Debug, Clone, Deserialize)]
    pub struct TomlSegment {
        pub seg: Option<String>,
        pub n: usize,
        pub k: usize,
        pub interval: [TomlNum; 2],
        pub samples: Option<usize>,
        #[serde(default)]
        pub acceptance: Option<AcceptanceSpec>,
    }

    #[derive(Debug, Clone, Deserialize)]
    pub struct TomlNkCandidate {
        pub n: usize,
        pub k: usize,
    }

    #[derive(Debug, Clone, Deserialize)]
    #[serde(rename_all = "lowercase")]
    pub enum TomlSplitStrategy {
        Mid,
        Worst,
    }

    #[derive(Debug, Clone, Deserialize)]
    pub struct TomlSweep {
        pub enabled: bool,
        pub interval: [TomlNum; 2],
        pub min_width: f64,
        pub max_segments: usize,
        pub split_strategy: TomlSplitStrategy,
        pub merge_pass: bool,
        pub samples: usize,
        #[serde(default)]
        pub candidates: Vec<TomlNkCandidate>,
        #[serde(default)]
        pub w_madds: Option<f64>,
        #[serde(default)]
        pub w_segments: Option<f64>,
        #[serde(default)]
        pub w_selectors: Option<f64>,
        #[serde(default)]
        pub acceptance: Option<AcceptanceSpec>,
    }

    /* -------------------------------- Loader -------------------------------- */

    fn resolve_rel(base: &Path, s: &str) -> PathBuf {
        let p = PathBuf::from(s);
        if p.is_relative() {
            base.join(p)
        } else {
            p
        }
    }

    fn to_num(n: &TomlNum) -> TomlNumber {
        match n {
            TomlNum::Num(v) => TomlNumber::Num(*v),
            TomlNum::Expr(s) => TomlNumber::Expr(s.clone()),
        }
    }

    fn load_library(path: &Path) -> TomlLibrary {
        let txt = std::fs::read_to_string(path)
            .unwrap_or_else(|_| panic!("failed to read {}", path.display()));
        toml::from_str(&txt).unwrap_or_else(|_| panic!("failed to parse {}", path.display()))
    }

    pub fn load_from_path(path: &Path) -> Option<FpxRemezConfig> {
        let base = path.parent()?;
        let txt = std::fs::read_to_string(path).ok()?;
        let manifest: TomlManifest =
            toml::from_str(&txt).expect("Could not deserialize aether-fpx-gen manifest.");

        let mut cfg = FpxRemezConfig::default();

        if let Some(s) = manifest.output.rs_out_dir {
            cfg.output.rs_out_dir = resolve_rel(base, &s);
        }
        if let Some(s) = manifest.output.json_out_dir {
            cfg.output.json_out_dir = resolve_rel(base, &s);
        }
        if let Some(s) = manifest.output.c_out_dir {
            cfg.output.c_out_dir = resolve_rel(base, &s);
        }
        if let Some(s) = manifest.output.gpu_out_dir {
            cfg.output.gpu_out_dir = resolve_rel(base, &s);
        }
        if let Some(s) = manifest.output.plots_out_dir {
            cfg.output.plots_out_dir = resolve_rel(base, &s);
        }

        if let Some(s) = manifest.defaults.samples {
            cfg.defaults.samples = s;
        }
        cfg.defaults.acceptance = manifest.defaults.acceptance;

        let mut funcs: BTreeMap<String, TomlFunction> = BTreeMap::new();
        for f in &manifest.inputs.files {
            println!("Function: {}", f);
            let lib = load_library(&resolve_rel(base, f));
            for (name, fun) in lib.functions {
                if funcs.insert(name.clone(), fun).is_some() {
                    panic!("duplicate function {name}");
                }
            }
        }

        funcs.retain(|k, _| manifest.enable.get(k).copied().unwrap_or(false));

        for (name, ov) in &manifest.override_ {
            let f = funcs
                .get_mut(name)
                .expect("override refers to unknown function");
            if let Some(a) = &ov.acceptance {
                f.acceptance = Some(a.clone());
            }
            if let Some(s) = &ov.sweep {
                f.sweep = Some(s.clone());
            }
        }

        cfg.functions = funcs
            .into_iter()
            .map(|(name, f)| {
                let sweep = f.sweep.map(|sw| SweepSpec {
                    enabled: sw.enabled,
                    interval: [to_num(&sw.interval[0]), to_num(&sw.interval[1])],
                    min_width: sw.min_width,
                    max_segments: sw.max_segments,
                    split_strategy: match sw.split_strategy {
                        TomlSplitStrategy::Mid => SplitStrategy::Mid,
                        TomlSplitStrategy::Worst => SplitStrategy::Worst,
                    },
                    merge_pass: sw.merge_pass,
                    samples: sw.samples,
                    candidates: sw
                        .candidates
                        .into_iter()
                        .map(|c| NkCandidate { n: c.n, k: c.k })
                        .collect(),
                    w_madds: sw.w_madds.unwrap_or(1.0),
                    w_segments: sw.w_segments.unwrap_or(3.0),
                    w_selectors: sw.w_selectors.unwrap_or(1.0),
                    acceptance: sw.acceptance,
                });

                let segments = f
                    .segments
                    .into_iter()
                    .map(|s| SegmentSpec {
                        seg: s.seg,
                        n: s.n,
                        k: s.k,
                        interval: [to_num(&s.interval[0]), to_num(&s.interval[1])],
                        samples: s.samples,
                        acceptance: s.acceptance,
                    })
                    .collect();

                (
                    name,
                    FunctionSpec {
                        file_prefix: f.file_prefix,
                        quadrants: f.quadrants,
                        mirror: f.mirror,
                        segments,
                        sweep,
                        acceptance: f.acceptance,
                    },
                )
            })
            .collect();

        Some(cfg)
    }
}

/* ================================================================================================
 * Public Entry Points
 * ================================================================================================ */

fn find_toml_upwards(name: &str) -> Option<PathBuf> {
    let mut dir = PathBuf::from(std::env::var_os("CARGO_MANIFEST_DIR")?);
    loop {
        let cand = dir.join(name);
        if cand.is_file() {
            return Some(cand);
        }
        let cand = dir.join("aether-fpx-gen").join(name);
        if cand.is_file() {
            return Some(cand);
        }
        if !dir.pop() {
            return None;
        }
    }
}

pub fn load_fpx_remez_toml() -> Option<FpxRemezConfig> {
    static FPX_REMEZ_CONFIG: OnceLock<Option<FpxRemezConfig>> = OnceLock::new();

    FPX_REMEZ_CONFIG
        .get_or_init(|| {
            let p = find_toml_upwards("fpx-remez.toml").expect("Could not find fpx-remez.toml");
            toml_decode::load_from_path(&p)
        })
        .clone()
}

pub fn effective_config() -> FpxRemezConfig {
    let mut cfg = FpxRemezConfig::default();
    if let Some(t) = load_fpx_remez_toml() {
        cfg = t;
    }
    cfg
}

pub fn out_dir_for_use(kind: &str, default_rel: &str) -> PathBuf {
    use std::env;

    let env_key = format!("FPX_REMEZ_{}_OUT_DIR", kind.to_ascii_uppercase());
    if let Some(v) = env::var_os(&env_key) {
        return PathBuf::from(v);
    }

    if let Some(cfg) = load_fpx_remez_toml() {
        return out_dir_for(kind, &cfg, default_rel);
    }

    if let Some(v) = env::var_os("FPX_REMEZ_OUT_DIR") {
        return PathBuf::from(v);
    }

    if let Some(out) = env::var_os("OUT_DIR") {
        return PathBuf::from(out).join("fpx-remez").join(kind);
    }

    PathBuf::from(default_rel)
}

pub fn out_dir_for(kind: &str, cfg: &FpxRemezConfig, default_rel: &str) -> PathBuf {
    use std::env;

    let env_key = format!("FPX_REMEZ_{}_OUT_DIR", kind.to_ascii_uppercase());
    if let Some(v) = env::var_os(&env_key) {
        return PathBuf::from(v);
    }

    let from_cfg = match kind {
        "rs" => Some(cfg.output.rs_out_dir.clone()),
        "c" => Some(cfg.output.c_out_dir.clone()),
        "gpu" => Some(cfg.output.gpu_out_dir.clone()),
        "json" => Some(cfg.output.json_out_dir.clone()),
        "plots" => Some(cfg.output.plots_out_dir.clone()),
        _ => None,
    };

    if let Some(p) = from_cfg {
        return p;
    }

    if let Some(v) = env::var_os("FPX_REMEZ_OUT_DIR") {
        return PathBuf::from(v);
    }

    PathBuf::from(default_rel)
}

fn ordered_i64(x: f64) -> i64 {
    let bits = x.to_bits() as i64;
    if bits < 0 {
        i64::MIN - bits
    } else {
        bits
    }
}

fn ordered_i32(x: f32) -> i32 {
    let bits = x.to_bits() as i32;
    if bits < 0 {
        i32::MIN - bits
    } else {
        bits
    }
}

fn ordered_u128(x: f128) -> u128 {
    let bits = x.to_bits();
    if (bits >> 127) != 0 {
        !bits
    } else {
        bits | (1u128 << 127)
    }
}

fn ulp_diff_f64(a: f64, b: f64) -> u64 {
    if a.is_nan() || b.is_nan() {
        return u64::MAX;
    }
    let ia = ordered_i64(a) as i128;
    let ib = ordered_i64(b) as i128;
    (ia - ib).abs() as u64
}

fn ulp_diff_f32(a: f32, b: f32) -> u32 {
    if a.is_nan() || b.is_nan() {
        return u32::MAX;
    }
    let ia = ordered_i32(a) as i64;
    let ib = ordered_i32(b) as i64;
    (ia - ib).abs() as u32
}

fn ulp_diff_f128(a: f128, b: f128) -> u128 {
    if a.is_nan() || b.is_nan() {
        return u128::MAX;
    }
    let ia = ordered_u128(a);
    let ib = ordered_u128(b);
    ia.abs_diff(ib)
}
