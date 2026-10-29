#![feature(f128)]
#![cfg(feature = "ocp")]

use std::env;
use std::fs;
use std::path::PathBuf;

use aether_viz::{plot_series_with_config, PlotConfig, PlotStyle, XYSeries};
use serde::Deserialize;

use aether_fpx_bench::{valid_ulp_dist_f128, valid_ulp_dist_f32, valid_ulp_dist_f64};

#[cfg(feature = "mx")]
use aether_fpx::ocp::fp4::Fp4E2M1;
#[cfg(feature = "mx")]
use aether_fpx::ocp::fp6::Fp6E3M2;
#[cfg(feature = "mx")]
use aether_fpx::ocp::fp8::{Fp8E4M3, Fp8E5M2};
use aether_fpx::trig::*;

#[derive(Deserialize)]
#[serde(rename_all = "snake_case")]
enum FunctionKind {
    Sin,
    Cos,
    Tan,
    Asin,
    Acos,
    Atan,
    Sinh,
    Cosh,
    Tanh,
    Exp,
    Expm1,
    Log,
    Log1p,
}

#[derive(Clone, Copy, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum Backend {
    Core,
    Libm,
    Fpx,
}

#[derive(Clone, Copy, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum Precision {
    F32,
    F64,
    F128,
    #[serde(alias = "fp4_e2m1")]
    Fp4E2M1,
    #[serde(alias = "fp6_e3m2")]
    Fp6E3M2,
    #[serde(alias = "fp8_e4m3")]
    Fp8E4M3,
    #[serde(alias = "fp8_e5m2")]
    Fp8E5M2,
}

#[derive(Clone, Copy, Deserialize)]
struct SeriesCfg {
    backend: Backend,
    precision: Precision,
}

#[derive(Deserialize)]
struct FuncCfg {
    name: String,
    kind: FunctionKind,
    x_min: f64,
    x_max: f64,
    samples: Option<usize>,
    baseline: Option<SeriesCfg>,
    series: Option<Vec<SeriesCfg>>,
}

#[derive(Deserialize)]
struct Config {
    functions: Vec<FuncCfg>,
}

struct SampleRow {
    x: f64,
    values: Vec<EvalValue>,
}

#[derive(Clone, Copy, Debug)]
enum EvalValue {
    F32(f32),
    F64(f64),
    F128(f128),
}

impl EvalValue {
    #[inline]
    fn as_f64_for_plot(self) -> f64 {
        match self {
            EvalValue::F32(v) => v as f64,
            EvalValue::F64(v) => v,
            EvalValue::F128(v) => v as f64,
        }
    }

    #[inline]
    fn as_f32_rounded(self) -> f32 {
        match self {
            EvalValue::F32(v) => v,
            EvalValue::F64(v) => v as f32,
            EvalValue::F128(v) => v as f32,
        }
    }

    #[inline]
    fn as_f64_rounded(self) -> f64 {
        match self {
            EvalValue::F32(v) => v as f64,
            EvalValue::F64(v) => v,
            EvalValue::F128(v) => v as f64,
        }
    }

    #[inline]
    fn as_f128_rounded(self) -> f128 {
        match self {
            EvalValue::F32(v) => v as f128,
            EvalValue::F64(v) => v as f128,
            EvalValue::F128(v) => v,
        }
    }
}

fn default_series() -> Vec<SeriesCfg> {
    let mut series = vec![
        SeriesCfg {
            backend: Backend::Core,
            precision: Precision::F32,
        },
        SeriesCfg {
            backend: Backend::Libm,
            precision: Precision::F32,
        },
        SeriesCfg {
            backend: Backend::Fpx,
            precision: Precision::F32,
        },
        SeriesCfg {
            backend: Backend::Core,
            precision: Precision::F64,
        },
        SeriesCfg {
            backend: Backend::Libm,
            precision: Precision::F64,
        },
        SeriesCfg {
            backend: Backend::Fpx,
            precision: Precision::F64,
        },
        SeriesCfg {
            backend: Backend::Core,
            precision: Precision::F128,
        },
        SeriesCfg {
            backend: Backend::Fpx,
            precision: Precision::F128,
        },
    ];

    #[cfg(feature = "mx")]
    {
        series.extend([
            SeriesCfg {
                backend: Backend::Fpx,
                precision: Precision::Fp4E2M1,
            },
            SeriesCfg {
                backend: Backend::Fpx,
                precision: Precision::Fp6E3M2,
            },
            SeriesCfg {
                backend: Backend::Fpx,
                precision: Precision::Fp8E4M3,
            },
            SeriesCfg {
                backend: Backend::Fpx,
                precision: Precision::Fp8E5M2,
            },
        ]);
    }

    series
}

fn default_baseline() -> SeriesCfg {
    SeriesCfg {
        backend: Backend::Core,
        precision: Precision::F128,
    }
}

fn backend_name(backend: Backend) -> &'static str {
    match backend {
        Backend::Core => "core",
        Backend::Libm => "libm",
        Backend::Fpx => "fpx",
    }
}

fn precision_name(precision: Precision) -> &'static str {
    match precision {
        Precision::F32 => "f32",
        Precision::F64 => "f64",
        Precision::F128 => "f128",
        Precision::Fp4E2M1 => "fp4_e2m1",
        Precision::Fp6E3M2 => "fp6_e3m2",
        Precision::Fp8E4M3 => "fp8_e4m3",
        Precision::Fp8E5M2 => "fp8_e5m2",
    }
}

fn comparison_precision(precision: Precision) -> Precision {
    match precision {
        Precision::Fp4E2M1 | Precision::Fp6E3M2 | Precision::Fp8E4M3 | Precision::Fp8E5M2 => {
            Precision::F32
        }
        other => other,
    }
}

fn series_label(series: SeriesCfg) -> String {
    format!(
        "{} {}",
        backend_name(series.backend),
        precision_name(series.precision)
    )
}

fn csv_label(series: SeriesCfg) -> String {
    format!(
        "{}_{}",
        backend_name(series.backend),
        precision_name(series.precision)
    )
}

fn baseline_label(baseline: SeriesCfg) -> String {
    format!(
        "{} {}",
        backend_name(baseline.backend),
        precision_name(baseline.precision)
    )
}

fn supports_precision(backend: Backend, precision: Precision) -> bool {
    match (backend, precision) {
        (Backend::Libm, Precision::F128) => false,
        (Backend::Core, Precision::Fp4E2M1) => false,
        (Backend::Core, Precision::Fp6E3M2) => false,
        (Backend::Core, Precision::Fp8E4M3) => false,
        (Backend::Core, Precision::Fp8E5M2) => false,
        (Backend::Libm, Precision::Fp4E2M1) => false,
        (Backend::Libm, Precision::Fp6E3M2) => false,
        (Backend::Libm, Precision::Fp8E4M3) => false,
        (Backend::Libm, Precision::Fp8E5M2) => false,
        #[cfg(not(feature = "mx"))]
        (Backend::Fpx, Precision::Fp4E2M1) => false,
        #[cfg(not(feature = "mx"))]
        (Backend::Fpx, Precision::Fp6E3M2) => false,
        #[cfg(not(feature = "mx"))]
        (Backend::Fpx, Precision::Fp8E4M3) => false,
        #[cfg(not(feature = "mx"))]
        (Backend::Fpx, Precision::Fp8E5M2) => false,
        _ => true,
    }
}

fn matched_baseline_for_precision(baseline: SeriesCfg, precision: Precision) -> SeriesCfg {
    let candidate = SeriesCfg {
        backend: baseline.backend,
        precision: comparison_precision(precision),
    };

    if supports_precision(candidate.backend, candidate.precision) {
        candidate
    } else {
        baseline
    }
}

fn matched_baseline_for_series(baseline: SeriesCfg, series: SeriesCfg) -> SeriesCfg {
    matched_baseline_for_precision(baseline, series.precision)
}

fn precision_groups(series_cfgs: &[SeriesCfg]) -> Vec<Precision> {
    let mut out = Vec::new();

    for precision in [
        Precision::F32,
        Precision::F64,
        Precision::F128,
        Precision::Fp4E2M1,
        Precision::Fp6E3M2,
        Precision::Fp8E4M3,
        Precision::Fp8E5M2,
    ] {
        if series_cfgs
            .iter()
            .any(|series| series.precision == precision)
        {
            out.push(precision);
        }
    }

    out
}

fn eval_core_f32(kind: &FunctionKind, x: f32) -> f32 {
    match kind {
        FunctionKind::Sin => x.sin(),
        FunctionKind::Cos => x.cos(),
        FunctionKind::Tan => x.tan(),
        FunctionKind::Asin => x.asin(),
        FunctionKind::Acos => x.acos(),
        FunctionKind::Atan => x.atan(),
        FunctionKind::Sinh => x.sinh(),
        FunctionKind::Cosh => x.cosh(),
        FunctionKind::Tanh => x.tanh(),
        FunctionKind::Exp => x.exp(),
        FunctionKind::Expm1 => x.exp_m1(),
        FunctionKind::Log => x.ln(),
        FunctionKind::Log1p => x.ln_1p(),
    }
}

fn eval_core_f64(kind: &FunctionKind, x: f64) -> f64 {
    match kind {
        FunctionKind::Sin => x.sin(),
        FunctionKind::Cos => x.cos(),
        FunctionKind::Tan => x.tan(),
        FunctionKind::Asin => x.asin(),
        FunctionKind::Acos => x.acos(),
        FunctionKind::Atan => x.atan(),
        FunctionKind::Sinh => x.sinh(),
        FunctionKind::Cosh => x.cosh(),
        FunctionKind::Tanh => x.tanh(),
        FunctionKind::Exp => x.exp(),
        FunctionKind::Expm1 => x.exp_m1(),
        FunctionKind::Log => x.ln(),
        FunctionKind::Log1p => x.ln_1p(),
    }
}

fn eval_core_f128(kind: &FunctionKind, x: f128) -> f128 {
    match kind {
        FunctionKind::Sin => x.sin(),
        FunctionKind::Cos => x.cos(),
        FunctionKind::Tan => x.tan(),
        FunctionKind::Asin => x.asin(),
        FunctionKind::Acos => x.acos(),
        FunctionKind::Atan => x.atan(),
        FunctionKind::Sinh => x.sinh(),
        FunctionKind::Cosh => x.cosh(),
        FunctionKind::Tanh => x.tanh(),
        FunctionKind::Exp => x.exp(),
        FunctionKind::Expm1 => x.exp_m1(),
        FunctionKind::Log => x.ln(),
        FunctionKind::Log1p => x.ln_1p(),
    }
}

fn eval_libm_f32(kind: &FunctionKind, x: f32) -> f32 {
    match kind {
        FunctionKind::Sin => libm::sinf(x),
        FunctionKind::Cos => libm::cosf(x),
        FunctionKind::Tan => libm::tanf(x),
        FunctionKind::Asin => libm::asinf(x),
        FunctionKind::Acos => libm::acosf(x),
        FunctionKind::Atan => libm::atanf(x),
        FunctionKind::Sinh => libm::sinhf(x),
        FunctionKind::Cosh => libm::coshf(x),
        FunctionKind::Tanh => libm::tanhf(x),
        FunctionKind::Exp => libm::expf(x),
        FunctionKind::Expm1 => libm::expm1f(x),
        FunctionKind::Log => libm::logf(x),
        FunctionKind::Log1p => libm::log1pf(x),
    }
}

fn eval_libm_f64(kind: &FunctionKind, x: f64) -> f64 {
    match kind {
        FunctionKind::Sin => libm::sin(x),
        FunctionKind::Cos => libm::cos(x),
        FunctionKind::Tan => libm::tan(x),
        FunctionKind::Asin => libm::asin(x),
        FunctionKind::Acos => libm::acos(x),
        FunctionKind::Atan => libm::atan(x),
        FunctionKind::Sinh => libm::sinh(x),
        FunctionKind::Cosh => libm::cosh(x),
        FunctionKind::Tanh => libm::tanh(x),
        FunctionKind::Exp => libm::exp(x),
        FunctionKind::Expm1 => libm::expm1(x),
        FunctionKind::Log => libm::log(x),
        FunctionKind::Log1p => libm::log1p(x),
    }
}

fn eval_fpx_f32(kind: &FunctionKind, x: f32) -> f32 {
    match kind {
        FunctionKind::Sin => sinf(x),
        FunctionKind::Cos => cosf(x),
        FunctionKind::Tan => tanf(x),
        FunctionKind::Asin => asinf(x),
        FunctionKind::Acos => acosf(x),
        FunctionKind::Atan => atanf(x),
        FunctionKind::Sinh => sinhf(x),
        FunctionKind::Cosh => coshf(x),
        FunctionKind::Tanh => tanhf(x),
        FunctionKind::Exp => expf(x),
        FunctionKind::Expm1 => expm1f(x),
        FunctionKind::Log => logf(x),
        FunctionKind::Log1p => log1pf(x),
    }
}

fn eval_fpx_f64(kind: &FunctionKind, x: f64) -> f64 {
    match kind {
        FunctionKind::Sin => sin(x),
        FunctionKind::Cos => cos(x),
        FunctionKind::Tan => tan(x),
        FunctionKind::Asin => asin(x),
        FunctionKind::Acos => acos(x),
        FunctionKind::Atan => atan(x),
        FunctionKind::Sinh => sinh(x),
        FunctionKind::Cosh => cosh(x),
        FunctionKind::Tanh => tanh(x),
        FunctionKind::Exp => exp(x),
        FunctionKind::Expm1 => expm1(x),
        FunctionKind::Log => log(x),
        FunctionKind::Log1p => log1p(x),
    }
}

fn eval_fpx_f128(kind: &FunctionKind, x: f128) -> f128 {
    match kind {
        FunctionKind::Sin => sinl(x),
        FunctionKind::Cos => cosl(x),
        FunctionKind::Tan => tanl(x),
        FunctionKind::Asin => asinl(x),
        FunctionKind::Acos => acosl(x),
        FunctionKind::Atan => atanl(x),
        FunctionKind::Sinh => sinhl(x),
        FunctionKind::Cosh => coshl(x),
        FunctionKind::Tanh => tanhl(x),
        FunctionKind::Exp => expl(x),
        FunctionKind::Expm1 => expm1l(x),
        FunctionKind::Log => logl(x),
        FunctionKind::Log1p => log1pl(x),
    }
}

#[cfg(feature = "mx")]
fn eval_fpx_fp4_e2m1(kind: &FunctionKind, x: f32) -> f32 {
    let x = Fp4E2M1::from(x);
    match kind {
        FunctionKind::Sin => sin_fp4_e2m1(x),
        FunctionKind::Cos => cos_fp4_e2m1(x),
        FunctionKind::Tan => tan_fp4_e2m1(x),
        FunctionKind::Asin => asin_fp4_e2m1(x),
        FunctionKind::Acos => acos_fp4_e2m1(x),
        FunctionKind::Atan => atan_fp4_e2m1(x),
        FunctionKind::Sinh => sinh_fp4_e2m1(x),
        FunctionKind::Cosh => cosh_fp4_e2m1(x),
        FunctionKind::Tanh => tanh_fp4_e2m1(x),
        FunctionKind::Exp => exp_fp4_e2m1(x),
        FunctionKind::Expm1 => expm1_fp4_e2m1(x),
        FunctionKind::Log => log_fp4_e2m1(x),
        FunctionKind::Log1p => log1p_fp4_e2m1(x),
    }
    .into()
}

#[cfg(feature = "mx")]
fn eval_fpx_fp6_e3m2(kind: &FunctionKind, x: f32) -> f32 {
    let x = Fp6E3M2::from(x);
    match kind {
        FunctionKind::Sin => sin_fp6_e3m2(x),
        FunctionKind::Cos => cos_fp6_e3m2(x),
        FunctionKind::Tan => tan_fp6_e3m2(x),
        FunctionKind::Asin => asin_fp6_e3m2(x),
        FunctionKind::Acos => acos_fp6_e3m2(x),
        FunctionKind::Atan => atan_fp6_e3m2(x),
        FunctionKind::Sinh => sinh_fp6_e3m2(x),
        FunctionKind::Cosh => cosh_fp6_e3m2(x),
        FunctionKind::Tanh => tanh_fp6_e3m2(x),
        FunctionKind::Exp => exp_fp6_e3m2(x),
        FunctionKind::Expm1 => expm1_fp6_e3m2(x),
        FunctionKind::Log => log_fp6_e3m2(x),
        FunctionKind::Log1p => log1p_fp6_e3m2(x),
    }
    .into()
}

#[cfg(feature = "mx")]
fn eval_fpx_fp8_e4m3(kind: &FunctionKind, x: f32) -> f32 {
    let x = Fp8E4M3::from(x);
    match kind {
        FunctionKind::Sin => sin_fp8_e4m3(x),
        FunctionKind::Cos => cos_fp8_e4m3(x),
        FunctionKind::Tan => tan_fp8_e4m3(x),
        FunctionKind::Asin => asin_fp8_e4m3(x),
        FunctionKind::Acos => acos_fp8_e4m3(x),
        FunctionKind::Atan => atan_fp8_e4m3(x),
        FunctionKind::Sinh => sinh_fp8_e4m3(x),
        FunctionKind::Cosh => cosh_fp8_e4m3(x),
        FunctionKind::Tanh => tanh_fp8_e4m3(x),
        FunctionKind::Exp => exp_fp8_e4m3(x),
        FunctionKind::Expm1 => expm1_fp8_e4m3(x),
        FunctionKind::Log => log_fp8_e4m3(x),
        FunctionKind::Log1p => log1p_fp8_e4m3(x),
    }
    .into()
}

#[cfg(feature = "mx")]
fn eval_fpx_fp8_e5m2(kind: &FunctionKind, x: f32) -> f32 {
    let x = Fp8E5M2::from(x);
    match kind {
        FunctionKind::Sin => sin_fp8_e5m2(x),
        FunctionKind::Cos => cos_fp8_e5m2(x),
        FunctionKind::Tan => tan_fp8_e5m2(x),
        FunctionKind::Asin => asin_fp8_e5m2(x),
        FunctionKind::Acos => acos_fp8_e5m2(x),
        FunctionKind::Atan => atan_fp8_e5m2(x),
        FunctionKind::Sinh => sinh_fp8_e5m2(x),
        FunctionKind::Cosh => cosh_fp8_e5m2(x),
        FunctionKind::Tanh => tanh_fp8_e5m2(x),
        FunctionKind::Exp => exp_fp8_e5m2(x),
        FunctionKind::Expm1 => expm1_fp8_e5m2(x),
        FunctionKind::Log => log_fp8_e5m2(x),
        FunctionKind::Log1p => log1p_fp8_e5m2(x),
    }
    .into()
}

#[cfg(not(feature = "mx"))]
fn eval_fpx_fp4_e2m1(_: &FunctionKind, _: f32) -> f32 {
    panic!("mx support is disabled; rerun with --features mx")
}

#[cfg(not(feature = "mx"))]
fn eval_fpx_fp6_e3m2(_: &FunctionKind, _: f32) -> f32 {
    panic!("mx support is disabled; rerun with --features mx")
}

#[cfg(not(feature = "mx"))]
fn eval_fpx_fp8_e4m3(_: &FunctionKind, _: f32) -> f32 {
    panic!("mx support is disabled; rerun with --features mx")
}

#[cfg(not(feature = "mx"))]
fn eval_fpx_fp8_e5m2(_: &FunctionKind, _: f32) -> f32 {
    panic!("mx support is disabled; rerun with --features mx")
}

fn eval_series_value(kind: &FunctionKind, series: SeriesCfg, x: f64) -> EvalValue {
    match (series.backend, series.precision) {
        (Backend::Core, Precision::F32) => EvalValue::F32(eval_core_f32(kind, x as f32)),
        (Backend::Core, Precision::F64) => EvalValue::F64(eval_core_f64(kind, x)),
        (Backend::Core, Precision::F128) => EvalValue::F128(eval_core_f128(kind, x as f128)),
        (Backend::Core, Precision::Fp4E2M1) => panic!("core fp4_e2m1 is not available"),
        (Backend::Core, Precision::Fp6E3M2) => panic!("core fp6_e3m2 is not available"),
        (Backend::Core, Precision::Fp8E4M3) => panic!("core fp8_e4m3 is not available"),
        (Backend::Core, Precision::Fp8E5M2) => panic!("core fp8_e5m2 is not available"),

        (Backend::Libm, Precision::F32) => EvalValue::F32(eval_libm_f32(kind, x as f32)),
        (Backend::Libm, Precision::F64) => EvalValue::F64(eval_libm_f64(kind, x)),
        (Backend::Libm, Precision::F128) => panic!("libm f128 is not available"),
        (Backend::Libm, Precision::Fp4E2M1) => panic!("libm fp4_e2m1 is not available"),
        (Backend::Libm, Precision::Fp6E3M2) => panic!("libm fp6_e3m2 is not available"),
        (Backend::Libm, Precision::Fp8E4M3) => panic!("libm fp8_e4m3 is not available"),
        (Backend::Libm, Precision::Fp8E5M2) => panic!("libm fp8_e5m2 is not available"),

        (Backend::Fpx, Precision::F32) => EvalValue::F32(eval_fpx_f32(kind, x as f32)),
        (Backend::Fpx, Precision::F64) => EvalValue::F64(eval_fpx_f64(kind, x)),
        (Backend::Fpx, Precision::F128) => EvalValue::F128(eval_fpx_f128(kind, x as f128)),
        (Backend::Fpx, Precision::Fp4E2M1) => EvalValue::F32(eval_fpx_fp4_e2m1(kind, x as f32)),
        (Backend::Fpx, Precision::Fp6E3M2) => EvalValue::F32(eval_fpx_fp6_e3m2(kind, x as f32)),
        (Backend::Fpx, Precision::Fp8E4M3) => EvalValue::F32(eval_fpx_fp8_e4m3(kind, x as f32)),
        (Backend::Fpx, Precision::Fp8E5M2) => EvalValue::F32(eval_fpx_fp8_e5m2(kind, x as f32)),
    }
}

fn baseline_input_for_series(series: SeriesCfg, x: f64) -> f64 {
    #[cfg(feature = "mx")]
    {
        let quantized = match series.precision {
            Precision::Fp4E2M1 => {
                let value: f32 = Fp4E2M1::from(x as f32).into();
                value
            }
            Precision::Fp6E3M2 => {
                let value: f32 = Fp6E3M2::from(x as f32).into();
                value
            }
            Precision::Fp8E4M3 => {
                let value: f32 = Fp8E4M3::from(x as f32).into();
                value
            }
            Precision::Fp8E5M2 => {
                let value: f32 = Fp8E5M2::from(x as f32).into();
                value
            }
            _ => return x,
        };

        return quantized as f64;
    }

    #[cfg(not(feature = "mx"))]
    {
        let _ = series;
        x
    }
}

fn eval_baseline_value(
    kind: &FunctionKind,
    baseline: SeriesCfg,
    series: SeriesCfg,
    x: f64,
) -> EvalValue {
    let matched_baseline = matched_baseline_for_series(baseline, series);
    let value = eval_series_value(kind, matched_baseline, baseline_input_for_series(series, x));

    #[cfg(feature = "mx")]
    {
        let rounded = value.as_f32_rounded();
        return match series.precision {
            Precision::Fp4E2M1 => EvalValue::F32(Fp4E2M1::from(rounded).into()),
            Precision::Fp6E3M2 => EvalValue::F32(Fp6E3M2::from(rounded).into()),
            Precision::Fp8E4M3 => EvalValue::F32(Fp8E4M3::from(rounded).into()),
            Precision::Fp8E5M2 => EvalValue::F32(Fp8E5M2::from(rounded).into()),
            _ => value,
        };
    }

    #[cfg(not(feature = "mx"))]
    value
}

fn line_style(series: SeriesCfg) -> PlotStyle {
    match (series.backend, series.precision) {
        (_, Precision::F128) => PlotStyle::Line,
        (Backend::Fpx, _) => PlotStyle::Line,
        _ => PlotStyle::Points,
    }
}

#[cfg(feature = "mx")]
fn sign_magnitude_ulp_dist(a: u8, b: u8, sign_mask: u8) -> u8 {
    let magnitude_mask = sign_mask - 1;
    let a_magnitude = a & magnitude_mask;
    let b_magnitude = b & magnitude_mask;

    if a_magnitude == 0 && b_magnitude == 0 {
        return 0;
    }

    if (a & sign_mask) != (b & sign_mask) {
        a_magnitude.saturating_add(b_magnitude)
    } else {
        a_magnitude.abs_diff(b_magnitude)
    }
}

#[cfg(feature = "mx")]
fn native_mx_ulp_error(precision: Precision, got: f32, baseline: f32) -> f64 {
    if got.is_nan() || baseline.is_nan() {
        return f64::NAN;
    }

    let distance = match precision {
        Precision::Fp4E2M1 => sign_magnitude_ulp_dist(
            Fp4E2M1::from(got).to_bits(),
            Fp4E2M1::from(baseline).to_bits(),
            0x08,
        ),
        Precision::Fp6E3M2 => sign_magnitude_ulp_dist(
            Fp6E3M2::from(got).to_bits(),
            Fp6E3M2::from(baseline).to_bits(),
            0x20,
        ),
        Precision::Fp8E4M3 => sign_magnitude_ulp_dist(
            Fp8E4M3::from(got).to_bits(),
            Fp8E4M3::from(baseline).to_bits(),
            0x80,
        ),
        Precision::Fp8E5M2 => sign_magnitude_ulp_dist(
            Fp8E5M2::from(got).to_bits(),
            Fp8E5M2::from(baseline).to_bits(),
            0x80,
        ),
        _ => unreachable!("native_mx_ulp_error requires an MX precision"),
    };

    distance as f64
}

fn ulp_error(series: SeriesCfg, got: EvalValue, baseline: EvalValue) -> f64 {
    match series.precision {
        Precision::F32 => valid_ulp_dist_f32(got.as_f32_rounded(), baseline.as_f32_rounded())
            .map(|value| value as f64)
            .unwrap_or(f64::NAN),

        #[cfg(feature = "mx")]
        Precision::Fp4E2M1 | Precision::Fp6E3M2 | Precision::Fp8E4M3 | Precision::Fp8E5M2 => {
            native_mx_ulp_error(
                series.precision,
                got.as_f32_rounded(),
                baseline.as_f32_rounded(),
            )
        }

        #[cfg(not(feature = "mx"))]
        Precision::Fp4E2M1 | Precision::Fp6E3M2 | Precision::Fp8E4M3 | Precision::Fp8E5M2 => {
            f64::NAN
        }

        Precision::F64 => valid_ulp_dist_f64(got.as_f64_rounded(), baseline.as_f64_rounded())
            .map(|value| value as f64)
            .unwrap_or(f64::NAN),

        Precision::F128 => valid_ulp_dist_f128(got.as_f128_rounded(), baseline.as_f128_rounded())
            .map(|value| value as f64)
            .unwrap_or(f64::NAN),
    }
}

fn abs_error_for_precision(series: SeriesCfg, got: EvalValue, baseline: EvalValue) -> f64 {
    match series.precision {
        Precision::F32
        | Precision::Fp4E2M1
        | Precision::Fp6E3M2
        | Precision::Fp8E4M3
        | Precision::Fp8E5M2 => (got.as_f32_rounded() - baseline.as_f32_rounded()).abs() as f64,
        Precision::F64 => (got.as_f64_rounded() - baseline.as_f64_rounded()).abs(),
        Precision::F128 => (got.as_f128_rounded() - baseline.as_f128_rounded()).abs() as f64,
    }
}

fn format_eval_value_csv(series: SeriesCfg, value: EvalValue) -> String {
    match series.precision {
        Precision::F32
        | Precision::Fp4E2M1
        | Precision::Fp6E3M2
        | Precision::Fp8E4M3
        | Precision::Fp8E5M2 => format!("{:.9e}", value.as_f32_rounded()),
        Precision::F64 => format!("{:.17e}", value.as_f64_rounded()),
        Precision::F128 => format!("{:?}", value.as_f128_rounded()),
    }
}

fn format_abs_error_csv(series: SeriesCfg, got: EvalValue, baseline: EvalValue) -> String {
    match series.precision {
        Precision::F32
        | Precision::Fp4E2M1
        | Precision::Fp6E3M2
        | Precision::Fp8E4M3
        | Precision::Fp8E5M2 => format!(
            "{:.9e}",
            (got.as_f32_rounded() - baseline.as_f32_rounded()).abs()
        ),
        Precision::F64 => format!(
            "{:.17e}",
            (got.as_f64_rounded() - baseline.as_f64_rounded()).abs()
        ),
        Precision::F128 => format!(
            "{:?}",
            (got.as_f128_rounded() - baseline.as_f128_rounded()).abs()
        ),
    }
}

fn format_eval_bits_csv(series: SeriesCfg, value: EvalValue) -> String {
    match series.precision {
        Precision::F32
        | Precision::Fp4E2M1
        | Precision::Fp6E3M2
        | Precision::Fp8E4M3
        | Precision::Fp8E5M2 => format!("0x{:08x}", value.as_f32_rounded().to_bits()),
        Precision::F64 => format!("0x{:016x}", value.as_f64_rounded().to_bits()),
        Precision::F128 => format!("0x{:032x}", value.as_f128_rounded().to_bits()),
    }
}

fn write_compare_csv(
    out_dir: &str,
    function_name: &str,
    series_cfgs: &[SeriesCfg],
    rows: &[SampleRow],
) -> String {
    let mut csv = String::from("x_f64_sample,x_f64_bits");

    for series in series_cfgs {
        csv.push(',');
        csv.push_str(&csv_label(*series));
        csv.push(',');
        csv.push_str(&format!("{}_bits", csv_label(*series)));
    }

    csv.push('\n');

    for row in rows {
        csv.push_str(&format!("{:.17e}", row.x));
        csv.push_str(&format!(",0x{:016x}", row.x.to_bits()));

        for (idx, value) in row.values.iter().copied().enumerate() {
            csv.push(',');
            csv.push_str(&format_eval_value_csv(series_cfgs[idx], value));
            csv.push(',');
            csv.push_str(&format_eval_bits_csv(series_cfgs[idx], value));
        }

        csv.push('\n');
    }

    let csv_path = format!("{}/{}_compare.csv", out_dir, function_name);
    fs::write(&csv_path, csv).expect("failed to write compare csv");
    csv_path
}

fn write_abs_error_csv(
    out_dir: &str,
    function_name: &str,
    kind: &FunctionKind,
    baseline: SeriesCfg,
    series_cfgs: &[SeriesCfg],
    rows: &[SampleRow],
) -> String {
    let mut csv = String::from("x_f64_sample,x_f64_bits");

    for series in series_cfgs {
        let matched_baseline = matched_baseline_for_series(baseline, *series);
        csv.push(',');
        csv.push_str(&format!(
            "abs_err_{}_vs_{}",
            csv_label(*series),
            csv_label(matched_baseline),
        ));
    }

    csv.push('\n');

    for row in rows {
        csv.push_str(&format!("{:.17e}", row.x));
        csv.push_str(&format!(",0x{:016x}", row.x.to_bits()));

        for (idx, value) in row.values.iter().copied().enumerate() {
            let baseline_value = eval_baseline_value(kind, baseline, series_cfgs[idx], row.x);
            csv.push(',');
            csv.push_str(&format_abs_error_csv(
                series_cfgs[idx],
                value,
                baseline_value,
            ));
        }

        csv.push('\n');
    }

    let csv_path = format!("{}/{}_abs_error.csv", out_dir, function_name);
    fs::write(&csv_path, csv).expect("failed to write abs error csv");
    csv_path
}

fn write_error_csv(
    out_dir: &str,
    function_name: &str,
    suffix: &str,
    column_prefix: &str,
    baseline: SeriesCfg,
    series_cfgs: &[SeriesCfg],
    xs: &[f64],
    series_ys: &[Vec<f64>],
) -> String {
    let mut csv = String::from("x");

    for series in series_cfgs {
        let matched_baseline = matched_baseline_for_series(baseline, *series);
        csv.push(',');
        csv.push_str(&format!(
            "{}_{}_vs_{}",
            column_prefix,
            csv_label(*series),
            csv_label(matched_baseline),
        ));
    }

    csv.push('\n');

    for (row_idx, x) in xs.iter().enumerate() {
        csv.push_str(&format!("{:.17}", x));

        for ys in series_ys {
            csv.push_str(&format!(",{:.17}", ys[row_idx]));
        }

        csv.push('\n');
    }

    let csv_path = format!("{}/{}_{}.csv", out_dir, function_name, suffix);
    fs::write(&csv_path, csv).expect("failed to write error csv");
    csv_path
}

fn render_all_response_plot(
    out_dir: &str,
    function_name: &str,
    xs: &[f64],
    series_cfgs: &[SeriesCfg],
    line_series_ys: &[Vec<f64>],
) -> String {
    let labels: Vec<String> = series_cfgs
        .iter()
        .map(|series| series_label(*series))
        .collect();

    let mut plot_series_vec = Vec::with_capacity(series_cfgs.len());

    for (idx, series) in series_cfgs.iter().copied().enumerate() {
        plot_series_vec.push(
            XYSeries::new(xs, &line_series_ys[idx])
                .with_label(&labels[idx])
                .with_style(line_style(series)),
        );
    }

    let svg_path = format!("{}/response_all.svg", out_dir);

    let title = format!("{} output by implementation and precision", function_name);
    let config = PlotConfig::publication(&title)
        .with_x_label("x")
        .with_y_label("y")
        .with_save_path(&svg_path);
    plot_series_with_config(&plot_series_vec, &config)
        .expect("failed to render all-response svg via aether_viz");

    svg_path
}

fn render_split_precision_plots(
    compare_dir: &str,
    abs_dir: &str,
    ulp_dir: &str,
    summary_dir: &str,
    function_name: &str,
    kind: &FunctionKind,
    xs: &[f64],
    baseline: SeriesCfg,
    series_cfgs: &[SeriesCfg],
    line_series_ys: &[Vec<f64>],
) -> Vec<String> {
    let mut written = Vec::new();

    for precision in precision_groups(series_cfgs) {
        let precision_name_str = precision_name(precision);
        let plot_baseline = matched_baseline_for_precision(baseline, precision);

        let selected: Vec<usize> = series_cfgs
            .iter()
            .enumerate()
            .filter_map(|(idx, series)| {
                if series.precision == precision {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect();

        if selected.is_empty() {
            continue;
        }

        // -----------------------------------------------------------------------------------------
        // Compare plot
        // -----------------------------------------------------------------------------------------
        let compare_labels: Vec<String> = selected
            .iter()
            .map(|&idx| series_label(series_cfgs[idx]))
            .collect();

        let mut compare_plot_series_vec = Vec::with_capacity(selected.len());

        for (label_idx, &idx) in selected.iter().enumerate() {
            let series = series_cfgs[idx];

            compare_plot_series_vec.push(
                XYSeries::new(xs, &line_series_ys[idx])
                    .with_label(&compare_labels[label_idx])
                    .with_style(line_style(series)),
            );
        }

        let svg_compare_path = format!(
            "{}/{}_compare_{}.svg",
            compare_dir, function_name, precision_name_str
        );

        let compare_title = format!(
            "{} {} output — baseline {}",
            function_name,
            precision_name_str,
            baseline_label(plot_baseline)
        );
        let compare_config = PlotConfig::publication(&compare_title)
            .with_x_label("x")
            .with_y_label("y")
            .with_save_path(&svg_compare_path);
        plot_series_with_config(&compare_plot_series_vec, &compare_config)
            .expect("failed to render precision comparison svg via aether_viz");

        written.push(svg_compare_path);

        // -----------------------------------------------------------------------------------------
        // Absolute error plot
        // -----------------------------------------------------------------------------------------
        let err_labels: Vec<String> = selected
            .iter()
            .map(|&idx| series_label(series_cfgs[idx]))
            .collect();

        let err_plot_data: Vec<Vec<f64>> = selected
            .iter()
            .map(|&idx| {
                let series = series_cfgs[idx];
                xs.iter()
                    .map(|&x| {
                        let value = eval_series_value(kind, series, x);
                        let baseline_value = eval_baseline_value(kind, baseline, series, x);
                        abs_error_for_precision(series, value, baseline_value)
                    })
                    .collect()
            })
            .collect();

        let mut err_plot_series_vec = Vec::with_capacity(selected.len());

        for (label_idx, &_idx) in selected.iter().enumerate() {
            err_plot_series_vec.push(
                XYSeries::new(xs, &err_plot_data[label_idx])
                    .with_label(&err_labels[label_idx])
                    .with_style(PlotStyle::Line),
            );
        }

        let svg_abs_path = format!(
            "{}/{}_abs_error_{}.svg",
            abs_dir, function_name, precision_name_str
        );

        let abs_title = format!(
            "{} {} absolute error — vs {}",
            function_name,
            precision_name_str,
            baseline_label(plot_baseline)
        );
        let abs_config = PlotConfig::publication(&abs_title)
            .with_x_label("x")
            .with_y_label("absolute error")
            .with_save_path(&svg_abs_path);
        plot_series_with_config(&err_plot_series_vec, &abs_config)
            .expect("failed to render precision absolute-error svg via aether_viz");

        written.push(svg_abs_path);

        // -----------------------------------------------------------------------------------------
        // ULP error plot
        // -----------------------------------------------------------------------------------------
        let ulp_labels: Vec<String> = selected
            .iter()
            .map(|&idx| series_label(series_cfgs[idx]))
            .collect();

        let ulp_plot_data: Vec<Vec<f64>> = selected
            .iter()
            .map(|&idx| {
                let series = series_cfgs[idx];
                xs.iter()
                    .map(|&x| {
                        let value = eval_series_value(kind, series, x);
                        let baseline_value = eval_baseline_value(kind, baseline, series, x);
                        ulp_error(series, value, baseline_value)
                    })
                    .collect()
            })
            .collect();

        let mut ulp_plot_series_vec = Vec::with_capacity(selected.len());

        for (label_idx, &_idx) in selected.iter().enumerate() {
            ulp_plot_series_vec.push(
                XYSeries::new(xs, &ulp_plot_data[label_idx])
                    .with_label(&ulp_labels[label_idx])
                    .with_style(PlotStyle::Line),
            );
        }

        let svg_ulp_path = format!(
            "{}/{}_ulp_error_{}.svg",
            ulp_dir, function_name, precision_name_str
        );

        let ulp_title = format!(
            "{} {} ULP error — vs {}",
            function_name,
            precision_name_str,
            baseline_label(plot_baseline)
        );
        let ulp_config = PlotConfig::publication(&ulp_title)
            .with_x_label("x")
            .with_y_label("ULP error")
            .with_save_path(&svg_ulp_path);
        plot_series_with_config(&ulp_plot_series_vec, &ulp_config)
            .expect("failed to render precision ULP svg via aether_viz");

        written.push(svg_ulp_path);

        // -----------------------------------------------------------------------------------------
        // Summary plot (absolute + ULP on one figure)
        // -----------------------------------------------------------------------------------------
        let summary_abs_labels: Vec<String> = err_labels
            .iter()
            .map(|label| format!("abs {label}"))
            .collect();
        let summary_ulp_labels: Vec<String> = ulp_labels
            .iter()
            .map(|label| format!("ulp {}", label))
            .collect();

        let mut summary_plot_series_vec = Vec::with_capacity(selected.len() * 2);

        for (label_idx, &_idx) in selected.iter().enumerate() {
            summary_plot_series_vec.push(
                XYSeries::new(xs, &err_plot_data[label_idx])
                    .with_label(&summary_abs_labels[label_idx])
                    .with_style(PlotStyle::Line),
            );

            summary_plot_series_vec.push(
                XYSeries::new(xs, &ulp_plot_data[label_idx])
                    .with_label(&summary_ulp_labels[label_idx])
                    .with_style(PlotStyle::Line),
            );
        }

        let svg_summary_path = format!(
            "{}/{}_summary_error_{}.svg",
            summary_dir, function_name, precision_name_str
        );

        let summary_title = format!(
            "{} {} absolute and ULP error — vs {}",
            function_name,
            precision_name_str,
            baseline_label(plot_baseline)
        );
        let summary_config = PlotConfig::publication(&summary_title)
            .with_x_label("x")
            .with_y_label("error")
            .with_save_path(&svg_summary_path);
        plot_series_with_config(&summary_plot_series_vec, &summary_config)
            .expect("failed to render precision summary svg via aether_viz");

        written.push(svg_summary_path);
    }

    written
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let toml_path = if args.len() > 1 {
        PathBuf::from(&args[1])
    } else {
        PathBuf::from("plot_inputs.toml")
    };

    if !toml_path.exists() {
        eprintln!("Missing config: {}", toml_path.display());
        eprintln!(
            "Usage: cargo run --manifest-path crates/aether-fpx/Cargo.toml --features ocp --bin plot_all -- <config.toml>"
        );
        std::process::exit(1);
    }

    let s = fs::read_to_string(&toml_path).expect("failed to read toml");
    let cfg: Config = toml::from_str(&s).expect("failed to parse toml");

    let out_root = "target/aether-fpx/plots";
    fs::create_dir_all(out_root).expect("failed to create output root directory");

    for f in cfg.functions.iter() {
        let n = f.samples.unwrap_or(1601usize);
        assert!(n >= 2, "plot sample count must be at least 2");

        let baseline = f.baseline.unwrap_or_else(default_baseline);
        let series_cfgs = f.series.clone().unwrap_or_else(default_series);

        let function_dir = format!("{}/{}", out_root, f.name);
        let compare_dir = format!("{}/compare", function_dir);
        let abs_dir = format!("{}/abs", function_dir);
        let ulp_dir = format!("{}/ulp", function_dir);
        let summary_dir = format!("{}/summary", function_dir);

        fs::create_dir_all(&function_dir).expect("failed to create function output directory");
        fs::create_dir_all(&compare_dir).expect("failed to create compare output directory");
        fs::create_dir_all(&abs_dir).expect("failed to create abs output directory");
        fs::create_dir_all(&ulp_dir).expect("failed to create ulp output directory");
        fs::create_dir_all(&summary_dir).expect("failed to create summary output directory");

        let mut xs = Vec::with_capacity(n);
        let mut rows = Vec::with_capacity(n);

        let mut error_series_ys: Vec<Vec<f64>> = (0..series_cfgs.len())
            .map(|_| Vec::with_capacity(n))
            .collect();

        let mut ulp_series_ys: Vec<Vec<f64>> = (0..series_cfgs.len())
            .map(|_| Vec::with_capacity(n))
            .collect();

        let mut line_series_ys: Vec<Vec<f64>> = (0..series_cfgs.len())
            .map(|_| Vec::with_capacity(n))
            .collect();

        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let x = f.x_min + (f.x_max - f.x_min) * t;

            let values: Vec<EvalValue> = series_cfgs
                .iter()
                .map(|series| eval_series_value(&f.kind, *series, x))
                .collect();

            for (idx, value) in values.iter().copied().enumerate() {
                let baseline_value = eval_baseline_value(&f.kind, baseline, series_cfgs[idx], x);
                line_series_ys[idx].push(value.as_f64_for_plot());
                error_series_ys[idx].push(abs_error_for_precision(
                    series_cfgs[idx],
                    value,
                    baseline_value,
                ));
                ulp_series_ys[idx].push(ulp_error(series_cfgs[idx], value, baseline_value));
            }

            xs.push(x);
            rows.push(SampleRow { x, values });
        }

        let csv_compare_path = write_compare_csv(&compare_dir, &f.name, &series_cfgs, &rows);

        let csv_abs_path =
            write_abs_error_csv(&abs_dir, &f.name, &f.kind, baseline, &series_cfgs, &rows);

        let csv_ulp_path = write_error_csv(
            &ulp_dir,
            &f.name,
            "ulp_error",
            "ulp_err",
            baseline,
            &series_cfgs,
            &xs,
            &ulp_series_ys,
        );

        let svg_response_path =
            render_all_response_plot(&function_dir, &f.name, &xs, &series_cfgs, &line_series_ys);

        let svg_paths = render_split_precision_plots(
            &compare_dir,
            &abs_dir,
            &ulp_dir,
            &summary_dir,
            &f.name,
            &f.kind,
            &xs,
            baseline,
            &series_cfgs,
            &line_series_ys,
        );

        println!("Wrote {}", csv_compare_path);
        println!("Wrote {}", csv_abs_path);
        println!("Wrote {}", csv_ulp_path);
        println!("Wrote {}", svg_response_path);

        for path in svg_paths {
            println!("Wrote {}", path);
        }
    }
}

#[cfg(all(test, feature = "mx"))]
mod tests {
    use super::sign_magnitude_ulp_dist;

    #[test]
    fn native_mx_ulp_distance_counts_adjacent_encodings() {
        assert_eq!(sign_magnitude_ulp_dist(0x01, 0x02, 0x08), 1);
        assert_eq!(sign_magnitude_ulp_dist(0x09, 0x0a, 0x08), 1);
        assert_eq!(sign_magnitude_ulp_dist(0x09, 0x01, 0x08), 2);
        assert_eq!(sign_magnitude_ulp_dist(0x08, 0x00, 0x08), 0);
    }
}
