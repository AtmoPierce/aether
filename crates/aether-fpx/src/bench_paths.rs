//! Batch implementations used to compare non-fused scalar, scalar FMA, and
//! explicitly packed SIMD Horner/Estrin arithmetic in the direct transcendental kernels.

use crate::bench;
use crate::simd;
use crate::trig::{acos, asin, atan, cos, exp, expm1, log, log1p, sin, tan};

pub const OPS: [&str; 13] = [
    "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "exp", "expm1", "log",
    "log1p",
];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ArithmeticPath {
    Scalar,
    Fma,
    Simd,
    SimdEstrin,
}

impl ArithmeticPath {
    pub const ALL: [Self; 4] = [Self::Scalar, Self::Fma, Self::Simd, Self::SimdEstrin];

    pub const fn name(self) -> &'static str {
        match self {
            Self::Scalar => "scalar",
            Self::Fma => "fma",
            Self::Simd => "simd",
            Self::SimdEstrin => "simd_estrin",
        }
    }
}

pub fn kernel_range_f32(op: &str) -> (f32, f32) {
    match op {
        "sin" => (sin::SIN_INTERVAL_MIN_F32, sin::SIN_INTERVAL_MAX_F32),
        "cos" => (cos::COS_INTERVAL_MIN_F32, cos::COS_INTERVAL_MAX_F32),
        "tan" => (tan::TAN_INTERVAL_MIN_F32, tan::TAN_INTERVAL_MAX_F32),
        "asin" => (asin::ASIN_INTERVAL_MIN_F32, asin::ASIN_INTERVAL_MAX_F32),
        "acos" => (acos::ACOS_INTERVAL_MIN_F32, acos::ACOS_INTERVAL_MAX_F32),
        "atan" => (atan::ATAN_INTERVAL_MIN_F32, atan::ATAN_INTERVAL_MAX_F32),
        "sinh" | "cosh" => (0.0, 0.5),
        "tanh" => (0.0, 0.25),
        "exp" => (exp::EXP_INTERVAL_MIN_F32, exp::EXP_INTERVAL_MAX_F32),
        "expm1" => (expm1::EXPM1_INTERVAL_MIN_F32, expm1::EXPM1_INTERVAL_MAX_F32),
        "log" => (log::LOG_INTERVAL_MIN_F32, log::LOG_INTERVAL_MAX_F32),
        "log1p" => (log1p::LOG1P_INTERVAL_MIN_F32, log1p::LOG1P_INTERVAL_MAX_F32),
        _ => panic!("unsupported direct kernel: {op}"),
    }
}

pub fn kernel_range_f64(op: &str) -> (f64, f64) {
    match op {
        "sin" => (sin::SIN_INTERVAL_MIN_F64, sin::SIN_INTERVAL_MAX_F64),
        "cos" => (cos::COS_INTERVAL_MIN_F64, cos::COS_INTERVAL_MAX_F64),
        "tan" => (tan::TAN_INTERVAL_MIN_F64, tan::TAN_INTERVAL_MAX_F64),
        "asin" => (asin::ASIN_INTERVAL_MIN_F64, asin::ASIN_INTERVAL_MAX_F64),
        "acos" => (acos::ACOS_INTERVAL_MIN_F64, acos::ACOS_INTERVAL_MAX_F64),
        "atan" => (atan::ATAN_INTERVAL_MIN_F64, atan::ATAN_INTERVAL_MAX_F64),
        "sinh" | "cosh" => (0.0, 0.125),
        "tanh" => (0.0, 0.0625),
        "exp" => (exp::EXP_INTERVAL_MIN_F64, exp::EXP_INTERVAL_MAX_F64),
        "expm1" => (expm1::EXPM1_INTERVAL_MIN_F64, expm1::EXPM1_INTERVAL_MAX_F64),
        "log" => (log::LOG_INTERVAL_MIN_F64, log::LOG_INTERVAL_MAX_F64),
        "log1p" => (log1p::LOG1P_INTERVAL_MIN_F64, log1p::LOG1P_INTERVAL_MAX_F64),
        _ => panic!("unsupported direct kernel: {op}"),
    }
}

pub(crate) fn segment_f32(op: &str, x: f32) -> usize {
    let first_upper_bound = |bounds: &[f32], value: f32| {
        bounds
            .iter()
            .position(|&bound| value <= bound)
            .unwrap_or(bounds.len() - 1)
    };
    match op {
        "sin" => first_upper_bound(&sin::SIN_SEG_BOUNDS_F32, x.abs()),
        "cos" => first_upper_bound(&cos::COS_SEG_BOUNDS_F32, x.abs()),
        "tan" => first_upper_bound(&tan::TAN_SEG_BOUNDS_F32, x.abs()),
        "asin" => first_upper_bound(&asin::ASIN_SEG_BOUNDS_F32, x),
        "acos" => first_upper_bound(&acos::ACOS_SEG_BOUNDS_F32, x),
        "atan" => first_upper_bound(&atan::ATAN_SEG_BOUNDS_F32, x),
        "exp" => first_upper_bound(&exp::EXP_SEG_BOUNDS_F32, x),
        "expm1" => first_upper_bound(&expm1::EXPM1_SEG_BOUNDS_F32, x),
        "log1p" => first_upper_bound(&log1p::LOG1P_SEG_BOUNDS_F32, x),
        "sinh" | "cosh" | "tanh" | "log" => 0,
        _ => panic!("unsupported direct kernel: {op}"),
    }
}

pub(crate) fn segment_f64(op: &str, x: f64) -> usize {
    let first_upper_bound = |bounds: &[f64], value: f64| {
        bounds
            .iter()
            .position(|&bound| value <= bound)
            .unwrap_or(bounds.len() - 1)
    };
    match op {
        "sin" => first_upper_bound(&sin::SIN_SEG_BOUNDS_F64, x.abs()),
        "cos" => first_upper_bound(&cos::COS_SEG_BOUNDS_F64, x.abs()),
        "tan" => first_upper_bound(&tan::TAN_SEG_BOUNDS_F64, x.abs()),
        "asin" => first_upper_bound(&asin::ASIN_SEG_BOUNDS_F64, x),
        "acos" => first_upper_bound(&acos::ACOS_SEG_BOUNDS_F64, x),
        "atan" => first_upper_bound(&atan::ATAN_SEG_BOUNDS_F64, x),
        "exp" => first_upper_bound(&exp::EXP_SEG_BOUNDS_F64, x),
        "expm1" => first_upper_bound(&expm1::EXPM1_SEG_BOUNDS_F64, x),
        "log1p" => first_upper_bound(&log1p::LOG1P_SEG_BOUNDS_F64, x),
        "sinh" | "cosh" | "tanh" | "log" => 0,
        _ => panic!("unsupported direct kernel: {op}"),
    }
}

pub(crate) fn segment_count(op: &str) -> usize {
    match op {
        "sin" | "atan" => 5,
        "cos" | "exp" | "expm1" => 6,
        "tan" | "log1p" => 3,
        "asin" => 4,
        "acos" => 2,
        "sinh" | "cosh" | "tanh" | "log" => 1,
        _ => panic!("unsupported direct kernel: {op}"),
    }
}

pub(crate) fn coefficients_f32(op: &str, segment: usize) -> &'static [f32] {
    match (op, segment) {
        ("sin", 0) => &sin::COEFFS32_SEG0,
        ("sin", 1) => &sin::COEFFS32_SEG1,
        ("sin", 2) => &sin::COEFFS32_SEG2,
        ("sin", 3) => &sin::COEFFS32_SEG3,
        ("sin", 4) => &sin::COEFFS32_SEG4,
        ("cos", 0) => &cos::COEFFS32_SEG0,
        ("cos", 1) => &cos::COEFFS32_SEG1,
        ("cos", 2) => &cos::COEFFS32_SEG2,
        ("cos", 3) => &cos::COEFFS32_SEG3,
        ("cos", 4) => &cos::COEFFS32_SEG4,
        ("cos", 5) => &cos::COEFFS32_SEG5,
        ("tan", 0) => &tan::COEFFS32_SEG0,
        ("tan", 1) => &tan::COEFFS32_SEG1,
        ("tan", 2) => &tan::COEFFS32_SEG2,
        ("asin", 0) => &asin::COEFFS32_SEG0,
        ("asin", 1) => &asin::COEFFS32_SEG1,
        ("asin", 2) => &asin::COEFFS32_SEG2,
        ("asin", 3) => &asin::COEFFS32_SEG3,
        ("acos", 0) => &acos::COEFFS32_SEG0,
        ("acos", 1) => &acos::COEFFS32_SEG1,
        ("atan", 0) => &atan::COEFFS32_SEG0,
        ("atan", 1) => &atan::COEFFS32_SEG1,
        ("atan", 2) => &atan::COEFFS32_SEG2,
        ("atan", 3) => &atan::COEFFS32_SEG3,
        ("atan", 4) => &atan::COEFFS32_SEG4,
        ("sinh", 0) => &bench::SINH_COEFFS_F32,
        ("cosh", 0) => &bench::COSH_COEFFS_F32,
        ("tanh", 0) => &bench::TANH_COEFFS_F32,
        ("exp", 0) => &exp::COEFFS32_SEG0,
        ("exp", 1) => &exp::COEFFS32_SEG1,
        ("exp", 2) => &exp::COEFFS32_SEG2,
        ("exp", 3) => &exp::COEFFS32_SEG3,
        ("exp", 4) => &exp::COEFFS32_SEG4,
        ("exp", 5) => &exp::COEFFS32_SEG5,
        ("expm1", 0) => &expm1::COEFFS32_SEG0,
        ("expm1", 1) => &expm1::COEFFS32_SEG1,
        ("expm1", 2) => &expm1::COEFFS32_SEG2,
        ("expm1", 3) => &expm1::COEFFS32_SEG3,
        ("expm1", 4) => &expm1::COEFFS32_SEG4,
        ("expm1", 5) => &expm1::COEFFS32_SEG5,
        ("log", 0) => &log::COEFFS32,
        ("log1p", 0) => &log1p::COEFFS32_SEG0,
        ("log1p", 1) => &log1p::COEFFS32_SEG1,
        ("log1p", 2) => &log1p::COEFFS32_SEG2,
        _ => panic!("unsupported direct kernel segment: {op}/{segment}"),
    }
}

pub(crate) fn coefficients_f64(op: &str, segment: usize) -> &'static [f64] {
    match (op, segment) {
        ("sin", 0) => &sin::COEFFS64_SEG0,
        ("sin", 1) => &sin::COEFFS64_SEG1,
        ("sin", 2) => &sin::COEFFS64_SEG2,
        ("sin", 3) => &sin::COEFFS64_SEG3,
        ("sin", 4) => &sin::COEFFS64_SEG4,
        ("cos", 0) => &cos::COEFFS64_SEG0,
        ("cos", 1) => &cos::COEFFS64_SEG1,
        ("cos", 2) => &cos::COEFFS64_SEG2,
        ("cos", 3) => &cos::COEFFS64_SEG3,
        ("cos", 4) => &cos::COEFFS64_SEG4,
        ("cos", 5) => &cos::COEFFS64_SEG5,
        ("tan", 0) => &tan::COEFFS64_SEG0,
        ("tan", 1) => &tan::COEFFS64_SEG1,
        ("tan", 2) => &tan::COEFFS64_SEG2,
        ("asin", 0) => &asin::COEFFS64_SEG0,
        ("asin", 1) => &asin::COEFFS64_SEG1,
        ("asin", 2) => &asin::COEFFS64_SEG2,
        ("asin", 3) => &asin::COEFFS64_SEG3,
        ("acos", 0) => &acos::COEFFS64_SEG0,
        ("acos", 1) => &acos::COEFFS64_SEG1,
        ("atan", 0) => &atan::COEFFS64_SEG0,
        ("atan", 1) => &atan::COEFFS64_SEG1,
        ("atan", 2) => &atan::COEFFS64_SEG2,
        ("atan", 3) => &atan::COEFFS64_SEG3,
        ("atan", 4) => &atan::COEFFS64_SEG4,
        ("sinh", 0) => &bench::SINH_COEFFS_F64,
        ("cosh", 0) => &bench::COSH_COEFFS_F64,
        ("tanh", 0) => &bench::TANH_COEFFS_F64,
        ("exp", 0) => &exp::COEFFS64_SEG0,
        ("exp", 1) => &exp::COEFFS64_SEG1,
        ("exp", 2) => &exp::COEFFS64_SEG2,
        ("exp", 3) => &exp::COEFFS64_SEG3,
        ("exp", 4) => &exp::COEFFS64_SEG4,
        ("exp", 5) => &exp::COEFFS64_SEG5,
        ("expm1", 0) => &expm1::COEFFS64_SEG0,
        ("expm1", 1) => &expm1::COEFFS64_SEG1,
        ("expm1", 2) => &expm1::COEFFS64_SEG2,
        ("expm1", 3) => &expm1::COEFFS64_SEG3,
        ("expm1", 4) => &expm1::COEFFS64_SEG4,
        ("expm1", 5) => &expm1::COEFFS64_SEG5,
        ("log", 0) => &log::COEFFS64,
        ("log1p", 0) => &log1p::COEFFS64_SEG0,
        ("log1p", 1) => &log1p::COEFFS64_SEG1,
        ("log1p", 2) => &log1p::COEFFS64_SEG2,
        _ => panic!("unsupported direct kernel segment: {op}/{segment}"),
    }
}

#[inline]
fn eval_simd_segment_f32(op: &str, segment: usize, values: &mut [f32]) {
    match (op, segment) {
        ("sin", 0) => simd::eval_horner_f32_in_place(values, &sin::COEFFS32_SEG0),
        ("sin", 1) => simd::eval_horner_f32_in_place(values, &sin::COEFFS32_SEG1),
        ("sin", 2) => simd::eval_horner_f32_in_place(values, &sin::COEFFS32_SEG2),
        ("sin", 3) => simd::eval_horner_f32_in_place(values, &sin::COEFFS32_SEG3),
        ("sin", 4) => simd::eval_horner_f32_in_place(values, &sin::COEFFS32_SEG4),
        ("cos", 0) => simd::eval_horner_f32_in_place(values, &cos::COEFFS32_SEG0),
        ("cos", 1) => simd::eval_horner_f32_in_place(values, &cos::COEFFS32_SEG1),
        ("cos", 2) => simd::eval_horner_f32_in_place(values, &cos::COEFFS32_SEG2),
        ("cos", 3) => simd::eval_horner_f32_in_place(values, &cos::COEFFS32_SEG3),
        ("cos", 4) => simd::eval_horner_f32_in_place(values, &cos::COEFFS32_SEG4),
        ("cos", 5) => simd::eval_horner_f32_in_place(values, &cos::COEFFS32_SEG5),
        ("tan", 0) => simd::eval_horner_f32_in_place(values, &tan::COEFFS32_SEG0),
        ("tan", 1) => simd::eval_horner_f32_in_place(values, &tan::COEFFS32_SEG1),
        ("tan", 2) => simd::eval_horner_f32_in_place(values, &tan::COEFFS32_SEG2),
        ("asin", 0) => simd::eval_horner_f32_in_place(values, &asin::COEFFS32_SEG0),
        ("asin", 1) => simd::eval_horner_f32_in_place(values, &asin::COEFFS32_SEG1),
        ("asin", 2) => simd::eval_horner_f32_in_place(values, &asin::COEFFS32_SEG2),
        ("asin", 3) => simd::eval_horner_f32_in_place(values, &asin::COEFFS32_SEG3),
        ("acos", 0) => simd::eval_horner_f32_in_place(values, &acos::COEFFS32_SEG0),
        ("acos", 1) => simd::eval_horner_f32_in_place(values, &acos::COEFFS32_SEG1),
        ("atan", 0) => simd::eval_horner_f32_in_place(values, &atan::COEFFS32_SEG0),
        ("atan", 1) => simd::eval_horner_f32_in_place(values, &atan::COEFFS32_SEG1),
        ("atan", 2) => simd::eval_horner_f32_in_place(values, &atan::COEFFS32_SEG2),
        ("atan", 3) => simd::eval_horner_f32_in_place(values, &atan::COEFFS32_SEG3),
        ("atan", 4) => simd::eval_horner_f32_in_place(values, &atan::COEFFS32_SEG4),
        ("sinh", 0) => simd::eval_horner_f32_in_place(values, &bench::SINH_COEFFS_F32),
        ("cosh", 0) => simd::eval_horner_f32_in_place(values, &bench::COSH_COEFFS_F32),
        ("tanh", 0) => simd::eval_horner_f32_in_place(values, &bench::TANH_COEFFS_F32),
        ("exp", 0) => simd::eval_horner_f32_in_place(values, &exp::COEFFS32_SEG0),
        ("exp", 1) => simd::eval_horner_f32_in_place(values, &exp::COEFFS32_SEG1),
        ("exp", 2) => simd::eval_horner_f32_in_place(values, &exp::COEFFS32_SEG2),
        ("exp", 3) => simd::eval_horner_f32_in_place(values, &exp::COEFFS32_SEG3),
        ("exp", 4) => simd::eval_horner_f32_in_place(values, &exp::COEFFS32_SEG4),
        ("exp", 5) => simd::eval_horner_f32_in_place(values, &exp::COEFFS32_SEG5),
        ("expm1", 0) => simd::eval_horner_f32_in_place(values, &expm1::COEFFS32_SEG0),
        ("expm1", 1) => simd::eval_horner_f32_in_place(values, &expm1::COEFFS32_SEG1),
        ("expm1", 2) => simd::eval_horner_f32_in_place(values, &expm1::COEFFS32_SEG2),
        ("expm1", 3) => simd::eval_horner_f32_in_place(values, &expm1::COEFFS32_SEG3),
        ("expm1", 4) => simd::eval_horner_f32_in_place(values, &expm1::COEFFS32_SEG4),
        ("expm1", 5) => simd::eval_horner_f32_in_place(values, &expm1::COEFFS32_SEG5),
        ("log", 0) => simd::eval_horner_f32_in_place(values, &log::COEFFS32),
        ("log1p", 0) => simd::eval_horner_f32_in_place(values, &log1p::COEFFS32_SEG0),
        ("log1p", 1) => simd::eval_horner_f32_in_place(values, &log1p::COEFFS32_SEG1),
        ("log1p", 2) => simd::eval_horner_f32_in_place(values, &log1p::COEFFS32_SEG2),
        _ => panic!("unsupported direct kernel segment: {op}/{segment}"),
    }
}

#[inline]
fn eval_simd_segment_f64(op: &str, segment: usize, values: &mut [f64]) {
    match (op, segment) {
        ("sin", 0) => simd::eval_horner_f64_in_place(values, &sin::COEFFS64_SEG0),
        ("sin", 1) => simd::eval_horner_f64_in_place(values, &sin::COEFFS64_SEG1),
        ("sin", 2) => simd::eval_horner_f64_in_place(values, &sin::COEFFS64_SEG2),
        ("sin", 3) => simd::eval_horner_f64_in_place(values, &sin::COEFFS64_SEG3),
        ("sin", 4) => simd::eval_horner_f64_in_place(values, &sin::COEFFS64_SEG4),
        ("cos", 0) => simd::eval_horner_f64_in_place(values, &cos::COEFFS64_SEG0),
        ("cos", 1) => simd::eval_horner_f64_in_place(values, &cos::COEFFS64_SEG1),
        ("cos", 2) => simd::eval_horner_f64_in_place(values, &cos::COEFFS64_SEG2),
        ("cos", 3) => simd::eval_horner_f64_in_place(values, &cos::COEFFS64_SEG3),
        ("cos", 4) => simd::eval_horner_f64_in_place(values, &cos::COEFFS64_SEG4),
        ("cos", 5) => simd::eval_horner_f64_in_place(values, &cos::COEFFS64_SEG5),
        ("tan", 0) => simd::eval_horner_f64_in_place(values, &tan::COEFFS64_SEG0),
        ("tan", 1) => simd::eval_horner_f64_in_place(values, &tan::COEFFS64_SEG1),
        ("tan", 2) => simd::eval_horner_f64_in_place(values, &tan::COEFFS64_SEG2),
        ("asin", 0) => simd::eval_horner_f64_in_place(values, &asin::COEFFS64_SEG0),
        ("asin", 1) => simd::eval_horner_f64_in_place(values, &asin::COEFFS64_SEG1),
        ("asin", 2) => simd::eval_horner_f64_in_place(values, &asin::COEFFS64_SEG2),
        ("asin", 3) => simd::eval_horner_f64_in_place(values, &asin::COEFFS64_SEG3),
        ("acos", 0) => simd::eval_horner_f64_in_place(values, &acos::COEFFS64_SEG0),
        ("acos", 1) => simd::eval_horner_f64_in_place(values, &acos::COEFFS64_SEG1),
        ("atan", 0) => simd::eval_horner_f64_in_place(values, &atan::COEFFS64_SEG0),
        ("atan", 1) => simd::eval_horner_f64_in_place(values, &atan::COEFFS64_SEG1),
        ("atan", 2) => simd::eval_horner_f64_in_place(values, &atan::COEFFS64_SEG2),
        ("atan", 3) => simd::eval_horner_f64_in_place(values, &atan::COEFFS64_SEG3),
        ("atan", 4) => simd::eval_horner_f64_in_place(values, &atan::COEFFS64_SEG4),
        ("sinh", 0) => simd::eval_horner_f64_in_place(values, &bench::SINH_COEFFS_F64),
        ("cosh", 0) => simd::eval_horner_f64_in_place(values, &bench::COSH_COEFFS_F64),
        ("tanh", 0) => simd::eval_horner_f64_in_place(values, &bench::TANH_COEFFS_F64),
        ("exp", 0) => simd::eval_horner_f64_in_place(values, &exp::COEFFS64_SEG0),
        ("exp", 1) => simd::eval_horner_f64_in_place(values, &exp::COEFFS64_SEG1),
        ("exp", 2) => simd::eval_horner_f64_in_place(values, &exp::COEFFS64_SEG2),
        ("exp", 3) => simd::eval_horner_f64_in_place(values, &exp::COEFFS64_SEG3),
        ("exp", 4) => simd::eval_horner_f64_in_place(values, &exp::COEFFS64_SEG4),
        ("exp", 5) => simd::eval_horner_f64_in_place(values, &exp::COEFFS64_SEG5),
        ("expm1", 0) => simd::eval_horner_f64_in_place(values, &expm1::COEFFS64_SEG0),
        ("expm1", 1) => simd::eval_horner_f64_in_place(values, &expm1::COEFFS64_SEG1),
        ("expm1", 2) => simd::eval_horner_f64_in_place(values, &expm1::COEFFS64_SEG2),
        ("expm1", 3) => simd::eval_horner_f64_in_place(values, &expm1::COEFFS64_SEG3),
        ("expm1", 4) => simd::eval_horner_f64_in_place(values, &expm1::COEFFS64_SEG4),
        ("expm1", 5) => simd::eval_horner_f64_in_place(values, &expm1::COEFFS64_SEG5),
        ("log", 0) => simd::eval_horner_f64_in_place(values, &log::COEFFS64),
        ("log1p", 0) => simd::eval_horner_f64_in_place(values, &log1p::COEFFS64_SEG0),
        ("log1p", 1) => simd::eval_horner_f64_in_place(values, &log1p::COEFFS64_SEG1),
        ("log1p", 2) => simd::eval_horner_f64_in_place(values, &log1p::COEFFS64_SEG2),
        _ => panic!("unsupported direct kernel segment: {op}/{segment}"),
    }
}

#[inline]
fn eval_simd_estrin_segment_f32(op: &str, segment: usize, values: &mut [f32]) {
    let coefficients = coefficients_f32(op, segment);
    match coefficients.len() {
        6 => simd::eval_estrin_f32_in_place::<6>(values, coefficients.try_into().unwrap()),
        7 => simd::eval_estrin_f32_in_place::<7>(values, coefficients.try_into().unwrap()),
        10 => simd::eval_estrin_f32_in_place::<10>(values, coefficients.try_into().unwrap()),
        12 => simd::eval_estrin_f32_in_place::<12>(values, coefficients.try_into().unwrap()),
        15 => simd::eval_estrin_f32_in_place::<15>(values, coefficients.try_into().unwrap()),
        18 => simd::eval_estrin_f32_in_place::<18>(values, coefficients.try_into().unwrap()),
        21 => simd::eval_estrin_f32_in_place::<21>(values, coefficients.try_into().unwrap()),
        len => panic!("unsupported f32 coefficient count: {len}"),
    }
}

#[inline]
fn eval_simd_estrin_segment_f64(op: &str, segment: usize, values: &mut [f64]) {
    let coefficients = coefficients_f64(op, segment);
    match coefficients.len() {
        6 => simd::eval_estrin_f64_in_place::<6>(values, coefficients.try_into().unwrap()),
        7 => simd::eval_estrin_f64_in_place::<7>(values, coefficients.try_into().unwrap()),
        10 => simd::eval_estrin_f64_in_place::<10>(values, coefficients.try_into().unwrap()),
        12 => simd::eval_estrin_f64_in_place::<12>(values, coefficients.try_into().unwrap()),
        15 => simd::eval_estrin_f64_in_place::<15>(values, coefficients.try_into().unwrap()),
        18 => simd::eval_estrin_f64_in_place::<18>(values, coefficients.try_into().unwrap()),
        21 => simd::eval_estrin_f64_in_place::<21>(values, coefficients.try_into().unwrap()),
        len => panic!("unsupported f64 coefficient count: {len}"),
    }
}

#[inline]
pub(crate) fn eval_simd_path_segment_f32(
    op: &str,
    segment: usize,
    path: ArithmeticPath,
    values: &mut [f32],
) {
    match path {
        ArithmeticPath::Simd => eval_simd_segment_f32(op, segment, values),
        ArithmeticPath::SimdEstrin => eval_simd_estrin_segment_f32(op, segment, values),
        _ => unreachable!("packed evaluator requires a SIMD arithmetic path"),
    }
}

#[inline]
pub(crate) fn eval_simd_path_segment_f64(
    op: &str,
    segment: usize,
    path: ArithmeticPath,
    values: &mut [f64],
) {
    match path {
        ArithmeticPath::Simd => eval_simd_segment_f64(op, segment, values),
        ArithmeticPath::SimdEstrin => eval_simd_estrin_segment_f64(op, segment, values),
        _ => unreachable!("packed evaluator requires a SIMD arithmetic path"),
    }
}

#[inline]
pub(crate) fn polynomial_input_f32(op: &str, x: f32) -> f32 {
    match op {
        "sin" | "cos" | "sinh" | "cosh" | "tanh" => x * x,
        "tan" => x.abs() * x.abs(),
        "log" => {
            let s = (x - 1.0) / (x + 1.0);
            s * s
        }
        _ => x,
    }
}

#[inline]
pub(crate) fn polynomial_input_f64(op: &str, x: f64) -> f64 {
    match op {
        "sin" | "cos" | "sinh" | "cosh" | "tanh" => x * x,
        "tan" => x.abs() * x.abs(),
        "log" => {
            let s = (x - 1.0) / (x + 1.0);
            s * s
        }
        _ => x,
    }
}

#[inline]
pub(crate) fn finish_f32(op: &str, path: ArithmeticPath, x: f32, polynomial_x: f32, p: f32) -> f32 {
    match op {
        "sin" | "sinh" | "tanh" => x * p,
        "cos" | "cosh" | "asin" | "acos" | "atan" | "exp" | "expm1" | "log1p" => p,
        "tan" => {
            let ax = x.abs();
            let t = if path == ArithmeticPath::Scalar {
                ax + ax * polynomial_x * p
            } else {
                (ax * polynomial_x).mul_add(p, ax)
            };
            if x.is_sign_negative() {
                -t
            } else {
                t
            }
        }
        "log" => {
            let s = (x - 1.0) / (x + 1.0);
            let correction = if path == ArithmeticPath::Scalar {
                1.0 + p * polynomial_x
            } else {
                p.mul_add(polynomial_x, 1.0)
            };
            2.0 * s * correction
        }
        _ => panic!("unsupported direct kernel: {op}"),
    }
}

#[inline]
pub(crate) fn finish_f64(op: &str, path: ArithmeticPath, x: f64, polynomial_x: f64, p: f64) -> f64 {
    match op {
        "sin" | "sinh" | "tanh" => x * p,
        "cos" | "cosh" | "asin" | "acos" | "atan" | "exp" | "expm1" | "log1p" => p,
        "tan" => {
            let ax = x.abs();
            let t = if path == ArithmeticPath::Scalar {
                ax + ax * polynomial_x * p
            } else {
                (ax * polynomial_x).mul_add(p, ax)
            };
            if x.is_sign_negative() {
                -t
            } else {
                t
            }
        }
        "log" => {
            let s = (x - 1.0) / (x + 1.0);
            let correction = if path == ArithmeticPath::Scalar {
                1.0 + p * polynomial_x
            } else {
                p.mul_add(polynomial_x, 1.0)
            };
            2.0 * s * correction
        }
        _ => panic!("unsupported direct kernel: {op}"),
    }
}

#[inline]
fn eval_scalar_path_f32(coefficients: &[f32], path: ArithmeticPath, values: &mut [f32]) {
    for value in values {
        let x = *value;
        let mut p = 0.0_f32;
        for &coefficient in coefficients.iter().rev() {
            p = if path == ArithmeticPath::Scalar {
                p * x + coefficient
            } else {
                p.mul_add(x, coefficient)
            };
        }
        *value = p;
    }
}

#[inline]
fn eval_scalar_path_f64(coefficients: &[f64], path: ArithmeticPath, values: &mut [f64]) {
    for value in values {
        let x = *value;
        let mut p = 0.0_f64;
        for &coefficient in coefficients.iter().rev() {
            p = if path == ArithmeticPath::Scalar {
                p * x + coefficient
            } else {
                p.mul_add(x, coefficient)
            };
        }
        *value = p;
    }
}

/// Evaluate one segment's polynomial on inputs that are already in polynomial space.
/// This excludes range selection, packing, scatter, and result reconstruction.
pub fn eval_polynomial_f32(op: &str, segment: usize, path: ArithmeticPath, values: &mut [f32]) {
    match path {
        ArithmeticPath::Scalar | ArithmeticPath::Fma => {
            eval_scalar_path_f32(coefficients_f32(op, segment), path, values)
        }
        ArithmeticPath::Simd | ArithmeticPath::SimdEstrin => {
            eval_simd_path_segment_f32(op, segment, path, values)
        }
    }
}

/// Evaluate one segment's polynomial on inputs that are already in polynomial space.
/// This excludes range selection, packing, scatter, and result reconstruction.
pub fn eval_polynomial_f64(op: &str, segment: usize, path: ArithmeticPath, values: &mut [f64]) {
    match path {
        ArithmeticPath::Scalar | ArithmeticPath::Fma => {
            eval_scalar_path_f64(coefficients_f64(op, segment), path, values)
        }
        ArithmeticPath::Simd | ArithmeticPath::SimdEstrin => {
            eval_simd_path_segment_f64(op, segment, path, values)
        }
    }
}

/// Specialized assembly-inspection entry point for the measured sine segment-zero workload.
/// It delegates to the same scalar and architecture-specific SIMD implementations as
/// [`eval_polynomial_f32`], without retaining the generic operation dispatcher in the dump.
#[doc(hidden)]
#[inline(always)]
pub fn eval_sin_segment0_polynomial_f32(path: ArithmeticPath, values: &mut [f32]) {
    match path {
        ArithmeticPath::Scalar | ArithmeticPath::Fma => {
            eval_scalar_path_f32(&sin::COEFFS32_SEG0, path, values)
        }
        ArithmeticPath::Simd => simd::eval_horner_f32_in_place(values, &sin::COEFFS32_SEG0),
        ArithmeticPath::SimdEstrin => simd::eval_estrin_f32_in_place(values, &sin::COEFFS32_SEG0),
    }
}

/// Binary64 counterpart to [`eval_sin_segment0_polynomial_f32`].
#[doc(hidden)]
#[inline(always)]
pub fn eval_sin_segment0_polynomial_f64(path: ArithmeticPath, values: &mut [f64]) {
    match path {
        ArithmeticPath::Scalar | ArithmeticPath::Fma => {
            eval_scalar_path_f64(&sin::COEFFS64_SEG0, path, values)
        }
        ArithmeticPath::Simd => simd::eval_horner_f64_in_place(values, &sin::COEFFS64_SEG0),
        ArithmeticPath::SimdEstrin => simd::eval_estrin_f64_in_place(values, &sin::COEFFS64_SEG0),
    }
}

const MAX_SEGMENTS: usize = 6;

#[derive(Default)]
pub struct WorkspaceF32 {
    originals: Vec<f32>,
    segments: Vec<u8>,
    grouped_values: Vec<f32>,
    grouped_indices: Vec<usize>,
}

impl WorkspaceF32 {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            originals: Vec::with_capacity(capacity),
            segments: Vec::with_capacity(capacity),
            grouped_values: Vec::with_capacity(capacity),
            grouped_indices: Vec::with_capacity(capacity),
        }
    }
}

#[derive(Default)]
pub struct WorkspaceF64 {
    originals: Vec<f64>,
    segments: Vec<u8>,
    grouped_values: Vec<f64>,
    grouped_indices: Vec<usize>,
}

impl WorkspaceF64 {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            originals: Vec::with_capacity(capacity),
            segments: Vec::with_capacity(capacity),
            grouped_values: Vec::with_capacity(capacity),
            grouped_indices: Vec::with_capacity(capacity),
        }
    }
}

pub fn eval_f32(op: &str, path: ArithmeticPath, values: &mut [f32]) {
    eval_f32_with_workspace(
        op,
        path,
        values,
        &mut WorkspaceF32::with_capacity(values.len()),
    );
}

pub fn eval_f32_with_workspace(
    op: &str,
    path: ArithmeticPath,
    values: &mut [f32],
    workspace: &mut WorkspaceF32,
) {
    workspace.originals.clear();
    workspace.originals.extend_from_slice(values);
    workspace.segments.clear();
    workspace.segments.extend(
        workspace
            .originals
            .iter()
            .map(|&x| segment_f32(op, x) as u8),
    );
    match path {
        ArithmeticPath::Scalar | ArithmeticPath::Fma => {
            for (value, &x) in values.iter_mut().zip(&workspace.originals) {
                *value = polynomial_input_f32(op, x);
            }
            for (index, value) in values.iter_mut().enumerate() {
                let x = *value;
                let mut p = 0.0_f32;
                for &coefficient in coefficients_f32(op, workspace.segments[index] as usize)
                    .iter()
                    .rev()
                {
                    p = if path == ArithmeticPath::Scalar {
                        p * x + coefficient
                    } else {
                        p.mul_add(x, coefficient)
                    };
                }
                *value = p;
            }
            for index in 0..values.len() {
                values[index] = finish_f32(
                    op,
                    path,
                    workspace.originals[index],
                    polynomial_input_f32(op, workspace.originals[index]),
                    values[index],
                );
            }
        }
        ArithmeticPath::Simd | ArithmeticPath::SimdEstrin => {
            let runs = workspace
                .segments
                .windows(2)
                .filter(|pair| pair[0] != pair[1])
                .count()
                + usize::from(!workspace.segments.is_empty());
            if runs.saturating_mul(32) <= values.len() {
                for (value, &original) in values.iter_mut().zip(&workspace.originals) {
                    *value = polynomial_input_f32(op, original);
                }
                let mut start = 0;
                while start < values.len() {
                    let segment = workspace.segments[start];
                    let mut end = start + 1;
                    while end < values.len() && workspace.segments[end] == segment {
                        end += 1;
                    }
                    eval_simd_path_segment_f32(op, segment as usize, path, &mut values[start..end]);
                    start = end;
                }
                for index in 0..values.len() {
                    let original = workspace.originals[index];
                    values[index] = finish_f32(
                        op,
                        path,
                        original,
                        polynomial_input_f32(op, original),
                        values[index],
                    );
                }
                return;
            }

            let segment_count = segment_count(op);
            let mut counts = [0_usize; MAX_SEGMENTS];
            for &segment in &workspace.segments {
                counts[segment as usize] += 1;
            }
            let mut offsets = [0_usize; MAX_SEGMENTS + 1];
            for segment in 0..segment_count {
                offsets[segment + 1] = offsets[segment] + counts[segment];
            }
            let mut cursors = offsets;
            workspace.grouped_values.clear();
            workspace.grouped_values.resize(values.len(), 0.0);
            workspace.grouped_indices.clear();
            workspace.grouped_indices.resize(values.len(), 0);
            for (index, (&original, &segment)) in workspace
                .originals
                .iter()
                .zip(&workspace.segments)
                .enumerate()
            {
                let segment = segment as usize;
                let grouped_index = cursors[segment];
                workspace.grouped_values[grouped_index] = polynomial_input_f32(op, original);
                workspace.grouped_indices[grouped_index] = index;
                cursors[segment] += 1;
            }
            for segment in 0..segment_count {
                eval_simd_path_segment_f32(
                    op,
                    segment,
                    path,
                    &mut workspace.grouped_values[offsets[segment]..offsets[segment + 1]],
                );
            }
            for (&index, &polynomial) in workspace
                .grouped_indices
                .iter()
                .zip(&workspace.grouped_values)
            {
                let original = workspace.originals[index];
                values[index] = finish_f32(
                    op,
                    path,
                    original,
                    polynomial_input_f32(op, original),
                    polynomial,
                );
            }
        }
    }
}

pub fn eval_f64(op: &str, path: ArithmeticPath, values: &mut [f64]) {
    eval_f64_with_workspace(
        op,
        path,
        values,
        &mut WorkspaceF64::with_capacity(values.len()),
    );
}

pub fn eval_f64_with_workspace(
    op: &str,
    path: ArithmeticPath,
    values: &mut [f64],
    workspace: &mut WorkspaceF64,
) {
    workspace.originals.clear();
    workspace.originals.extend_from_slice(values);
    workspace.segments.clear();
    workspace.segments.extend(
        workspace
            .originals
            .iter()
            .map(|&x| segment_f64(op, x) as u8),
    );
    match path {
        ArithmeticPath::Scalar | ArithmeticPath::Fma => {
            for (value, &x) in values.iter_mut().zip(&workspace.originals) {
                *value = polynomial_input_f64(op, x);
            }
            for (index, value) in values.iter_mut().enumerate() {
                let x = *value;
                let mut p = 0.0_f64;
                for &coefficient in coefficients_f64(op, workspace.segments[index] as usize)
                    .iter()
                    .rev()
                {
                    p = if path == ArithmeticPath::Scalar {
                        p * x + coefficient
                    } else {
                        p.mul_add(x, coefficient)
                    };
                }
                *value = p;
            }
            for index in 0..values.len() {
                values[index] = finish_f64(
                    op,
                    path,
                    workspace.originals[index],
                    polynomial_input_f64(op, workspace.originals[index]),
                    values[index],
                );
            }
        }
        ArithmeticPath::Simd | ArithmeticPath::SimdEstrin => {
            let runs = workspace
                .segments
                .windows(2)
                .filter(|pair| pair[0] != pair[1])
                .count()
                + usize::from(!workspace.segments.is_empty());
            if runs.saturating_mul(32) <= values.len() {
                for (value, &original) in values.iter_mut().zip(&workspace.originals) {
                    *value = polynomial_input_f64(op, original);
                }
                let mut start = 0;
                while start < values.len() {
                    let segment = workspace.segments[start];
                    let mut end = start + 1;
                    while end < values.len() && workspace.segments[end] == segment {
                        end += 1;
                    }
                    eval_simd_path_segment_f64(op, segment as usize, path, &mut values[start..end]);
                    start = end;
                }
                for index in 0..values.len() {
                    let original = workspace.originals[index];
                    values[index] = finish_f64(
                        op,
                        path,
                        original,
                        polynomial_input_f64(op, original),
                        values[index],
                    );
                }
                return;
            }

            let segment_count = segment_count(op);
            let mut counts = [0_usize; MAX_SEGMENTS];
            for &segment in &workspace.segments {
                counts[segment as usize] += 1;
            }
            let mut offsets = [0_usize; MAX_SEGMENTS + 1];
            for segment in 0..segment_count {
                offsets[segment + 1] = offsets[segment] + counts[segment];
            }
            let mut cursors = offsets;
            workspace.grouped_values.clear();
            workspace.grouped_values.resize(values.len(), 0.0);
            workspace.grouped_indices.clear();
            workspace.grouped_indices.resize(values.len(), 0);
            for (index, (&original, &segment)) in workspace
                .originals
                .iter()
                .zip(&workspace.segments)
                .enumerate()
            {
                let segment = segment as usize;
                let grouped_index = cursors[segment];
                workspace.grouped_values[grouped_index] = polynomial_input_f64(op, original);
                workspace.grouped_indices[grouped_index] = index;
                cursors[segment] += 1;
            }
            for segment in 0..segment_count {
                eval_simd_path_segment_f64(
                    op,
                    segment,
                    path,
                    &mut workspace.grouped_values[offsets[segment]..offsets[segment + 1]],
                );
            }
            for (&index, &polynomial) in workspace
                .grouped_indices
                .iter()
                .zip(&workspace.grouped_values)
            {
                let original = workspace.originals[index];
                values[index] = finish_f64(
                    op,
                    path,
                    original,
                    polynomial_input_f64(op, original),
                    polynomial,
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scalar_estrin_f32(x: f32, coefficients: &[f32]) -> f32 {
        match coefficients.len() {
            6 => {
                crate::trig::common::eval_estrin_f32_slice::<6>(x, coefficients.try_into().unwrap())
            }
            7 => {
                crate::trig::common::eval_estrin_f32_slice::<7>(x, coefficients.try_into().unwrap())
            }
            10 => crate::trig::common::eval_estrin_f32_slice::<10>(
                x,
                coefficients.try_into().unwrap(),
            ),
            12 => crate::trig::common::eval_estrin_f32_slice::<12>(
                x,
                coefficients.try_into().unwrap(),
            ),
            15 => crate::trig::common::eval_estrin_f32_slice::<15>(
                x,
                coefficients.try_into().unwrap(),
            ),
            18 => crate::trig::common::eval_estrin_f32_slice::<18>(
                x,
                coefficients.try_into().unwrap(),
            ),
            21 => crate::trig::common::eval_estrin_f32_slice::<21>(
                x,
                coefficients.try_into().unwrap(),
            ),
            len => panic!("unsupported f32 coefficient count: {len}"),
        }
    }

    fn scalar_estrin_f64(x: f64, coefficients: &[f64]) -> f64 {
        match coefficients.len() {
            6 => {
                crate::trig::common::eval_estrin_f64_slice::<6>(x, coefficients.try_into().unwrap())
            }
            7 => {
                crate::trig::common::eval_estrin_f64_slice::<7>(x, coefficients.try_into().unwrap())
            }
            10 => crate::trig::common::eval_estrin_f64_slice::<10>(
                x,
                coefficients.try_into().unwrap(),
            ),
            12 => crate::trig::common::eval_estrin_f64_slice::<12>(
                x,
                coefficients.try_into().unwrap(),
            ),
            15 => crate::trig::common::eval_estrin_f64_slice::<15>(
                x,
                coefficients.try_into().unwrap(),
            ),
            18 => crate::trig::common::eval_estrin_f64_slice::<18>(
                x,
                coefficients.try_into().unwrap(),
            ),
            21 => crate::trig::common::eval_estrin_f64_slice::<21>(
                x,
                coefficients.try_into().unwrap(),
            ),
            len => panic!("unsupported f64 coefficient count: {len}"),
        }
    }

    #[test]
    fn packed_and_scalar_fma_paths_match_for_all_direct_kernels() {
        for op in OPS {
            let (lo, hi) = kernel_range_f64(op);
            let input = (0..257)
                .map(|index| lo + (hi - lo) * index as f64 / 256.0)
                .collect::<Vec<_>>();
            let mut fma = input.clone();
            let mut simd = input;
            eval_f64(op, ArithmeticPath::Fma, &mut fma);
            eval_f64(op, ArithmeticPath::Simd, &mut simd);
            assert_eq!(fma, simd, "f64 {op}");

            let (lo, hi) = kernel_range_f32(op);
            let input = (0..257)
                .map(|index| lo + (hi - lo) * index as f32 / 256.0)
                .collect::<Vec<_>>();
            let mut fma = input.clone();
            let mut simd = input;
            eval_f32(op, ArithmeticPath::Fma, &mut fma);
            eval_f32(op, ArithmeticPath::Simd, &mut simd);
            assert_eq!(fma, simd, "f32 {op}");
        }
    }

    #[test]
    fn mixed_segment_partition_matches_scalar_fma() {
        for op in OPS {
            let (lo, hi) = kernel_range_f64(op);
            let input = (0..257)
                .map(|index| lo + (hi - lo) * ((index * 137) % 257) as f64 / 256.0)
                .collect::<Vec<_>>();
            let mut fma = input.clone();
            let mut simd = input;
            eval_f64(op, ArithmeticPath::Fma, &mut fma);
            eval_f64(op, ArithmeticPath::Simd, &mut simd);
            assert_eq!(fma, simd, "mixed f64 {op}");

            let (lo, hi) = kernel_range_f32(op);
            let input = (0..257)
                .map(|index| lo + (hi - lo) * ((index * 137) % 257) as f32 / 256.0)
                .collect::<Vec<_>>();
            let mut fma = input.clone();
            let mut simd = input;
            eval_f32(op, ArithmeticPath::Fma, &mut fma);
            eval_f32(op, ArithmeticPath::Simd, &mut simd);
            assert_eq!(fma, simd, "mixed f32 {op}");
        }
    }

    #[test]
    fn packed_estrin_matches_scalar_estrin_for_every_segment() {
        for op in OPS {
            for segment in 0..segment_count(op) {
                let input_f32 = (0..19).map(|index| index as f32 / 36.0).collect::<Vec<_>>();
                let expected_f32 = input_f32
                    .iter()
                    .map(|&x| scalar_estrin_f32(x, coefficients_f32(op, segment)))
                    .collect::<Vec<_>>();
                let mut packed_f32 = input_f32;
                eval_polynomial_f32(op, segment, ArithmeticPath::SimdEstrin, &mut packed_f32);
                for (&packed, &expected) in packed_f32.iter().zip(&expected_f32) {
                    let tolerance = 4.0 * f32::EPSILON * expected.abs().max(1.0);
                    assert!((packed - expected).abs() <= tolerance, "f32 {op}/{segment}");
                }

                let input_f64 = (0..11).map(|index| index as f64 / 20.0).collect::<Vec<_>>();
                let expected_f64 = input_f64
                    .iter()
                    .map(|&x| scalar_estrin_f64(x, coefficients_f64(op, segment)))
                    .collect::<Vec<_>>();
                let mut packed_f64 = input_f64;
                eval_polynomial_f64(op, segment, ArithmeticPath::SimdEstrin, &mut packed_f64);
                for (&packed, &expected) in packed_f64.iter().zip(&expected_f64) {
                    let tolerance = 4.0 * f64::EPSILON * expected.abs().max(1.0);
                    assert!((packed - expected).abs() <= tolerance, "f64 {op}/{segment}");
                }
            }
        }
    }
}
