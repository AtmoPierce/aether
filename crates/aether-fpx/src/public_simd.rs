//! Public FPX functions with horizontal SIMD polynomial evaluation.
//!
//! Each call accepts one scalar input. Argument reduction and reconstruction
//! follow the public scalar functions. During the polynomial stage, SIMD lanes
//! hold distinct coefficients or partial terms for that one input.

use crate::bench_paths::{self, ArithmeticPath};
use crate::policy::{apply_policy_tan_f32, apply_policy_tan_f64};
use crate::trig;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum HorizontalMethod {
    BlockedHorner,
    Estrin,
}

#[inline]
fn polynomial_array_f32<const N: usize>(
    x: f32,
    coefficients: &[f32; N],
    method: HorizontalMethod,
) -> f32 {
    match method {
        HorizontalMethod::BlockedHorner => crate::simd::eval_blocked_horner_f32(x, coefficients),
        HorizontalMethod::Estrin => crate::simd::eval_horizontal_estrin_f32(x, coefficients),
    }
}

#[inline]
fn polynomial_f32(x: f32, coefficients: &[f32], method: HorizontalMethod) -> f32 {
    match coefficients.len() {
        6 => polynomial_array_f32::<6>(x, coefficients.try_into().unwrap(), method),
        7 => polynomial_array_f32::<7>(x, coefficients.try_into().unwrap(), method),
        10 => polynomial_array_f32::<10>(x, coefficients.try_into().unwrap(), method),
        12 => polynomial_array_f32::<12>(x, coefficients.try_into().unwrap(), method),
        15 => polynomial_array_f32::<15>(x, coefficients.try_into().unwrap(), method),
        18 => polynomial_array_f32::<18>(x, coefficients.try_into().unwrap(), method),
        21 => polynomial_array_f32::<21>(x, coefficients.try_into().unwrap(), method),
        len => panic!("unsupported f32 coefficient count: {len}"),
    }
}

#[inline]
fn polynomial_array_f64<const N: usize>(
    x: f64,
    coefficients: &[f64; N],
    method: HorizontalMethod,
) -> f64 {
    match method {
        HorizontalMethod::BlockedHorner => crate::simd::eval_blocked_horner_f64(x, coefficients),
        HorizontalMethod::Estrin => crate::simd::eval_horizontal_estrin_f64(x, coefficients),
    }
}

#[inline]
fn polynomial_f64(x: f64, coefficients: &[f64], method: HorizontalMethod) -> f64 {
    match coefficients.len() {
        6 => polynomial_array_f64::<6>(x, coefficients.try_into().unwrap(), method),
        7 => polynomial_array_f64::<7>(x, coefficients.try_into().unwrap(), method),
        10 => polynomial_array_f64::<10>(x, coefficients.try_into().unwrap(), method),
        12 => polynomial_array_f64::<12>(x, coefficients.try_into().unwrap(), method),
        15 => polynomial_array_f64::<15>(x, coefficients.try_into().unwrap(), method),
        18 => polynomial_array_f64::<18>(x, coefficients.try_into().unwrap(), method),
        21 => polynomial_array_f64::<21>(x, coefficients.try_into().unwrap(), method),
        len => panic!("unsupported f64 coefficient count: {len}"),
    }
}

#[inline]
fn kernel_f32(op: &str, x: f32, method: HorizontalMethod) -> f32 {
    let segment = bench_paths::segment_f32(op, x);
    let polynomial_x = bench_paths::polynomial_input_f32(op, x);
    let p = polynomial_f32(
        polynomial_x,
        bench_paths::coefficients_f32(op, segment),
        method,
    );
    bench_paths::finish_f32(op, ArithmeticPath::Fma, x, polynomial_x, p)
}

#[inline]
fn kernel_f64(op: &str, x: f64, method: HorizontalMethod) -> f64 {
    let segment = bench_paths::segment_f64(op, x);
    let polynomial_x = bench_paths::polynomial_input_f64(op, x);
    let p = polynomial_f64(
        polynomial_x,
        bench_paths::coefficients_f64(op, segment),
        method,
    );
    bench_paths::finish_f64(op, ArithmeticPath::Fma, x, polynomial_x, p)
}

/// Evaluate one binary32 public function. SIMD is confined to the polynomial
/// coefficient/partial-term stage for this single input.
#[inline]
pub fn eval_f32(op: &str, method: HorizontalMethod, x: f32) -> f32 {
    match op {
        "sin" => sin_cos_f32(x, false, method),
        "cos" => sin_cos_f32(x, true, method),
        "tan" => tan_f32(x, method),
        "asin" => asin_f32(x, method),
        "acos" => acos_f32(x, method),
        "atan" => atan_f32(x, method),
        "sinh" => sinh_f32(x, method),
        "cosh" => cosh_f32(x, method),
        "tanh" => tanh_f32(x, method),
        "exp" => exp_f32(x, method),
        "expm1" => expm1_f32(x, method),
        "log" => log_f32(x, method),
        "log1p" => log1p_f32(x, method),
        _ => panic!("unsupported horizontal SIMD function: {op}"),
    }
}

/// Evaluate one binary64 public function. SIMD is confined to the polynomial
/// coefficient/partial-term stage for this single input.
#[inline]
pub fn eval_f64(op: &str, method: HorizontalMethod, x: f64) -> f64 {
    match op {
        "sin" => sin_cos_f64(x, false, method),
        "cos" => sin_cos_f64(x, true, method),
        "tan" => tan_f64(x, method),
        "asin" => asin_f64(x, method),
        "acos" => acos_f64(x, method),
        "atan" => atan_f64(x, method),
        "sinh" => sinh_f64(x, method),
        "cosh" => cosh_f64(x, method),
        "tanh" => tanh_f64(x, method),
        "exp" => exp_f64(x, method),
        "expm1" => expm1_f64(x, method),
        "log" => log_f64(x, method),
        "log1p" => log1p_f64(x, method),
        _ => panic!("unsupported horizontal SIMD function: {op}"),
    }
}

#[inline]
fn sin_cos_f32(x: f32, cosine: bool, method: HorizontalMethod) -> f32 {
    if !x.is_finite() {
        return f32::NAN;
    }
    let (r, q) = trig::rem_pio2_f32(x);
    match (cosine, q) {
        (false, 0) => kernel_f32("sin", r, method),
        (false, 1) => kernel_f32("cos", r, method),
        (false, 2) => -kernel_f32("sin", r, method),
        (false, 3) => -kernel_f32("cos", r, method),
        (true, 0) => kernel_f32("cos", r, method),
        (true, 1) => -kernel_f32("sin", r, method),
        (true, 2) => -kernel_f32("cos", r, method),
        (true, 3) => kernel_f32("sin", r, method),
        _ => unreachable!(),
    }
}

#[inline]
fn sin_cos_f64(x: f64, cosine: bool, method: HorizontalMethod) -> f64 {
    if !x.is_finite() {
        return f64::NAN;
    }
    let (r, q) = trig::rem_pio2_f64(x);
    match (cosine, q) {
        (false, 0) => kernel_f64("sin", r, method),
        (false, 1) => kernel_f64("cos", r, method),
        (false, 2) => -kernel_f64("sin", r, method),
        (false, 3) => -kernel_f64("cos", r, method),
        (true, 0) => kernel_f64("cos", r, method),
        (true, 1) => -kernel_f64("sin", r, method),
        (true, 2) => -kernel_f64("cos", r, method),
        (true, 3) => kernel_f64("sin", r, method),
        _ => unreachable!(),
    }
}

#[inline]
fn tan_f32(x: f32, method: HorizontalMethod) -> f32 {
    if !x.is_finite() {
        return f32::NAN;
    }
    let (r, q) = trig::rem_pio2_tan_f32(x);
    let t = kernel_f32("tan", apply_policy_tan_f32(r), method);
    trig::tan::reconstruct_tan_f32(q, t)
}

#[inline]
fn tan_f64(x: f64, method: HorizontalMethod) -> f64 {
    if !x.is_finite() {
        return f64::NAN;
    }
    let (r, q) = trig::rem_pio2_tan_f64(x);
    let r = apply_policy_tan_f64(r);
    if q & 1 == 0 && r.abs() <= 0.003 {
        return trig::tan::reconstruct_tan_f64(q, kernel_f64("tan", r, method));
    }
    let s = kernel_f64("sin", r, method);
    let c = kernel_f64("cos", r, method);
    if q & 1 == 0 {
        s / c
    } else {
        -c / s
    }
}

#[inline]
fn asin_f32(x: f32, method: HorizontalMethod) -> f32 {
    if !x.is_finite() || x.abs() > 1.0 {
        return f32::NAN;
    }
    if x == 0.0 {
        return x;
    }
    let ax = x.abs();
    let y = if ax <= trig::asin::ASIN_INTERVAL_MAX_F32 {
        kernel_f32("asin", ax, method)
    } else {
        let t = ((1.0 - ax) * 0.5).sqrt();
        core::f32::consts::FRAC_PI_2 - 2.0 * kernel_f32("asin", t, method)
    };
    if x.is_sign_negative() {
        -y
    } else {
        y
    }
}

#[inline]
fn asin_f64(x: f64, method: HorizontalMethod) -> f64 {
    if !x.is_finite() || x.abs() > 1.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return x;
    }
    if x.abs() == 1.0 {
        return if x.is_sign_negative() {
            -core::f64::consts::FRAC_PI_2
        } else {
            core::f64::consts::FRAC_PI_2
        };
    }
    atan_f64(x / ((1.0 - x) * (1.0 + x)).sqrt(), method)
}

#[inline]
fn acos_f32(x: f32, method: HorizontalMethod) -> f32 {
    if !x.is_finite() || x.abs() > 1.0 {
        return f32::NAN;
    }
    if x.is_sign_negative() {
        return core::f32::consts::PI - acos_f32(-x, method);
    }
    if x <= trig::acos::ACOS_INTERVAL_MAX_F32 {
        return kernel_f32("acos", x, method);
    }
    2.0 * asin_f32(((1.0 - x) * 0.5).sqrt(), method)
}

#[inline]
fn acos_f64(x: f64, method: HorizontalMethod) -> f64 {
    if !x.is_finite() || x.abs() > 1.0 {
        return f64::NAN;
    }
    if x > 0.5 {
        2.0 * asin_f64(((1.0 - x) * 0.5).sqrt(), method)
    } else if x < -0.5 {
        core::f64::consts::PI - 2.0 * asin_f64(((1.0 + x) * 0.5).sqrt(), method)
    } else {
        core::f64::consts::FRAC_PI_2 - asin_f64(x, method)
    }
}

const ATAN_COEFFICIENTS_F32: [f32; 10] = [
    1.0,
    1.0 / 3.0,
    1.0 / 5.0,
    1.0 / 7.0,
    1.0 / 9.0,
    1.0 / 11.0,
    1.0 / 13.0,
    1.0 / 15.0,
    1.0 / 17.0,
    1.0 / 19.0,
];

const ATAN_COEFFICIENTS_F64: [f64; 18] = [
    1.0,
    1.0 / 3.0,
    1.0 / 5.0,
    1.0 / 7.0,
    1.0 / 9.0,
    1.0 / 11.0,
    1.0 / 13.0,
    1.0 / 15.0,
    1.0 / 17.0,
    1.0 / 19.0,
    1.0 / 21.0,
    1.0 / 23.0,
    1.0 / 25.0,
    1.0 / 27.0,
    1.0 / 29.0,
    1.0 / 31.0,
    1.0 / 33.0,
    1.0 / 35.0,
];

#[inline]
fn atan_f32(x: f32, method: HorizontalMethod) -> f32 {
    if x.is_nan() {
        return f32::NAN;
    }
    if x == f32::INFINITY {
        return core::f32::consts::FRAC_PI_2;
    }
    if x == f32::NEG_INFINITY {
        return -core::f32::consts::FRAC_PI_2;
    }
    if x == 0.0 {
        return x;
    }
    let (r, scale, reciprocal) = trig::atan::reduce_atan_positive_f32(x.abs());
    let p = polynomial_f32(-(r * r), &ATAN_COEFFICIENTS_F32, method);
    let y = trig::atan::finish_atan_positive_f32(r * p, scale, reciprocal);
    if x.is_sign_negative() {
        -y
    } else {
        y
    }
}

#[inline]
fn atan_f64(x: f64, method: HorizontalMethod) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == f64::INFINITY {
        return core::f64::consts::FRAC_PI_2;
    }
    if x == f64::NEG_INFINITY {
        return -core::f64::consts::FRAC_PI_2;
    }
    if x == 0.0 {
        return x;
    }
    let (r, scale, reciprocal) = trig::atan::reduce_atan_positive_f64(x.abs());
    let p = polynomial_f64(-(r * r), &ATAN_COEFFICIENTS_F64, method);
    let y = trig::atan::finish_atan_positive_f64(r * p, scale, reciprocal);
    if x.is_sign_negative() {
        -y
    } else {
        y
    }
}

#[inline]
fn sinh_f32(x: f32, method: HorizontalMethod) -> f32 {
    if x.is_nan() || x.is_infinite() || x == 0.0 {
        return x;
    }
    if x.abs() <= 0.5 {
        return kernel_f32("sinh", x, method);
    }
    let e = expm1_f32(x.abs(), method);
    let y = e * (e + 2.0) / (2.0 * (e + 1.0));
    if x.is_sign_negative() {
        -y
    } else {
        y
    }
}

#[inline]
fn sinh_f64(x: f64, method: HorizontalMethod) -> f64 {
    if x.is_nan() || x.is_infinite() || x == 0.0 {
        return x;
    }
    if x.abs() <= 0.125 {
        return kernel_f64("sinh", x, method);
    }
    let e = expm1_f64(x.abs(), method);
    let y = e * (e + 2.0) / (2.0 * (e + 1.0));
    if x.is_sign_negative() {
        -y
    } else {
        y
    }
}

#[inline]
fn cosh_f32(x: f32, method: HorizontalMethod) -> f32 {
    if x.is_nan() {
        return f32::NAN;
    }
    if x.is_infinite() {
        return f32::INFINITY;
    }
    if x.abs() <= 0.5 {
        return kernel_f32("cosh", x, method);
    }
    let e = expm1_f32(x.abs(), method);
    1.0 + e * e / (2.0 * (e + 1.0))
}

#[inline]
fn cosh_f64(x: f64, method: HorizontalMethod) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    if x.abs() <= 0.125 {
        return kernel_f64("cosh", x, method);
    }
    let e = expm1_f64(x.abs(), method);
    1.0 + e * e / (2.0 * (e + 1.0))
}

#[inline]
fn tanh_f32(x: f32, method: HorizontalMethod) -> f32 {
    if x.is_nan() {
        return f32::NAN;
    }
    if x == f32::INFINITY {
        return 1.0;
    }
    if x == f32::NEG_INFINITY {
        return -1.0;
    }
    if x == 0.0 {
        return x;
    }
    if x.abs() <= 0.25 {
        return kernel_f32("tanh", x, method);
    }
    if x.abs() >= 9.0 {
        return if x.is_sign_negative() { -1.0 } else { 1.0 };
    }
    let e = expm1_f32(2.0 * x.abs(), method);
    let y = e / (e + 2.0);
    if x.is_sign_negative() {
        -y
    } else {
        y
    }
}

#[inline]
fn tanh_f64(x: f64, method: HorizontalMethod) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == f64::INFINITY {
        return 1.0;
    }
    if x == f64::NEG_INFINITY {
        return -1.0;
    }
    if x == 0.0 {
        return x;
    }
    if x.abs() <= 0.0625 {
        return kernel_f64("tanh", x, method);
    }
    if x.abs() >= 19.0 {
        return if x.is_sign_negative() { -1.0 } else { 1.0 };
    }
    let e = expm1_f64(2.0 * x.abs(), method);
    let y = e / (e + 2.0);
    if x.is_sign_negative() {
        -y
    } else {
        y
    }
}

#[inline]
fn exp_f32(x: f32, method: HorizontalMethod) -> f32 {
    if x.is_nan() {
        return f32::NAN;
    }
    if x == f32::INFINITY {
        return f32::INFINITY;
    }
    if x == f32::NEG_INFINITY || x < -150.0 * core::f32::consts::LN_2 {
        return 0.0;
    }
    if x > 128.0 * core::f32::consts::LN_2 {
        return f32::INFINITY;
    }
    let (r, k) = trig::exp::reduce_exp_f32(x);
    kernel_f32("exp", r, method) * trig::exp::scale_pow2_f32(k)
}

#[inline]
fn exp_f64(x: f64, method: HorizontalMethod) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == f64::INFINITY {
        return f64::INFINITY;
    }
    if x == f64::NEG_INFINITY || x < -1075.0 * core::f64::consts::LN_2 {
        return 0.0;
    }
    if x > 1024.0 * core::f64::consts::LN_2 {
        return f64::INFINITY;
    }
    let (r, k) = trig::exp::reduce_exp_f64(x);
    kernel_f64("exp", r, method) * trig::exp::scale_pow2_f64(k)
}

#[inline]
fn expm1_f32(x: f32, method: HorizontalMethod) -> f32 {
    if x.is_nan() {
        return f32::NAN;
    }
    if x == f32::INFINITY {
        return f32::INFINITY;
    }
    if x == f32::NEG_INFINITY {
        return -1.0;
    }
    if x == 0.0 || x.abs() <= f32::EPSILON {
        return x;
    }
    if x.abs() <= 0.03125 {
        return trig::expm1f(x);
    }
    if x.abs() <= trig::expm1::EXPM1_INTERVAL_MAX_F32 {
        return kernel_f32("expm1", x, method);
    }
    let (r, k) = trig::exp::reduce_exp_f32(x);
    let scale = trig::exp::scale_pow2_f32(k);
    scale * kernel_f32("expm1", r, method) + (scale - 1.0)
}

#[inline]
fn expm1_f64(x: f64, method: HorizontalMethod) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == f64::INFINITY {
        return f64::INFINITY;
    }
    if x == f64::NEG_INFINITY {
        return -1.0;
    }
    if x == 0.0 || x.abs() <= f64::EPSILON {
        return x;
    }
    if x.abs() <= 0.5 {
        return trig::expm1(x);
    }
    let (r, k) = trig::exp::reduce_exp_f64(x);
    trig::exp::scale_pow2_f64(k) * kernel_f64("exp", r, method) - 1.0
}

#[inline]
fn log_f32(x: f32, method: HorizontalMethod) -> f32 {
    if x.is_nan() || x < 0.0 {
        return f32::NAN;
    }
    if x == 0.0 {
        return f32::NEG_INFINITY;
    }
    if x == f32::INFINITY {
        return f32::INFINITY;
    }
    if x == 1.0 {
        return 0.0;
    }
    let (m, e) = trig::log::reduce_log_f32(x);
    trig::log::reconstruct_log_f32(kernel_f32("log", m, method), e)
}

#[inline]
fn log_f64(x: f64, method: HorizontalMethod) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x == f64::INFINITY {
        return f64::INFINITY;
    }
    if x == 1.0 {
        return 0.0;
    }
    let (m, e) = trig::log::reduce_log_f64(x);
    trig::log::reconstruct_log_f64(kernel_f64("log", m, method), e)
}

#[inline]
fn log1p_f32(x: f32, method: HorizontalMethod) -> f32 {
    if x.is_nan() || x < -1.0 {
        return f32::NAN;
    }
    if x == -1.0 {
        return f32::NEG_INFINITY;
    }
    if x == f32::INFINITY {
        return f32::INFINITY;
    }
    if x == 0.0 || x.abs() <= f32::EPSILON {
        return x;
    }
    if x.abs() <= f32::EPSILON.sqrt() {
        return x - 0.5 * x * x;
    }
    if x >= trig::log1p::LOG1P_INTERVAL_MIN_F32 && x <= trig::log1p::LOG1P_INTERVAL_MAX_F32 {
        return kernel_f32("log1p", x, method);
    }
    let y = 1.0 + x;
    log_f32(y, method) + (x - (y - 1.0)) / y
}

#[inline]
fn log1p_f64(x: f64, method: HorizontalMethod) -> f64 {
    if x.is_nan() || x < -1.0 {
        return f64::NAN;
    }
    if x == -1.0 {
        return f64::NEG_INFINITY;
    }
    if x == f64::INFINITY {
        return f64::INFINITY;
    }
    if x == 0.0 || x.abs() <= f64::EPSILON {
        return x;
    }
    if x.abs() <= 0.5 {
        return trig::log1p(x);
    }
    let y = 1.0 + x;
    log_f64(y, method) + (x - (y - 1.0)) / y
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ordered_f32(value: f32) -> u32 {
        let bits = value.to_bits();
        if bits & (1 << 31) != 0 {
            !bits
        } else {
            bits | (1 << 31)
        }
    }

    fn ordered_f64(value: f64) -> u64 {
        let bits = value.to_bits();
        if bits & (1 << 63) != 0 {
            !bits
        } else {
            bits | (1 << 63)
        }
    }

    fn ulp_f32(a: f32, b: f32) -> u32 {
        if a.is_nan() && b.is_nan() {
            0
        } else {
            ordered_f32(a).abs_diff(ordered_f32(b))
        }
    }

    fn ulp_f64(a: f64, b: f64) -> u64 {
        if a.is_nan() && b.is_nan() {
            0
        } else {
            ordered_f64(a).abs_diff(ordered_f64(b))
        }
    }

    #[test]
    fn horizontal_paths_match_scalar_public_functions() {
        let cases32: [(&str, f32, f32, fn(f32) -> f32); 13] = [
            ("sin", -8.0, 8.0, trig::sinf),
            ("cos", -8.0, 8.0, trig::cosf),
            ("tan", -1.25, 1.25, trig::tanf),
            ("asin", -1.0, 1.0, trig::asinf),
            ("acos", -1.0, 1.0, trig::acosf),
            ("atan", -8.0, 8.0, trig::atanf),
            ("sinh", -4.0, 4.0, trig::sinhf),
            ("cosh", -4.0, 4.0, trig::coshf),
            ("tanh", -4.0, 4.0, trig::tanhf),
            ("exp", -4.0, 4.0, trig::expf),
            ("expm1", -4.0, 4.0, trig::expm1f),
            ("log", 0.25, 16.0, trig::logf),
            ("log1p", -0.875, 8.0, trig::log1pf),
        ];
        let cases64: [(&str, f64, f64, fn(f64) -> f64); 13] = [
            ("sin", -8.0, 8.0, trig::sin),
            ("cos", -8.0, 8.0, trig::cos),
            ("tan", -1.25, 1.25, trig::tan),
            ("asin", -1.0, 1.0, trig::asin),
            ("acos", -1.0, 1.0, trig::acos),
            ("atan", -8.0, 8.0, trig::atan),
            ("sinh", -4.0, 4.0, trig::sinh),
            ("cosh", -4.0, 4.0, trig::cosh),
            ("tanh", -4.0, 4.0, trig::tanh),
            ("exp", -4.0, 4.0, trig::exp),
            ("expm1", -4.0, 4.0, trig::expm1),
            ("log", 0.25, 16.0, trig::log),
            ("log1p", -0.875, 8.0, trig::log1p),
        ];
        for method in [HorizontalMethod::BlockedHorner, HorizontalMethod::Estrin] {
            for (op, lo, hi, scalar) in cases32 {
                let worst = (0..=4096)
                    .map(|i| lo + (hi - lo) * (i as f32 / 4096.0))
                    .map(|x| ulp_f32(eval_f32(op, method, x), scalar(x)))
                    .max()
                    .unwrap();
                assert!(worst <= 6, "{op}/{method:?}: {worst} ULP");
            }
            for (op, lo, hi, scalar) in cases64 {
                let worst = (0..=4096)
                    .map(|i| lo + (hi - lo) * (i as f64 / 4096.0))
                    .map(|x| ulp_f64(eval_f64(op, method, x), scalar(x)))
                    .max()
                    .unwrap();
                assert!(worst <= 6, "{op}/{method:?}: {worst} ULP");
            }
        }
    }
}
