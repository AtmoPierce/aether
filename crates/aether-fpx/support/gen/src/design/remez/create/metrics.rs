#![allow(incomplete_features)]

//! Metrics + acceptance helpers for Remez exports.
//!
//! - ULP acceptance helpers (used by macros + generators)
//! - Abs-error diagnostics in f128 space (used by export/reporting)

use crate::remez::RemezApprox;
use aether_fpx_bench::ulp::{ulp_dist_f128, ulp_dist_f32, ulp_dist_f64};

use core::f128;

pub const MAX_ULP_ALLOWED_F128: u128 = 5;
pub const MAX_ULP_ALLOWED_F64: u64 = 5;
pub const MAX_ULP_ALLOWED_F32: u32 = 5;
pub const ABS_TOL_NEAR_ZERO_F64: f64 = 3.5e-12;
pub const ABS_TOL_NEAR_ZERO_F32: f32 = 8.0e-7;

/* ------------------------------------------------------------------------------------------------
 * Helpers: absolute error eval in f128 space (diagnostics / reporting)
 * ------------------------------------------------------------------------------------------------ */

#[inline(never)]
pub(crate) fn max_abs_errors_for_coeffs_f64<const N: usize, const K: usize, F>(
    rz: &RemezApprox<f128, N, K>,
    coeffs: &[f64],
    f: &F,
    a: f128,
    b: f128,
    samples: usize,
) -> (f64, f64)
where
    F: Fn(f128) -> f128,
{
    let n = samples.max(1);
    let mut max_poly = 0.0_f64;
    let mut max_func = 0.0_f64;

    for j in 0..=n {
        let t = j as f128 / n as f128;
        let x = a + (b - a) * t;

        let p_ref = rz.eval_poly(x); // f128
        let fx = f(x); // f128

        // Horner with coeffs interpreted as f64 -> cast to f128
        let mut p_q = 0.0_f128;
        for &c in coeffs.iter().rev() {
            p_q = p_q * x + (c as f128);
        }

        let e_poly = (p_q - p_ref).abs() as f64;
        let e_func = (p_q - fx).abs() as f64;

        if e_poly > max_poly {
            max_poly = e_poly;
        }
        if e_func > max_func {
            max_func = e_func;
        }
    }

    (max_poly, max_func)
}

#[inline(never)]
pub(crate) fn max_abs_errors_for_coeffs_f32<const N: usize, const K: usize, F>(
    rz: &RemezApprox<f128, N, K>,
    coeffs: &[f32],
    f: &F,
    a: f128,
    b: f128,
    samples: usize,
) -> (f64, f64)
where
    F: Fn(f128) -> f128,
{
    let n = samples.max(1);
    let mut max_poly = 0.0_f64;
    let mut max_func = 0.0_f64;

    for j in 0..=n {
        let t = j as f128 / n as f128;
        let x = a + (b - a) * t;

        let p_ref = rz.eval_poly(x);
        let fx = f(x);

        let mut p_q = 0.0_f128;
        for &c in coeffs.iter().rev() {
            p_q = p_q * x + (c as f128);
        }

        let e_poly = (p_q - p_ref).abs() as f64;
        let e_func = (p_q - fx).abs() as f64;

        if e_poly > max_poly {
            max_poly = e_poly;
        }
        if e_func > max_func {
            max_func = e_func;
        }
    }

    (max_poly, max_func)
}

#[inline(never)]
pub(crate) fn max_abs_errors_for_coeffs_f16<const N: usize, const K: usize, F>(
    rz: &RemezApprox<f128, N, K>,
    coeffs: &[f16; N],
    f: &F,
    a: f128,
    b: f128,
    samples: usize,
) -> (f64, f64)
where
    F: Fn(f128) -> f128,
{
    let n = samples.max(1);
    let mut max_poly = 0.0_f64;
    let mut max_func = 0.0_f64;

    for j in 0..=n {
        let t = j as f128 / n as f128;
        let x = a + (b - a) * t;

        let p_ref = rz.eval_poly(x);
        let fx = f(x);

        let mut p_q = 0.0_f128;
        for &c in coeffs.iter().rev() {
            let c128 = (c as f32) as f128;
            p_q = p_q * x + c128;
        }

        let e_poly = (p_q - p_ref).abs() as f64;
        let e_func = (p_q - fx).abs() as f64;

        if e_poly > max_poly {
            max_poly = e_poly;
        }
        if e_func > max_func {
            max_func = e_func;
        }
    }

    (max_poly, max_func)
}

/* ------------------------------------------------------------------------------------------------
 * Helpers: ULP eval in target float space (this is the *acceptance* metric)
 * NOTE: macros depend on these staying callable under create::... (via re-export).
 * ------------------------------------------------------------------------------------------------ */

#[inline(never)]
pub fn max_ulps_for_coeffs_f128<const N: usize, const K: usize, F>(
    rz: &RemezApprox<f128, N, K>,
    coeffs: &[f128],
    f: &F,
    a: f128,
    b: f128,
    samples: usize,
) -> (u128, u128)
where
    F: Fn(f128) -> f128,
{
    let n = samples.max(1);
    let mut max_poly = 0u128;
    let mut max_func = 0u128;

    for j in 0..=n {
        let t = j as f128 / n as f128;
        let x = a + (b - a) * t;
        let mut pq = 0.0_f128;
        for &c in coeffs.iter().rev() {
            pq = pq * x + c;
        }

        max_poly = max_poly.max(ulp_dist_f128(pq, rz.eval_poly(x)));
        max_func = max_func.max(ulp_dist_f128(pq, f(x)));
    }

    (max_poly, max_func)
}

#[inline(never)]
pub fn max_ulps_for_coeffs_f64<const N: usize, const K: usize, F>(
    rz: &RemezApprox<f128, N, K>,
    coeffs: &[f64],
    f: &F,
    a: f128,
    b: f128,
    samples: usize,
) -> (u64, u64)
where
    F: Fn(f128) -> f128,
{
    let n = samples.max(1);
    let mut max_poly: u64 = 0;
    let mut max_func: u64 = 0;

    for j in 0..=n {
        let t = j as f128 / n as f128;
        let x128 = a + (b - a) * t;

        let x = x128 as f64;
        let mut pq = 0.0_f64;
        for &c in coeffs.iter().rev() {
            pq = pq * x + c;
        }

        let pref = rz.eval_poly(x128) as f64;
        let fx = f(x128) as f64;

        let ulp_poly = ulp_dist_f64(pq, pref);
        if ulp_poly > max_poly {
            max_poly = ulp_poly;
            eprintln!(
                "[remez:f64] new max POLY ulp = {} at x = {:.17e}",
                ulp_poly, x
            );
        }

        let ulp_func = ulp_dist_f64(pq, fx);
        if ulp_func > max_func {
            max_func = ulp_func;
            eprintln!(
                "[remez:f64] new max FUNC ulp = {} at x = {:.17e}",
                ulp_func, x
            );
        }
    }

    (max_poly, max_func)
}

#[inline(never)]
pub fn max_ulps_for_coeffs_f32<const N: usize, const K: usize, F>(
    rz: &RemezApprox<f128, N, K>,
    coeffs: &[f32],
    f: &F,
    a: f128,
    b: f128,
    samples: usize,
) -> (u32, u32)
where
    F: Fn(f128) -> f128,
{
    let n = samples.max(1);
    let mut max_poly: u32 = 0;
    let mut max_func: u32 = 0;

    for j in 0..=n {
        let t = j as f128 / n as f128;
        let x128 = a + (b - a) * t;

        let x = x128 as f32;
        let mut pq = 0.0_f32;
        for &c in coeffs.iter().rev() {
            pq = pq * x + c;
        }

        let pref = rz.eval_poly(x128) as f32;
        let fx = f(x128) as f32;

        let ulp_poly = ulp_dist_f32(pq, pref);
        if ulp_poly > max_poly {
            max_poly = ulp_poly;
            eprintln!(
                "[remez:f32] new max POLY ulp = {} at x = {:.9e}",
                ulp_poly, x
            );
        }

        let ulp_func = ulp_dist_f32(pq, fx);
        if ulp_func > max_func {
            max_func = ulp_func;
            eprintln!(
                "[remez:f32] new max FUNC ulp = {} at x = {:.9e}",
                ulp_func, x
            );
        }
    }

    (max_poly, max_func)
}
