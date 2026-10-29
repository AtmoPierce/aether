use aether_core::real::real::{Real, RealCast};

use crate::coefficients::PolyDomain;
use crate::config::SegmentSpec;
use crate::remez::RemezApprox;

/* ================================================================================================
 * Runtime / solver mapping
 * ================================================================================================ */

pub fn runtime_f_df<T>(name: &str) -> (fn(T) -> T, fn(T) -> T)
where
    T: Real + RealCast + Copy,
{
    match name {
        "sin" => (T::sin, cos_df::<T>),
        "cos" => (T::cos, cos_df::<T>),
        "tan" => (T::tan, tan_df::<T>),

        "sinh" => (T::sinh, sinh_df::<T>),
        "cosh" => (T::cosh, cosh_df::<T>),
        "tanh" => (T::tanh, tanh_df::<T>),

        "asin" => (T::asin, asin_df::<T>),
        "acos" => (T::acos, acos_df::<T>),
        "atan" => (T::atan, atan_df::<T>),

        "log" => (T::ln, log_df::<T>),
        "log1p" => (T::ln_1p, log1p_df::<T>),

        "exp" => (T::exp, T::exp),
        "expm1" => (T::exp_m1, T::exp),

        _ => panic!("no runtime mapping for {name}"),
    }
}

pub fn solve_f_df<T>(dom: PolyDomain, name: &str) -> (fn(T) -> T, fn(T) -> T)
where
    T: Real + RealCast + Copy,
{
    match dom {
        PolyDomain::X => match name {
            "log" => (log_atanh_tail_solve::<T>, log_atanh_tail_solve_df::<T>),
            _ => runtime_f_df::<T>(name),
        },

        PolyDomain::Z => match name {
            "sin" => (zsin_solve::<T>, zsin_solve_df::<T>),
            "cos" => (zcos_solve::<T>, zcos_solve_df::<T>),
            "tan" => (ztan_solve::<T>, ztan_solve_df::<T>),
            _ => panic!("Z-domain solve for non sin/cos"),
        },
    }
}

/* ================================================================================================
 * Derivatives
 * ================================================================================================ */

pub fn cos_df<T>(x: T) -> T
where
    T: Real + RealCast + Copy,
{
    -T::sin(x)
}

pub fn tan_df<T>(x: T) -> T
where
    T: Real + RealCast + Copy,
{
    let c = T::cos(x);
    T::ONE / (c * c)
}

pub fn sinh_df<T>(x: T) -> T
where
    T: Real + RealCast + Copy,
{
    T::cosh(x)
}

pub fn cosh_df<T>(x: T) -> T
where
    T: Real + RealCast + Copy,
{
    T::sinh(x)
}

pub fn tanh_df<T>(x: T) -> T
where
    T: Real + RealCast + Copy,
{
    let t = T::tanh(x);
    T::ONE - t * t
}

pub fn asin_df<T>(x: T) -> T
where
    T: Real + RealCast + Copy,
{
    T::ONE / T::sqrt(T::ONE - x * x)
}

pub fn acos_df<T>(x: T) -> T
where
    T: Real + RealCast + Copy,
{
    -T::ONE / T::sqrt(T::ONE - x * x)
}

pub fn atan_df<T>(x: T) -> T
where
    T: Real + RealCast + Copy,
{
    T::ONE / (T::ONE + x * x)
}

pub fn log_df<T>(x: T) -> T
where
    T: Real + RealCast + Copy,
{
    T::ONE / x
}

pub fn log1p_df<T>(x: T) -> T
where
    T: Real + RealCast + Copy,
{
    T::ONE / (T::ONE + x)
}

/// Tail in the atanh reduction
///
/// `atanh(s) / s = 1 + z Q(z)`, where `z = s^2`.  Fitting `Q` instead of the
/// full quotient fixes the leading one structurally in the generated kernel.
/// The series is used only while designing the polynomial; generated runtime
/// code evaluates the resulting Remez coefficients.
pub fn log_atanh_tail_solve<T>(z: T) -> T
where
    T: Real + RealCast + Copy,
{
    let mut sum = T::ZERO;
    let mut z_pow = T::ONE;

    // Q(z) = 1/3 + z/5 + z^2/7 + ... .  The configured log interval has
    // |z| < 0.03, so 32 terms leave ample margin even for f128 design.
    for n in 0..32 {
        sum = sum + z_pow / T::from_usize(2 * n + 3);
        z_pow = z_pow * z;
    }

    sum
}

pub fn log_atanh_tail_solve_df<T>(z: T) -> T
where
    T: Real + RealCast + Copy,
{
    let mut sum = T::ZERO;
    let mut z_pow = T::ONE;

    // Q'(z) = 1/5 + 2z/7 + 3z^2/9 + ... .
    for n in 0..32 {
        sum = sum + T::from_usize(n + 1) * z_pow / T::from_usize(2 * n + 5);
        z_pow = z_pow * z;
    }

    sum
}

pub fn solve_interval_from_x_for_name<T>(name: &str, dom: PolyDomain, ax: T, bx: T) -> (T, T)
where
    T: Real + RealCast + Copy,
{
    if name == "log" {
        let s_a = (ax - T::ONE) / (ax + T::ONE);
        let s_b = (bx - T::ONE) / (bx + T::ONE);
        let s_a2 = s_a * s_a;
        let s_b2 = s_b * s_b;

        if ax <= T::ONE && bx >= T::ONE {
            (T::ZERO, if s_a2 >= s_b2 { s_a2 } else { s_b2 })
        } else if s_a2 <= s_b2 {
            (s_a2, s_b2)
        } else {
            (s_b2, s_a2)
        }
    } else if name == "tan" {
        let a2 = ax * ax;
        let b2 = bx * bx;

        if ax <= T::ZERO && bx >= T::ZERO {
            (T::ZERO, if a2 >= b2 { a2 } else { b2 })
        } else if a2 <= b2 {
            (a2, b2)
        } else {
            (b2, a2)
        }
    } else {
        solve_interval_from_x(dom, ax, bx)
    }
}

/* ================================================================================================
 * Z-domain sin/cos kernels
 * ================================================================================================ */

pub fn zsin_solve<T>(z: T) -> T
where
    T: Real + RealCast + Copy,
{
    let x = T::sqrt(z);

    if x == T::ZERO {
        T::ONE
    } else {
        T::sin(x) / x
    }
}

pub fn zsin_solve_df<T>(z: T) -> T
where
    T: Real + RealCast + Copy,
{
    let x = T::sqrt(z);

    if x == T::ZERO {
        -T::ONE / T::from_f64(6.0)
    } else {
        (x * T::cos(x) - T::sin(x)) / (T::from_f64(2.0) * x * x * x)
    }
}

pub fn zcos_solve<T>(z: T) -> T
where
    T: Real + RealCast + Copy,
{
    T::cos(T::sqrt(z))
}

pub fn zcos_solve_df<T>(z: T) -> T
where
    T: Real + RealCast + Copy,
{
    let x = T::sqrt(z);

    if x == T::ZERO {
        -T::ONE / T::from_f64(2.0)
    } else {
        -T::sin(x) / (T::from_f64(2.0) * x)
    }
}

pub fn ztan_solve<T>(z: T) -> T
where
    T: Real + RealCast + Copy,
{
    let x = T::sqrt(z);

    if x == T::ZERO {
        T::from_f64(1.0 / 3.0)
    } else {
        (T::tan(x) - x) / (x * x * x)
    }
}

pub fn ztan_solve_df<T>(z: T) -> T
where
    T: Real + RealCast + Copy,
{
    let x = T::sqrt(z);

    if x == T::ZERO {
        T::from_f64(2.0 / 15.0)
    } else {
        let tan_x = T::tan(x);
        let numerator = x * tan_x * tan_x - T::from_f64(3.0) * (tan_x - x);
        numerator / (T::from_f64(2.0) * x * x * x * x * x)
    }
}

/* ================================================================================================
 * Math helpers
 * ================================================================================================ */

pub fn seg_interval<T>(seg: &SegmentSpec) -> (T, T)
where
    T: Real + RealCast + Copy,
{
    (
        T::from_f64(seg.interval[0].as_f64()),
        T::from_f64(seg.interval[1].as_f64()),
    )
}

pub fn seg_interval_f128(seg: &SegmentSpec) -> (f128, f128) {
    seg_interval::<f128>(seg)
}

pub fn solve_interval_from_x<T>(dom: PolyDomain, ax: T, bx: T) -> (T, T)
where
    T: Real + RealCast + Copy,
{
    match dom {
        PolyDomain::X => (ax, bx),

        PolyDomain::Z => {
            let za = ax * ax;
            let zb = bx * bx;

            if za <= zb {
                (za, zb)
            } else {
                (zb, za)
            }
        }
    }
}

/* ================================================================================================
 * Solvers
 * ================================================================================================ */

pub fn solve_remez<T, const N: usize, const K: usize>(
    dom: PolyDomain,
    name: &str,
    interval: (T, T),
) -> RemezApprox<T, N, K>
where
    T: Real + RealCast + Copy,
{
    let (f, df) = solve_f_df::<T>(dom, name);

    let mut rz = RemezApprox::<T, N, K>::new(interval.0, interval.1);
    rz.solve(f, df);
    rz
}

/* ================================================================================================
 * Backward-compatible f128 aliases
 * ================================================================================================ */

pub fn runtime_f_df_f128(name: &str) -> (fn(f128) -> f128, fn(f128) -> f128) {
    runtime_f_df::<f128>(name)
}

pub fn solve_f_df_f128(dom: PolyDomain, name: &str) -> (fn(f128) -> f128, fn(f128) -> f128) {
    solve_f_df::<f128>(dom, name)
}

pub fn solve_remez_f128<const N: usize, const K: usize>(
    dom: PolyDomain,
    name: &str,
    interval: (f128, f128),
) -> RemezApprox<f128, N, K> {
    solve_remez::<f128, N, K>(dom, name, interval)
}
