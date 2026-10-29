#![feature(f16, f128)]

use core::{f128, f16};
use std::hint::black_box;

/* ---------------- config ---------------- */

#[derive(Clone, Copy, Debug)]
enum BenchKind {
    Fixed,
    Bounded,
}

const DEFAULT_BENCH_KIND: BenchKind = BenchKind::Bounded;

/* ---------------- aether-fpx wrappers ---------------- */

// f32
#[inline(never)]
fn fpx_sin_f32(x: f32) -> f32 {
    aether_fpx::sinf(x)
}
#[inline(never)]
fn fpx_cos_f32(x: f32) -> f32 {
    aether_fpx::cosf(x)
}
#[inline(never)]
fn fpx_tan_f32(x: f32) -> f32 {
    aether_fpx::tanf(x)
}
#[inline(never)]
fn fpx_asin_f32(x: f32) -> f32 {
    aether_fpx::asinf(x)
}
#[inline(never)]
fn fpx_acos_f32(x: f32) -> f32 {
    aether_fpx::acosf(x)
}
#[inline(never)]
fn fpx_atan_f32(x: f32) -> f32 {
    aether_fpx::atanf(x)
}
#[inline(never)]
fn fpx_sinh_f32(x: f32) -> f32 {
    aether_fpx::sinhf(x)
}
#[inline(never)]
fn fpx_cosh_f32(x: f32) -> f32 {
    aether_fpx::coshf(x)
}
#[inline(never)]
fn fpx_tanh_f32(x: f32) -> f32 {
    aether_fpx::tanhf(x)
}
#[inline(never)]
fn fpx_exp_f32(x: f32) -> f32 {
    aether_fpx::expf(x)
}
#[inline(never)]
fn fpx_log_f32(x: f32) -> f32 {
    aether_fpx::logf(x)
}
#[inline(never)]
fn fpx_log1p_f32(x: f32) -> f32 {
    aether_fpx::log1pf(x)
}
#[inline(never)]
fn fpx_expm1_f32(x: f32) -> f32 {
    aether_fpx::expm1f(x)
}

// f64
#[inline(never)]
fn fpx_sin_f64(x: f64) -> f64 {
    aether_fpx::sin(x)
}
#[inline(never)]
fn fpx_cos_f64(x: f64) -> f64 {
    aether_fpx::cos(x)
}
#[inline(never)]
fn fpx_tan_f64(x: f64) -> f64 {
    aether_fpx::tan(x)
}
#[inline(never)]
fn fpx_asin_f64(x: f64) -> f64 {
    aether_fpx::asin(x)
}
#[inline(never)]
fn fpx_acos_f64(x: f64) -> f64 {
    aether_fpx::acos(x)
}
#[inline(never)]
fn fpx_atan_f64(x: f64) -> f64 {
    aether_fpx::atan(x)
}
#[inline(never)]
fn fpx_sinh_f64(x: f64) -> f64 {
    aether_fpx::sinh(x)
}
#[inline(never)]
fn fpx_cosh_f64(x: f64) -> f64 {
    aether_fpx::cosh(x)
}
#[inline(never)]
fn fpx_tanh_f64(x: f64) -> f64 {
    aether_fpx::tanh(x)
}
#[inline(never)]
fn fpx_exp_f64(x: f64) -> f64 {
    aether_fpx::exp(x)
}
#[inline(never)]
fn fpx_log_f64(x: f64) -> f64 {
    aether_fpx::log(x)
}
#[inline(never)]
fn fpx_log1p_f64(x: f64) -> f64 {
    aether_fpx::log1p(x)
}
#[inline(never)]
fn fpx_expm1_f64(x: f64) -> f64 {
    aether_fpx::expm1(x)
}

macro_rules! kernel_wrapper {
    ($name:ident, $target:path, $ty:ty) => {
        #[inline(never)]
        fn $name(x: $ty) -> $ty {
            $target(x)
        }
    };
}

// f32 kernels
kernel_wrapper!(
    fpx_sin_kernel_horner_f32,
    aether_fpx::bench::sin_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_cos_kernel_horner_f32,
    aether_fpx::bench::cos_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_tan_kernel_horner_f32,
    aether_fpx::bench::tan_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_asin_kernel_horner_f32,
    aether_fpx::bench::asin_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_acos_kernel_horner_f32,
    aether_fpx::bench::acos_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_atan_kernel_horner_f32,
    aether_fpx::bench::atan_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_sinh_kernel_horner_f32,
    aether_fpx::bench::sinh_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_cosh_kernel_horner_f32,
    aether_fpx::bench::cosh_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_tanh_kernel_horner_f32,
    aether_fpx::bench::tanh_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_exp_kernel_horner_f32,
    aether_fpx::bench::exp_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_log_kernel_horner_f32,
    aether_fpx::bench::log_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_log1p_kernel_horner_f32,
    aether_fpx::bench::log1p_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_expm1_kernel_horner_f32,
    aether_fpx::bench::expm1_kernel_horner_f32,
    f32
);
kernel_wrapper!(
    fpx_sin_kernel_estrin_f32,
    aether_fpx::bench::sin_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_cos_kernel_estrin_f32,
    aether_fpx::bench::cos_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_tan_kernel_estrin_f32,
    aether_fpx::bench::tan_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_asin_kernel_estrin_f32,
    aether_fpx::bench::asin_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_acos_kernel_estrin_f32,
    aether_fpx::bench::acos_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_atan_kernel_estrin_f32,
    aether_fpx::bench::atan_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_sinh_kernel_estrin_f32,
    aether_fpx::bench::sinh_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_cosh_kernel_estrin_f32,
    aether_fpx::bench::cosh_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_tanh_kernel_estrin_f32,
    aether_fpx::bench::tanh_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_exp_kernel_estrin_f32,
    aether_fpx::bench::exp_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_log_kernel_estrin_f32,
    aether_fpx::bench::log_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_log1p_kernel_estrin_f32,
    aether_fpx::bench::log1p_kernel_estrin_f32,
    f32
);
kernel_wrapper!(
    fpx_expm1_kernel_estrin_f32,
    aether_fpx::bench::expm1_kernel_estrin_f32,
    f32
);

// f64 kernels
kernel_wrapper!(
    fpx_sin_kernel_horner_f64,
    aether_fpx::bench::sin_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_cos_kernel_horner_f64,
    aether_fpx::bench::cos_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_tan_kernel_horner_f64,
    aether_fpx::bench::tan_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_asin_kernel_horner_f64,
    aether_fpx::bench::asin_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_acos_kernel_horner_f64,
    aether_fpx::bench::acos_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_atan_kernel_horner_f64,
    aether_fpx::bench::atan_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_sinh_kernel_horner_f64,
    aether_fpx::bench::sinh_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_cosh_kernel_horner_f64,
    aether_fpx::bench::cosh_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_tanh_kernel_horner_f64,
    aether_fpx::bench::tanh_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_exp_kernel_horner_f64,
    aether_fpx::bench::exp_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_log_kernel_horner_f64,
    aether_fpx::bench::log_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_log1p_kernel_horner_f64,
    aether_fpx::bench::log1p_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_expm1_kernel_horner_f64,
    aether_fpx::bench::expm1_kernel_horner_f64,
    f64
);
kernel_wrapper!(
    fpx_sin_kernel_estrin_f64,
    aether_fpx::bench::sin_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_cos_kernel_estrin_f64,
    aether_fpx::bench::cos_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_tan_kernel_estrin_f64,
    aether_fpx::bench::tan_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_asin_kernel_estrin_f64,
    aether_fpx::bench::asin_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_acos_kernel_estrin_f64,
    aether_fpx::bench::acos_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_atan_kernel_estrin_f64,
    aether_fpx::bench::atan_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_sinh_kernel_estrin_f64,
    aether_fpx::bench::sinh_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_cosh_kernel_estrin_f64,
    aether_fpx::bench::cosh_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_tanh_kernel_estrin_f64,
    aether_fpx::bench::tanh_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_exp_kernel_estrin_f64,
    aether_fpx::bench::exp_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_log_kernel_estrin_f64,
    aether_fpx::bench::log_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_log1p_kernel_estrin_f64,
    aether_fpx::bench::log1p_kernel_estrin_f64,
    f64
);
kernel_wrapper!(
    fpx_expm1_kernel_estrin_f64,
    aether_fpx::bench::expm1_kernel_estrin_f64,
    f64
);

// f16
#[inline(never)]
fn fpx_sin_f16(x: f16) -> f16 {
    aether_fpx::sins(x)
}
#[inline(never)]
fn fpx_cos_f16(x: f16) -> f16 {
    aether_fpx::coss(x)
}
#[inline(never)]
fn fpx_tan_f16(x: f16) -> f16 {
    aether_fpx::tans(x)
}
#[inline(never)]
fn fpx_asin_f16(x: f16) -> f16 {
    aether_fpx::asins(x)
}
#[inline(never)]
fn fpx_acos_f16(x: f16) -> f16 {
    aether_fpx::acoss(x)
}
#[inline(never)]
fn fpx_atan_f16(x: f16) -> f16 {
    aether_fpx::atans(x)
}
#[inline(never)]
fn fpx_sinh_f16(x: f16) -> f16 {
    aether_fpx::sinhs(x)
}
#[inline(never)]
fn fpx_cosh_f16(x: f16) -> f16 {
    aether_fpx::coshs(x)
}
#[inline(never)]
fn fpx_tanh_f16(x: f16) -> f16 {
    aether_fpx::tanhs(x)
}
#[inline(never)]
fn fpx_exp_f16(x: f16) -> f16 {
    aether_fpx::exps(x)
}
#[inline(never)]
fn fpx_log_f16(x: f16) -> f16 {
    aether_fpx::logs(x)
}
#[inline(never)]
fn fpx_log1p_f16(x: f16) -> f16 {
    aether_fpx::log1ps(x)
}
#[inline(never)]
fn fpx_expm1_f16(x: f16) -> f16 {
    aether_fpx::expm1s(x)
}

// f128
#[inline(never)]
fn fpx_sin_f128(x: f128) -> f128 {
    aether_fpx::trig::sinl(x)
}
#[inline(never)]
fn fpx_cos_f128(x: f128) -> f128 {
    aether_fpx::trig::cosl(x)
}
#[inline(never)]
fn fpx_tan_f128(x: f128) -> f128 {
    aether_fpx::trig::tanl(x)
}
#[inline(never)]
fn fpx_asin_f128(x: f128) -> f128 {
    aether_fpx::trig::asinl(x)
}
#[inline(never)]
fn fpx_acos_f128(x: f128) -> f128 {
    aether_fpx::trig::acosl(x)
}
#[inline(never)]
fn fpx_atan_f128(x: f128) -> f128 {
    aether_fpx::trig::atanl(x)
}
#[inline(never)]
fn fpx_sinh_f128(x: f128) -> f128 {
    aether_fpx::trig::sinhl(x)
}
#[inline(never)]
fn fpx_cosh_f128(x: f128) -> f128 {
    aether_fpx::trig::coshl(x)
}
#[inline(never)]
fn fpx_tanh_f128(x: f128) -> f128 {
    aether_fpx::trig::tanhl(x)
}
#[inline(never)]
fn fpx_exp_f128(x: f128) -> f128 {
    aether_fpx::trig::expl(x)
}
#[inline(never)]
fn fpx_log_f128(x: f128) -> f128 {
    aether_fpx::trig::logl(x)
}
#[inline(never)]
fn fpx_log1p_f128(x: f128) -> f128 {
    aether_fpx::trig::log1pl(x)
}
#[inline(never)]
fn fpx_expm1_f128(x: f128) -> f128 {
    aether_fpx::trig::expm1l(x)
}

/* ---------------- libm wrappers ---------------- */

// f32
#[inline(never)]
fn libm_sin_f32(x: f32) -> f32 {
    libm::sinf(x)
}
#[inline(never)]
fn libm_cos_f32(x: f32) -> f32 {
    libm::cosf(x)
}
#[inline(never)]
fn libm_tan_f32(x: f32) -> f32 {
    libm::tanf(x)
}
#[inline(never)]
fn libm_asin_f32(x: f32) -> f32 {
    libm::asinf(x)
}
#[inline(never)]
fn libm_acos_f32(x: f32) -> f32 {
    libm::acosf(x)
}
#[inline(never)]
fn libm_atan_f32(x: f32) -> f32 {
    libm::atanf(x)
}
#[inline(never)]
fn libm_sinh_f32(x: f32) -> f32 {
    libm::sinhf(x)
}
#[inline(never)]
fn libm_cosh_f32(x: f32) -> f32 {
    libm::coshf(x)
}
#[inline(never)]
fn libm_tanh_f32(x: f32) -> f32 {
    libm::tanhf(x)
}
#[inline(never)]
fn libm_exp_f32(x: f32) -> f32 {
    libm::expf(x)
}
#[inline(never)]
fn libm_log_f32(x: f32) -> f32 {
    libm::logf(x)
}
#[inline(never)]
fn libm_log1p_f32(x: f32) -> f32 {
    libm::log1pf(x)
}
#[inline(never)]
fn libm_expm1_f32(x: f32) -> f32 {
    libm::expm1f(x)
}

// f64
#[inline(never)]
fn libm_sin_f64(x: f64) -> f64 {
    libm::sin(x)
}
#[inline(never)]
fn libm_cos_f64(x: f64) -> f64 {
    libm::cos(x)
}
#[inline(never)]
fn libm_tan_f64(x: f64) -> f64 {
    libm::tan(x)
}
#[inline(never)]
fn libm_asin_f64(x: f64) -> f64 {
    libm::asin(x)
}
#[inline(never)]
fn libm_acos_f64(x: f64) -> f64 {
    libm::acos(x)
}
#[inline(never)]
fn libm_atan_f64(x: f64) -> f64 {
    libm::atan(x)
}
#[inline(never)]
fn libm_sinh_f64(x: f64) -> f64 {
    libm::sinh(x)
}
#[inline(never)]
fn libm_cosh_f64(x: f64) -> f64 {
    libm::cosh(x)
}
#[inline(never)]
fn libm_tanh_f64(x: f64) -> f64 {
    libm::tanh(x)
}
#[inline(never)]
fn libm_exp_f64(x: f64) -> f64 {
    libm::exp(x)
}
#[inline(never)]
fn libm_log_f64(x: f64) -> f64 {
    libm::log(x)
}
#[inline(never)]
fn libm_log1p_f64(x: f64) -> f64 {
    libm::log1p(x)
}
#[inline(never)]
fn libm_expm1_f64(x: f64) -> f64 {
    libm::expm1(x)
}

// f16 via f32 libm
#[inline(never)]
fn libm_sin_f16(x: f16) -> f16 {
    libm::sinf(x as f32) as f16
}
#[inline(never)]
fn libm_cos_f16(x: f16) -> f16 {
    libm::cosf(x as f32) as f16
}
#[inline(never)]
fn libm_tan_f16(x: f16) -> f16 {
    libm::tanf(x as f32) as f16
}
#[inline(never)]
fn libm_asin_f16(x: f16) -> f16 {
    libm::asinf(x as f32) as f16
}
#[inline(never)]
fn libm_acos_f16(x: f16) -> f16 {
    libm::acosf(x as f32) as f16
}
#[inline(never)]
fn libm_atan_f16(x: f16) -> f16 {
    libm::atanf(x as f32) as f16
}
#[inline(never)]
fn libm_sinh_f16(x: f16) -> f16 {
    libm::sinhf(x as f32) as f16
}
#[inline(never)]
fn libm_cosh_f16(x: f16) -> f16 {
    libm::coshf(x as f32) as f16
}
#[inline(never)]
fn libm_tanh_f16(x: f16) -> f16 {
    libm::tanhf(x as f32) as f16
}
#[inline(never)]
fn libm_exp_f16(x: f16) -> f16 {
    libm::expf(x as f32) as f16
}
#[inline(never)]
fn libm_log_f16(x: f16) -> f16 {
    libm::logf(x as f32) as f16
}
#[inline(never)]
fn libm_log1p_f16(x: f16) -> f16 {
    libm::log1pf(x as f32) as f16
}
#[inline(never)]
fn libm_expm1_f16(x: f16) -> f16 {
    libm::expm1f(x as f32) as f16
}

/* ---------------- native wrappers ---------------- */

// f32
#[inline(never)]
fn nat_sin_f32(x: f32) -> f32 {
    x.sin()
}
#[inline(never)]
fn nat_cos_f32(x: f32) -> f32 {
    x.cos()
}
#[inline(never)]
fn nat_tan_f32(x: f32) -> f32 {
    x.tan()
}
#[inline(never)]
fn nat_asin_f32(x: f32) -> f32 {
    x.asin()
}
#[inline(never)]
fn nat_acos_f32(x: f32) -> f32 {
    x.acos()
}
#[inline(never)]
fn nat_atan_f32(x: f32) -> f32 {
    x.atan()
}
#[inline(never)]
fn nat_sinh_f32(x: f32) -> f32 {
    x.sinh()
}
#[inline(never)]
fn nat_cosh_f32(x: f32) -> f32 {
    x.cosh()
}
#[inline(never)]
fn nat_tanh_f32(x: f32) -> f32 {
    x.tanh()
}
#[inline(never)]
fn nat_exp_f32(x: f32) -> f32 {
    x.exp()
}
#[inline(never)]
fn nat_log_f32(x: f32) -> f32 {
    x.ln()
}
#[inline(never)]
fn nat_log1p_f32(x: f32) -> f32 {
    x.ln_1p()
}
#[inline(never)]
fn nat_expm1_f32(x: f32) -> f32 {
    x.exp_m1()
}

// f64
#[inline(never)]
fn nat_sin_f64(x: f64) -> f64 {
    x.sin()
}
#[inline(never)]
fn nat_cos_f64(x: f64) -> f64 {
    x.cos()
}
#[inline(never)]
fn nat_tan_f64(x: f64) -> f64 {
    x.tan()
}
#[inline(never)]
fn nat_asin_f64(x: f64) -> f64 {
    x.asin()
}
#[inline(never)]
fn nat_acos_f64(x: f64) -> f64 {
    x.acos()
}
#[inline(never)]
fn nat_atan_f64(x: f64) -> f64 {
    x.atan()
}
#[inline(never)]
fn nat_sinh_f64(x: f64) -> f64 {
    x.sinh()
}
#[inline(never)]
fn nat_cosh_f64(x: f64) -> f64 {
    x.cosh()
}
#[inline(never)]
fn nat_tanh_f64(x: f64) -> f64 {
    x.tanh()
}
#[inline(never)]
fn nat_exp_f64(x: f64) -> f64 {
    x.exp()
}
#[inline(never)]
fn nat_log_f64(x: f64) -> f64 {
    x.ln()
}
#[inline(never)]
fn nat_log1p_f64(x: f64) -> f64 {
    x.ln_1p()
}
#[inline(never)]
fn nat_expm1_f64(x: f64) -> f64 {
    x.exp_m1()
}

// f16 via f32 native
#[inline(never)]
fn nat_sin_f16(x: f16) -> f16 {
    (x as f32).sin() as f16
}
#[inline(never)]
fn nat_cos_f16(x: f16) -> f16 {
    (x as f32).cos() as f16
}
#[inline(never)]
fn nat_tan_f16(x: f16) -> f16 {
    (x as f32).tan() as f16
}
#[inline(never)]
fn nat_asin_f16(x: f16) -> f16 {
    (x as f32).asin() as f16
}
#[inline(never)]
fn nat_acos_f16(x: f16) -> f16 {
    (x as f32).acos() as f16
}
#[inline(never)]
fn nat_atan_f16(x: f16) -> f16 {
    (x as f32).atan() as f16
}
#[inline(never)]
fn nat_sinh_f16(x: f16) -> f16 {
    (x as f32).sinh() as f16
}
#[inline(never)]
fn nat_cosh_f16(x: f16) -> f16 {
    (x as f32).cosh() as f16
}
#[inline(never)]
fn nat_tanh_f16(x: f16) -> f16 {
    (x as f32).tanh() as f16
}
#[inline(never)]
fn nat_exp_f16(x: f16) -> f16 {
    (x as f32).exp() as f16
}
#[inline(never)]
fn nat_log_f16(x: f16) -> f16 {
    (x as f32).ln() as f16
}
#[inline(never)]
fn nat_log1p_f16(x: f16) -> f16 {
    (x as f32).ln_1p() as f16
}
#[inline(never)]
fn nat_expm1_f16(x: f16) -> f16 {
    (x as f32).exp_m1() as f16
}

// musl does not provide the native f128 transcendental symbols.
#[cfg(not(target_env = "musl"))]
mod native_f128 {
    use super::*;

    #[inline(never)]
    fn sin(x: f128) -> f128 {
        x.sin()
    }
    #[inline(never)]
    fn cos(x: f128) -> f128 {
        x.cos()
    }
    #[inline(never)]
    fn tan(x: f128) -> f128 {
        x.tan()
    }
    #[inline(never)]
    fn asin(x: f128) -> f128 {
        x.asin()
    }
    #[inline(never)]
    fn acos(x: f128) -> f128 {
        x.acos()
    }
    #[inline(never)]
    fn atan(x: f128) -> f128 {
        x.atan()
    }
    #[inline(never)]
    fn sinh(x: f128) -> f128 {
        x.sinh()
    }
    #[inline(never)]
    fn cosh(x: f128) -> f128 {
        x.cosh()
    }
    #[inline(never)]
    fn tanh(x: f128) -> f128 {
        x.tanh()
    }
    #[inline(never)]
    fn exp(x: f128) -> f128 {
        x.exp()
    }
    #[inline(never)]
    fn log(x: f128) -> f128 {
        x.ln()
    }
    #[inline(never)]
    fn log1p(x: f128) -> f128 {
        x.ln_1p()
    }
    #[inline(never)]
    fn expm1(x: f128) -> f128 {
        x.exp_m1()
    }

    pub(super) fn run(kind: BenchKind, op: &str, iters: usize) -> BenchOut {
        let function = match op {
            "sin" => sin,
            "cos" => cos,
            "tan" => tan,
            "asin" => asin,
            "acos" => acos,
            "atan" => atan,
            "sinh" => sinh,
            "cosh" => cosh,
            "tanh" => tanh,
            "exp" => exp,
            "log" => log,
            "log1p" => log1p,
            "expm1" => expm1,
            _ => usage(),
        };
        run_f128(kind, op, iters, function)
    }
}

/* ---------------- x87 wrappers ---------------- */

#[cfg(fpx_x87)]
type X87Op = unsafe extern "C" fn(f64) -> f64;

#[cfg(fpx_x87)]
unsafe extern "C" {
    fn fpx_sin_x87(x: f64) -> f64;
    fn fpx_cos_x87(x: f64) -> f64;
    fn fpx_tan_x87(x: f64) -> f64;
    fn fpx_sqrt_x87(x: f64) -> f64;
}

#[derive(Clone, Copy)]
struct RangeSpecF16 {
    lo: f16,
    hi: f16,
    step: f16,
}

#[derive(Clone, Copy)]
struct RangeSpecF32 {
    lo: f32,
    hi: f32,
    step: f32,
}

#[derive(Clone, Copy)]
struct RangeSpecF64 {
    lo: f64,
    hi: f64,
    step: f64,
}

#[derive(Clone, Copy)]
struct RangeSpecF128 {
    lo: f128,
    hi: f128,
    step: f128,
}

#[inline(always)]
fn range_spec_f16(lo: f32, hi: f32, points: f32) -> RangeSpecF16 {
    RangeSpecF16 {
        lo: lo as f16,
        hi: hi as f16,
        step: ((hi - lo) / points) as f16,
    }
}

#[inline(always)]
fn range_spec_f32(lo: f32, hi: f32, points: f32) -> RangeSpecF32 {
    RangeSpecF32 {
        lo,
        hi,
        step: (hi - lo) / points,
    }
}

#[inline(always)]
fn range_spec_f64(lo: f64, hi: f64, points: f64) -> RangeSpecF64 {
    RangeSpecF64 {
        lo,
        hi,
        step: (hi - lo) / points,
    }
}

#[inline(always)]
fn range_spec_f128(lo: f64, hi: f64, points: f64) -> RangeSpecF128 {
    let lo = lo as f128;
    let hi = hi as f128;
    RangeSpecF128 {
        lo,
        hi,
        step: (hi - lo) / points as f128,
    }
}

#[inline(always)]
fn public_range_f16(op: &str) -> RangeSpecF16 {
    match op {
        "sin" | "cos" => range_spec_f16(-8.0, 8.0, 1024.0),
        "tan" => range_spec_f16(-1.25, 1.25, 1024.0),
        "asin" | "acos" => range_spec_f16(-1.0, 1.0, 1024.0),
        "atan" => range_spec_f16(-8.0, 8.0, 1024.0),
        "sinh" | "cosh" | "tanh" | "exp" | "expm1" => range_spec_f16(-4.0, 4.0, 1024.0),
        "log" => range_spec_f16(0.25, 16.0, 1024.0),
        "log1p" => range_spec_f16(-0.875, 8.0, 1024.0),
        _ => range_spec_f16(0.5, 0.75, 256.0),
    }
}

#[inline(always)]
fn public_range_f32(op: &str) -> RangeSpecF32 {
    match op {
        "sin" | "cos" => range_spec_f32(-8.0, 8.0, 65536.0),
        "tan" => range_spec_f32(-1.25, 1.25, 65536.0),
        "asin" | "acos" => range_spec_f32(-1.0, 1.0, 65536.0),
        "atan" => range_spec_f32(-8.0, 8.0, 65536.0),
        "sinh" | "cosh" | "tanh" | "exp" | "expm1" => range_spec_f32(-4.0, 4.0, 65536.0),
        "log" => range_spec_f32(0.25, 16.0, 65536.0),
        "log1p" => range_spec_f32(-0.875, 8.0, 65536.0),
        _ => range_spec_f32(0.5, 0.75, 4096.0),
    }
}

#[inline(always)]
fn public_range_f64(op: &str) -> RangeSpecF64 {
    match op {
        "sin" | "cos" => range_spec_f64(-8.0, 8.0, 65536.0),
        "tan" => range_spec_f64(-1.25, 1.25, 65536.0),
        "asin" | "acos" => range_spec_f64(-1.0, 1.0, 65536.0),
        "atan" => range_spec_f64(-8.0, 8.0, 65536.0),
        "sinh" | "cosh" | "tanh" | "exp" | "expm1" => range_spec_f64(-4.0, 4.0, 65536.0),
        "log" => range_spec_f64(0.25, 16.0, 65536.0),
        "log1p" => range_spec_f64(-0.875, 8.0, 65536.0),
        _ => range_spec_f64(0.5, 0.75, 4096.0),
    }
}

#[inline(always)]
fn public_range_f128(op: &str) -> RangeSpecF128 {
    match op {
        "sin" | "cos" => range_spec_f128(-8.0, 8.0, 65536.0),
        "tan" => range_spec_f128(-1.25, 1.25, 65536.0),
        "asin" | "acos" => range_spec_f128(-1.0, 1.0, 65536.0),
        "atan" => range_spec_f128(-8.0, 8.0, 65536.0),
        "sinh" | "cosh" | "tanh" | "exp" | "expm1" => range_spec_f128(-4.0, 4.0, 65536.0),
        "log" => range_spec_f128(0.25, 16.0, 65536.0),
        "log1p" => range_spec_f128(-0.875, 8.0, 65536.0),
        _ => range_spec_f128(0.5, 0.75, 4096.0),
    }
}

#[inline(always)]
fn kernel_range_f32(op: &str) -> RangeSpecF32 {
    match op {
        "sin" => range_spec_f32(
            aether_fpx::trig::sin::SIN_INTERVAL_MIN_F32,
            aether_fpx::trig::sin::SIN_INTERVAL_MAX_F32,
            65536.0,
        ),
        "cos" => range_spec_f32(
            aether_fpx::trig::cos::COS_INTERVAL_MIN_F32,
            aether_fpx::trig::cos::COS_INTERVAL_MAX_F32,
            65536.0,
        ),
        "tan" => range_spec_f32(
            aether_fpx::trig::tan::TAN_INTERVAL_MIN_F32,
            aether_fpx::trig::tan::TAN_INTERVAL_MAX_F32,
            65536.0,
        ),
        "asin" => range_spec_f32(
            aether_fpx::trig::asin::ASIN_INTERVAL_MIN_F32,
            aether_fpx::trig::asin::ASIN_INTERVAL_MAX_F32,
            65536.0,
        ),
        "acos" => range_spec_f32(
            aether_fpx::trig::acos::ACOS_INTERVAL_MIN_F32,
            aether_fpx::trig::acos::ACOS_INTERVAL_MAX_F32,
            65536.0,
        ),
        "atan" => range_spec_f32(
            aether_fpx::trig::atan::ATAN_INTERVAL_MIN_F32,
            aether_fpx::trig::atan::ATAN_INTERVAL_MAX_F32,
            65536.0,
        ),
        "sinh" | "cosh" => range_spec_f32(0.0, 0.5, 65536.0),
        "tanh" => range_spec_f32(0.0, 0.25, 65536.0),
        "exp" => range_spec_f32(
            aether_fpx::trig::exp::EXP_INTERVAL_MIN_F32,
            aether_fpx::trig::exp::EXP_INTERVAL_MAX_F32,
            65536.0,
        ),
        "expm1" => range_spec_f32(
            aether_fpx::trig::expm1::EXPM1_INTERVAL_MIN_F32,
            aether_fpx::trig::expm1::EXPM1_INTERVAL_MAX_F32,
            65536.0,
        ),
        "log" => range_spec_f32(
            aether_fpx::trig::log::LOG_INTERVAL_MIN_F32,
            aether_fpx::trig::log::LOG_INTERVAL_MAX_F32,
            65536.0,
        ),
        "log1p" => range_spec_f32(
            aether_fpx::trig::log1p::LOG1P_INTERVAL_MIN_F32,
            aether_fpx::trig::log1p::LOG1P_INTERVAL_MAX_F32,
            65536.0,
        ),
        _ => range_spec_f32(-0.7853982, 0.7853982, 65536.0),
    }
}

#[inline(always)]
fn kernel_range_f64(op: &str) -> RangeSpecF64 {
    match op {
        "sin" => range_spec_f64(
            aether_fpx::trig::sin::SIN_INTERVAL_MIN_F64,
            aether_fpx::trig::sin::SIN_INTERVAL_MAX_F64,
            65536.0,
        ),
        "cos" => range_spec_f64(
            aether_fpx::trig::cos::COS_INTERVAL_MIN_F64,
            aether_fpx::trig::cos::COS_INTERVAL_MAX_F64,
            65536.0,
        ),
        "tan" => range_spec_f64(
            aether_fpx::trig::tan::TAN_INTERVAL_MIN_F64,
            aether_fpx::trig::tan::TAN_INTERVAL_MAX_F64,
            65536.0,
        ),
        "asin" => range_spec_f64(
            aether_fpx::trig::asin::ASIN_INTERVAL_MIN_F64,
            aether_fpx::trig::asin::ASIN_INTERVAL_MAX_F64,
            65536.0,
        ),
        "acos" => range_spec_f64(
            aether_fpx::trig::acos::ACOS_INTERVAL_MIN_F64,
            aether_fpx::trig::acos::ACOS_INTERVAL_MAX_F64,
            65536.0,
        ),
        "atan" => range_spec_f64(
            aether_fpx::trig::atan::ATAN_INTERVAL_MIN_F64,
            aether_fpx::trig::atan::ATAN_INTERVAL_MAX_F64,
            65536.0,
        ),
        "sinh" | "cosh" => range_spec_f64(0.0, 0.125, 65536.0),
        "tanh" => range_spec_f64(0.0, 0.0625, 65536.0),
        "exp" => range_spec_f64(
            aether_fpx::trig::exp::EXP_INTERVAL_MIN_F64,
            aether_fpx::trig::exp::EXP_INTERVAL_MAX_F64,
            65536.0,
        ),
        "expm1" => range_spec_f64(
            aether_fpx::trig::expm1::EXPM1_INTERVAL_MIN_F64,
            aether_fpx::trig::expm1::EXPM1_INTERVAL_MAX_F64,
            65536.0,
        ),
        "log" => range_spec_f64(
            aether_fpx::trig::log::LOG_INTERVAL_MIN_F64,
            aether_fpx::trig::log::LOG_INTERVAL_MAX_F64,
            65536.0,
        ),
        "log1p" => range_spec_f64(
            aether_fpx::trig::log1p::LOG1P_INTERVAL_MIN_F64,
            aether_fpx::trig::log1p::LOG1P_INTERVAL_MAX_F64,
            65536.0,
        ),
        _ => range_spec_f64(-0.7853981633974483, 0.7853981633974483, 65536.0),
    }
}

/* ---------------- bench loops ---------------- */

#[inline(always)]
fn wrap_f32(x: f32, lo: f32, hi: f32) -> f32 {
    if x > hi {
        lo
    } else {
        x
    }
}

#[inline(always)]
fn wrap_f64(x: f64, lo: f64, hi: f64) -> f64 {
    if x > hi {
        lo
    } else {
        x
    }
}

#[inline(always)]
fn wrap_f16(x: f16, lo: f16, hi: f16) -> f16 {
    if x > hi {
        lo
    } else {
        x
    }
}

#[inline(always)]
fn wrap_f128(x: f128, lo: f128, hi: f128) -> f128 {
    if x > hi {
        lo
    } else {
        x
    }
}

#[inline(never)]
fn bench_f32_fixed(x: f32, iters: usize, f: fn(f32) -> f32) -> f32 {
    let x = black_box(x);
    let mut acc = 0.0f32;

    for _ in 0..iters {
        acc += black_box(f(black_box(x)));
    }

    black_box(acc)
}

#[inline(never)]
fn bench_f64_fixed(x: f64, iters: usize, f: fn(f64) -> f64) -> f64 {
    let x = black_box(x);
    let mut acc = 0.0f64;

    for _ in 0..iters {
        acc += black_box(f(black_box(x)));
    }

    black_box(acc)
}

#[inline(never)]
fn bench_f16_fixed(x: f16, iters: usize, f: fn(f16) -> f16) -> f16 {
    let x = black_box(x);
    let mut acc = f16::from_bits(0);

    for _ in 0..iters {
        acc += black_box(f(black_box(x)));
    }

    black_box(acc)
}

#[inline(never)]
fn bench_f128_fixed(x: f128, iters: usize, f: fn(f128) -> f128) -> f128 {
    let x = black_box(x);
    let mut acc = 0.0f128;

    for _ in 0..iters {
        acc += black_box(f(black_box(x)));
    }

    black_box(acc)
}

#[inline(never)]
#[cfg(fpx_x87)]
fn bench_x87_fixed(x: f64, iters: usize, f: X87Op) -> f64 {
    let x = black_box(x);
    let mut acc = 0.0f64;

    for _ in 0..iters {
        acc += black_box(unsafe { f(black_box(x)) });
    }

    black_box(acc)
}

#[inline(never)]
fn bench_f32_bounded(mut x: f32, iters: usize, f: fn(f32) -> f32, spec: RangeSpecF32) -> f32 {
    let mut acc = 0.0f32;

    for _ in 0..iters {
        acc += black_box(f(black_box(x)));
        x = wrap_f32(x + spec.step, spec.lo, spec.hi);
    }

    black_box(acc)
}

#[inline(never)]
fn bench_f64_bounded(mut x: f64, iters: usize, f: fn(f64) -> f64, spec: RangeSpecF64) -> f64 {
    let mut acc = 0.0f64;

    for _ in 0..iters {
        acc += black_box(f(black_box(x)));
        x = wrap_f64(x + spec.step, spec.lo, spec.hi);
    }

    black_box(acc)
}

#[inline(never)]
fn bench_f16_bounded(mut x: f16, iters: usize, f: fn(f16) -> f16, spec: RangeSpecF16) -> f16 {
    let mut acc = f16::from_bits(0);

    for _ in 0..iters {
        acc += black_box(f(black_box(x)));
        x = wrap_f16(x + spec.step, spec.lo, spec.hi);
    }

    black_box(acc)
}

#[inline(never)]
fn bench_f128_bounded(mut x: f128, iters: usize, f: fn(f128) -> f128, spec: RangeSpecF128) -> f128 {
    let mut acc = 0.0f128;

    for _ in 0..iters {
        acc += black_box(f(black_box(x)));
        x = wrap_f128(x + spec.step, spec.lo, spec.hi);
    }

    black_box(acc)
}

#[inline(never)]
#[cfg(fpx_x87)]
fn bench_x87_bounded(mut x: f64, iters: usize, f: X87Op, spec: RangeSpecF64) -> f64 {
    let mut acc = 0.0f64;

    for _ in 0..iters {
        acc += black_box(unsafe { f(black_box(x)) });
        x = wrap_f64(x + spec.step, spec.lo, spec.hi);
    }

    black_box(acc)
}

/* ---------------- output sink ---------------- */

#[derive(Clone, Copy)]
enum BenchOut {
    F16(u16),
    F32(u32),
    F64(u64),
    F128(u128),
}

#[inline(never)]
fn consume(out: BenchOut) {
    match black_box(out) {
        BenchOut::F16(bits) => {
            eprintln!("sink=f16:{bits:04x}");
        }
        BenchOut::F32(bits) => {
            eprintln!("sink=f32:{bits:08x}");
        }
        BenchOut::F64(bits) => {
            eprintln!("sink=f64:{bits:016x}");
        }
        BenchOut::F128(bits) => {
            eprintln!("sink=f128:{bits:032x}");
        }
    }
}

/* ---------------- dispatch helpers ---------------- */

fn fpx_kernel_horner_f32(op: &str) -> Option<fn(f32) -> f32> {
    Some(match op {
        "sin" => fpx_sin_kernel_horner_f32,
        "cos" => fpx_cos_kernel_horner_f32,
        "tan" => fpx_tan_kernel_horner_f32,
        "asin" => fpx_asin_kernel_horner_f32,
        "acos" => fpx_acos_kernel_horner_f32,
        "atan" => fpx_atan_kernel_horner_f32,
        "sinh" => fpx_sinh_kernel_horner_f32,
        "cosh" => fpx_cosh_kernel_horner_f32,
        "tanh" => fpx_tanh_kernel_horner_f32,
        "exp" => fpx_exp_kernel_horner_f32,
        "expm1" => fpx_expm1_kernel_horner_f32,
        "log" => fpx_log_kernel_horner_f32,
        "log1p" => fpx_log1p_kernel_horner_f32,
        _ => return None,
    })
}

fn fpx_kernel_estrin_f32(op: &str) -> Option<fn(f32) -> f32> {
    Some(match op {
        "sin" => fpx_sin_kernel_estrin_f32,
        "cos" => fpx_cos_kernel_estrin_f32,
        "tan" => fpx_tan_kernel_estrin_f32,
        "asin" => fpx_asin_kernel_estrin_f32,
        "acos" => fpx_acos_kernel_estrin_f32,
        "atan" => fpx_atan_kernel_estrin_f32,
        "sinh" => fpx_sinh_kernel_estrin_f32,
        "cosh" => fpx_cosh_kernel_estrin_f32,
        "tanh" => fpx_tanh_kernel_estrin_f32,
        "exp" => fpx_exp_kernel_estrin_f32,
        "expm1" => fpx_expm1_kernel_estrin_f32,
        "log" => fpx_log_kernel_estrin_f32,
        "log1p" => fpx_log1p_kernel_estrin_f32,
        _ => return None,
    })
}

fn fpx_kernel_horner_f64(op: &str) -> Option<fn(f64) -> f64> {
    Some(match op {
        "sin" => fpx_sin_kernel_horner_f64,
        "cos" => fpx_cos_kernel_horner_f64,
        "tan" => fpx_tan_kernel_horner_f64,
        "asin" => fpx_asin_kernel_horner_f64,
        "acos" => fpx_acos_kernel_horner_f64,
        "atan" => fpx_atan_kernel_horner_f64,
        "sinh" => fpx_sinh_kernel_horner_f64,
        "cosh" => fpx_cosh_kernel_horner_f64,
        "tanh" => fpx_tanh_kernel_horner_f64,
        "exp" => fpx_exp_kernel_horner_f64,
        "expm1" => fpx_expm1_kernel_horner_f64,
        "log" => fpx_log_kernel_horner_f64,
        "log1p" => fpx_log1p_kernel_horner_f64,
        _ => return None,
    })
}

fn fpx_kernel_estrin_f64(op: &str) -> Option<fn(f64) -> f64> {
    Some(match op {
        "sin" => fpx_sin_kernel_estrin_f64,
        "cos" => fpx_cos_kernel_estrin_f64,
        "tan" => fpx_tan_kernel_estrin_f64,
        "asin" => fpx_asin_kernel_estrin_f64,
        "acos" => fpx_acos_kernel_estrin_f64,
        "atan" => fpx_atan_kernel_estrin_f64,
        "sinh" => fpx_sinh_kernel_estrin_f64,
        "cosh" => fpx_cosh_kernel_estrin_f64,
        "tanh" => fpx_tanh_kernel_estrin_f64,
        "exp" => fpx_exp_kernel_estrin_f64,
        "expm1" => fpx_expm1_kernel_estrin_f64,
        "log" => fpx_log_kernel_estrin_f64,
        "log1p" => fpx_log1p_kernel_estrin_f64,
        _ => return None,
    })
}

#[inline(never)]
fn run_f32(kind: BenchKind, op: &str, iters: usize, f: fn(f32) -> f32) -> BenchOut {
    let spec = public_range_f32(op);
    let fixed = spec.lo + (spec.hi - spec.lo) * 0.5;
    let y = match kind {
        BenchKind::Fixed => bench_f32_fixed(fixed, iters, f),
        BenchKind::Bounded => bench_f32_bounded(spec.lo, iters, f, spec),
    };
    BenchOut::F32(y.to_bits())
}

#[inline(never)]
fn run_f64(kind: BenchKind, op: &str, iters: usize, f: fn(f64) -> f64) -> BenchOut {
    let spec = public_range_f64(op);
    let fixed = spec.lo + (spec.hi - spec.lo) * 0.5;
    let y = match kind {
        BenchKind::Fixed => bench_f64_fixed(fixed, iters, f),
        BenchKind::Bounded => bench_f64_bounded(spec.lo, iters, f, spec),
    };
    BenchOut::F64(y.to_bits())
}

#[inline(never)]
fn run_f16(kind: BenchKind, op: &str, iters: usize, f: fn(f16) -> f16) -> BenchOut {
    let spec = public_range_f16(op);
    let fixed = spec.lo + (spec.hi - spec.lo) * (0.5f32 as f16);
    let y = match kind {
        BenchKind::Fixed => bench_f16_fixed(fixed, iters, f),
        BenchKind::Bounded => bench_f16_bounded(spec.lo, iters, f, spec),
    };
    BenchOut::F16(y.to_bits())
}

#[inline(never)]
fn run_f128(kind: BenchKind, op: &str, iters: usize, f: fn(f128) -> f128) -> BenchOut {
    let spec = public_range_f128(op);
    let fixed = spec.lo + (spec.hi - spec.lo) * 0.5f128;
    let y = match kind {
        BenchKind::Fixed => bench_f128_fixed(fixed, iters, f),
        BenchKind::Bounded => bench_f128_bounded(spec.lo, iters, f, spec),
    };
    BenchOut::F128(y.to_bits())
}

#[inline(never)]
fn run_f32_kernel(kind: BenchKind, op: &str, iters: usize, f: fn(f32) -> f32) -> BenchOut {
    let spec = kernel_range_f32(op);
    let fixed = spec.lo + (spec.hi - spec.lo) * 0.5;
    let y = match kind {
        BenchKind::Fixed => bench_f32_fixed(fixed, iters, f),
        BenchKind::Bounded => bench_f32_bounded(spec.lo, iters, f, spec),
    };
    BenchOut::F32(y.to_bits())
}

#[inline(never)]
fn run_f64_kernel(kind: BenchKind, op: &str, iters: usize, f: fn(f64) -> f64) -> BenchOut {
    let spec = kernel_range_f64(op);
    let fixed = spec.lo + (spec.hi - spec.lo) * 0.5;
    let y = match kind {
        BenchKind::Fixed => bench_f64_fixed(fixed, iters, f),
        BenchKind::Bounded => bench_f64_bounded(spec.lo, iters, f, spec),
    };
    BenchOut::F64(y.to_bits())
}

#[cfg(feature = "simd")]
#[inline(never)]
fn run_f32_arithmetic_path(
    kind: BenchKind,
    op: &str,
    iters: usize,
    path: aether_fpx::bench_paths::ArithmeticPath,
) -> BenchOut {
    let spec = kernel_range_f32(op);
    let fixed = spec.lo + (spec.hi - spec.lo) * 0.5;
    let mut x = match kind {
        BenchKind::Fixed => fixed,
        BenchKind::Bounded => spec.lo,
    };
    let mut values = Vec::with_capacity(iters);
    for _ in 0..iters {
        values.push(black_box(x));
        if matches!(kind, BenchKind::Bounded) {
            x = wrap_f32(x + spec.step, spec.lo, spec.hi);
        }
    }
    let mut workspace = aether_fpx::bench_paths::WorkspaceF32::with_capacity(iters);
    aether_fpx::bench_paths::eval_f32_with_workspace(
        op,
        path,
        black_box(&mut values),
        black_box(&mut workspace),
    );
    let sink = values
        .iter()
        .fold(0_u32, |bits, value| bits ^ black_box(value.to_bits()));
    BenchOut::F32(sink)
}

#[cfg(feature = "simd")]
#[inline(never)]
fn run_f64_arithmetic_path(
    kind: BenchKind,
    op: &str,
    iters: usize,
    path: aether_fpx::bench_paths::ArithmeticPath,
) -> BenchOut {
    let spec = kernel_range_f64(op);
    let fixed = spec.lo + (spec.hi - spec.lo) * 0.5;
    let mut x = match kind {
        BenchKind::Fixed => fixed,
        BenchKind::Bounded => spec.lo,
    };
    let mut values = Vec::with_capacity(iters);
    for _ in 0..iters {
        values.push(black_box(x));
        if matches!(kind, BenchKind::Bounded) {
            x = wrap_f64(x + spec.step, spec.lo, spec.hi);
        }
    }
    let mut workspace = aether_fpx::bench_paths::WorkspaceF64::with_capacity(iters);
    aether_fpx::bench_paths::eval_f64_with_workspace(
        op,
        path,
        black_box(&mut values),
        black_box(&mut workspace),
    );
    let sink = values
        .iter()
        .fold(0_u64, |bits, value| bits ^ black_box(value.to_bits()));
    BenchOut::F64(sink)
}

#[cfg(feature = "simd")]
#[inline(never)]
fn run_f32_public_simd(
    kind: BenchKind,
    op: &str,
    iters: usize,
    method: aether_fpx::public_simd::HorizontalMethod,
) -> BenchOut {
    let spec = public_range_f32(op);
    let fixed = spec.lo + (spec.hi - spec.lo) * 0.5;
    let mut x = if matches!(kind, BenchKind::Fixed) {
        fixed
    } else {
        spec.lo
    };
    let mut acc = 0.0_f32;
    for _ in 0..iters {
        acc += black_box(aether_fpx::public_simd::eval_f32(op, method, black_box(x)));
        if matches!(kind, BenchKind::Bounded) {
            x = wrap_f32(x + spec.step, spec.lo, spec.hi);
        }
    }
    BenchOut::F32(black_box(acc).to_bits())
}

#[cfg(feature = "simd")]
#[inline(never)]
fn run_f64_public_simd(
    kind: BenchKind,
    op: &str,
    iters: usize,
    method: aether_fpx::public_simd::HorizontalMethod,
) -> BenchOut {
    let spec = public_range_f64(op);
    let fixed = spec.lo + (spec.hi - spec.lo) * 0.5;
    let mut x = if matches!(kind, BenchKind::Fixed) {
        fixed
    } else {
        spec.lo
    };
    let mut acc = 0.0_f64;
    for _ in 0..iters {
        acc += black_box(aether_fpx::public_simd::eval_f64(op, method, black_box(x)));
        if matches!(kind, BenchKind::Bounded) {
            x = wrap_f64(x + spec.step, spec.lo, spec.hi);
        }
    }
    BenchOut::F64(black_box(acc).to_bits())
}

#[cfg(feature = "simd")]
#[inline(never)]
fn run_f32_raw_polynomial(
    kind: BenchKind,
    op: &str,
    iters: usize,
    path: aether_fpx::bench_paths::ArithmeticPath,
) -> BenchOut {
    let mut x = if matches!(kind, BenchKind::Fixed) {
        0.25
    } else {
        0.0
    };
    let mut values = Vec::with_capacity(iters);
    for _ in 0..iters {
        values.push(black_box(x));
        if matches!(kind, BenchKind::Bounded) {
            x += 1.0 / 4096.0;
            if x > 0.5 {
                x = 0.0;
            }
        }
    }
    aether_fpx::bench_paths::eval_polynomial_f32(op, 0, path, black_box(&mut values));
    let sink = values
        .iter()
        .fold(0_u32, |bits, value| bits ^ black_box(value.to_bits()));
    BenchOut::F32(sink)
}

#[cfg(feature = "simd")]
#[inline(never)]
fn run_f64_raw_polynomial(
    kind: BenchKind,
    op: &str,
    iters: usize,
    path: aether_fpx::bench_paths::ArithmeticPath,
) -> BenchOut {
    let mut x = if matches!(kind, BenchKind::Fixed) {
        0.25
    } else {
        0.0
    };
    let mut values = Vec::with_capacity(iters);
    for _ in 0..iters {
        values.push(black_box(x));
        if matches!(kind, BenchKind::Bounded) {
            x += 1.0 / 4096.0;
            if x > 0.5 {
                x = 0.0;
            }
        }
    }
    aether_fpx::bench_paths::eval_polynomial_f64(op, 0, path, black_box(&mut values));
    let sink = values
        .iter()
        .fold(0_u64, |bits, value| bits ^ black_box(value.to_bits()));
    BenchOut::F64(sink)
}

#[inline(never)]
#[cfg(fpx_x87)]
fn run_x87(kind: BenchKind, op: &str, iters: usize, f: X87Op) -> BenchOut {
    let spec = public_range_f64(op);
    let fixed = spec.lo + (spec.hi - spec.lo) * 0.5;
    let y = match kind {
        BenchKind::Fixed => bench_x87_fixed(fixed, iters, f),
        BenchKind::Bounded => bench_x87_bounded(spec.lo, iters, f, spec),
    };
    BenchOut::F64(y.to_bits())
}

fn usage() -> ! {
    eprintln!(
        "usage: trig_probe <sin|cos|tan|asin|acos|atan|sinh|cosh|tanh|exp|log|log1p|expm1|sqrt> \
         <fpx|fpx_blocked_horner|fpx_horizontal_estrin|fpx_kernel_horner|fpx_kernel_estrin|fpx_kernel_scalar|fpx_kernel_fma|fpx_kernel_simd|fpx_kernel_simd_estrin|fpx_kernel_raw_scalar|fpx_kernel_raw_fma|fpx_kernel_raw_simd|fpx_kernel_raw_simd_estrin|libm|nat|x87> \
         <f16|f32|f64|f128> [iters] [fixed|bounded]"
    );
    std::process::exit(2);
}

/* ---------------- main ---------------- */

fn main() {
    let mut args = std::env::args().skip(1);

    let op = args.next().unwrap_or_else(|| usage());
    let mode = args.next().unwrap_or_else(|| usage());
    let prec = args.next().unwrap_or_else(|| usage());

    let iters: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(100_000);

    let kind = match args.next().as_deref() {
        Some("fixed") => BenchKind::Fixed,
        Some("bounded") => BenchKind::Bounded,
        None => DEFAULT_BENCH_KIND,
        Some(_) => usage(),
    };

    let op_name = op.as_str();

    let out = match (op_name, mode.as_str(), prec.as_str()) {
        #[cfg(feature = "simd")]
        (_, "fpx_blocked_horner", "f32") => run_f32_public_simd(
            kind,
            op_name,
            iters,
            aether_fpx::public_simd::HorizontalMethod::BlockedHorner,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_blocked_horner", "f64") => run_f64_public_simd(
            kind,
            op_name,
            iters,
            aether_fpx::public_simd::HorizontalMethod::BlockedHorner,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_horizontal_estrin", "f32") => run_f32_public_simd(
            kind,
            op_name,
            iters,
            aether_fpx::public_simd::HorizontalMethod::Estrin,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_horizontal_estrin", "f64") => run_f64_public_simd(
            kind,
            op_name,
            iters,
            aether_fpx::public_simd::HorizontalMethod::Estrin,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_raw_scalar", "f32") => run_f32_raw_polynomial(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Scalar,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_raw_fma", "f32") => run_f32_raw_polynomial(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Fma,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_raw_simd", "f32") => run_f32_raw_polynomial(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Simd,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_raw_simd_estrin", "f32") => run_f32_raw_polynomial(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::SimdEstrin,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_raw_scalar", "f64") => run_f64_raw_polynomial(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Scalar,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_raw_fma", "f64") => run_f64_raw_polynomial(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Fma,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_raw_simd", "f64") => run_f64_raw_polynomial(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Simd,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_raw_simd_estrin", "f64") => run_f64_raw_polynomial(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::SimdEstrin,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_scalar", "f32") => run_f32_arithmetic_path(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Scalar,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_fma", "f32") => run_f32_arithmetic_path(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Fma,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_simd", "f32") => run_f32_arithmetic_path(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Simd,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_simd_estrin", "f32") => run_f32_arithmetic_path(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::SimdEstrin,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_scalar", "f64") => run_f64_arithmetic_path(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Scalar,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_fma", "f64") => run_f64_arithmetic_path(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Fma,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_simd", "f64") => run_f64_arithmetic_path(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::Simd,
        ),
        #[cfg(feature = "simd")]
        (_, "fpx_kernel_simd_estrin", "f64") => run_f64_arithmetic_path(
            kind,
            op_name,
            iters,
            aether_fpx::bench_paths::ArithmeticPath::SimdEstrin,
        ),

        // ---------- FPX f32 ----------
        ("sin", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_sin_f32),
        ("cos", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_cos_f32),
        ("tan", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_tan_f32),
        ("asin", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_asin_f32),
        ("acos", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_acos_f32),
        ("atan", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_atan_f32),
        ("sinh", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_sinh_f32),
        ("cosh", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_cosh_f32),
        ("tanh", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_tanh_f32),
        ("exp", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_exp_f32),
        ("log", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_log_f32),
        ("log1p", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_log1p_f32),
        ("expm1", "fpx", "f32") => run_f32(kind, op_name, iters, fpx_expm1_f32),

        (_, "fpx_kernel_horner", "f32") => {
            let f = fpx_kernel_horner_f32(op_name).unwrap_or_else(|| usage());
            run_f32_kernel(kind, op_name, iters, f)
        }
        (_, "fpx_kernel_estrin", "f32") => {
            let f = fpx_kernel_estrin_f32(op_name).unwrap_or_else(|| usage());
            run_f32_kernel(kind, op_name, iters, f)
        }

        // ---------- FPX f64 ----------
        ("sin", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_sin_f64),
        ("cos", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_cos_f64),
        ("tan", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_tan_f64),
        ("asin", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_asin_f64),
        ("acos", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_acos_f64),
        ("atan", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_atan_f64),
        ("sinh", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_sinh_f64),
        ("cosh", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_cosh_f64),
        ("tanh", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_tanh_f64),
        ("exp", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_exp_f64),
        ("log", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_log_f64),
        ("log1p", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_log1p_f64),
        ("expm1", "fpx", "f64") => run_f64(kind, op_name, iters, fpx_expm1_f64),

        (_, "fpx_kernel_horner", "f64") => {
            let f = fpx_kernel_horner_f64(op_name).unwrap_or_else(|| usage());
            run_f64_kernel(kind, op_name, iters, f)
        }
        (_, "fpx_kernel_estrin", "f64") => {
            let f = fpx_kernel_estrin_f64(op_name).unwrap_or_else(|| usage());
            run_f64_kernel(kind, op_name, iters, f)
        }

        // ---------- FPX f16 ----------
        ("sin", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_sin_f16),
        ("cos", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_cos_f16),
        ("tan", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_tan_f16),
        ("asin", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_asin_f16),
        ("acos", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_acos_f16),
        ("atan", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_atan_f16),
        ("sinh", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_sinh_f16),
        ("cosh", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_cosh_f16),
        ("tanh", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_tanh_f16),
        ("exp", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_exp_f16),
        ("log", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_log_f16),
        ("log1p", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_log1p_f16),
        ("expm1", "fpx", "f16") => run_f16(kind, op_name, iters, fpx_expm1_f16),

        // ---------- FPX f128 ----------
        ("sin", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_sin_f128),
        ("cos", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_cos_f128),
        ("tan", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_tan_f128),
        ("asin", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_asin_f128),
        ("acos", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_acos_f128),
        ("atan", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_atan_f128),
        ("sinh", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_sinh_f128),
        ("cosh", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_cosh_f128),
        ("tanh", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_tanh_f128),
        ("exp", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_exp_f128),
        ("log", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_log_f128),
        ("log1p", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_log1p_f128),
        ("expm1", "fpx", "f128") => run_f128(kind, op_name, iters, fpx_expm1_f128),

        // ---------- libm f32 ----------
        ("sin", "libm", "f32") => run_f32(kind, op_name, iters, libm_sin_f32),
        ("cos", "libm", "f32") => run_f32(kind, op_name, iters, libm_cos_f32),
        ("tan", "libm", "f32") => run_f32(kind, op_name, iters, libm_tan_f32),
        ("asin", "libm", "f32") => run_f32(kind, op_name, iters, libm_asin_f32),
        ("acos", "libm", "f32") => run_f32(kind, op_name, iters, libm_acos_f32),
        ("atan", "libm", "f32") => run_f32(kind, op_name, iters, libm_atan_f32),
        ("sinh", "libm", "f32") => run_f32(kind, op_name, iters, libm_sinh_f32),
        ("cosh", "libm", "f32") => run_f32(kind, op_name, iters, libm_cosh_f32),
        ("tanh", "libm", "f32") => run_f32(kind, op_name, iters, libm_tanh_f32),
        ("exp", "libm", "f32") => run_f32(kind, op_name, iters, libm_exp_f32),
        ("log", "libm", "f32") => run_f32(kind, op_name, iters, libm_log_f32),
        ("log1p", "libm", "f32") => run_f32(kind, op_name, iters, libm_log1p_f32),
        ("expm1", "libm", "f32") => run_f32(kind, op_name, iters, libm_expm1_f32),

        // ---------- libm f64 ----------
        ("sin", "libm", "f64") => run_f64(kind, op_name, iters, libm_sin_f64),
        ("cos", "libm", "f64") => run_f64(kind, op_name, iters, libm_cos_f64),
        ("tan", "libm", "f64") => run_f64(kind, op_name, iters, libm_tan_f64),
        ("asin", "libm", "f64") => run_f64(kind, op_name, iters, libm_asin_f64),
        ("acos", "libm", "f64") => run_f64(kind, op_name, iters, libm_acos_f64),
        ("atan", "libm", "f64") => run_f64(kind, op_name, iters, libm_atan_f64),
        ("sinh", "libm", "f64") => run_f64(kind, op_name, iters, libm_sinh_f64),
        ("cosh", "libm", "f64") => run_f64(kind, op_name, iters, libm_cosh_f64),
        ("tanh", "libm", "f64") => run_f64(kind, op_name, iters, libm_tanh_f64),
        ("exp", "libm", "f64") => run_f64(kind, op_name, iters, libm_exp_f64),
        ("log", "libm", "f64") => run_f64(kind, op_name, iters, libm_log_f64),
        ("log1p", "libm", "f64") => run_f64(kind, op_name, iters, libm_log1p_f64),
        ("expm1", "libm", "f64") => run_f64(kind, op_name, iters, libm_expm1_f64),

        // ---------- libm f16 ----------
        ("sin", "libm", "f16") => run_f16(kind, op_name, iters, libm_sin_f16),
        ("cos", "libm", "f16") => run_f16(kind, op_name, iters, libm_cos_f16),
        ("tan", "libm", "f16") => run_f16(kind, op_name, iters, libm_tan_f16),
        ("asin", "libm", "f16") => run_f16(kind, op_name, iters, libm_asin_f16),
        ("acos", "libm", "f16") => run_f16(kind, op_name, iters, libm_acos_f16),
        ("atan", "libm", "f16") => run_f16(kind, op_name, iters, libm_atan_f16),
        ("sinh", "libm", "f16") => run_f16(kind, op_name, iters, libm_sinh_f16),
        ("cosh", "libm", "f16") => run_f16(kind, op_name, iters, libm_cosh_f16),
        ("tanh", "libm", "f16") => run_f16(kind, op_name, iters, libm_tanh_f16),
        ("exp", "libm", "f16") => run_f16(kind, op_name, iters, libm_exp_f16),
        ("log", "libm", "f16") => run_f16(kind, op_name, iters, libm_log_f16),
        ("log1p", "libm", "f16") => run_f16(kind, op_name, iters, libm_log1p_f16),
        ("expm1", "libm", "f16") => run_f16(kind, op_name, iters, libm_expm1_f16),

        // ---------- native f32 ----------
        ("sin", "nat", "f32") => run_f32(kind, op_name, iters, nat_sin_f32),
        ("cos", "nat", "f32") => run_f32(kind, op_name, iters, nat_cos_f32),
        ("tan", "nat", "f32") => run_f32(kind, op_name, iters, nat_tan_f32),
        ("asin", "nat", "f32") => run_f32(kind, op_name, iters, nat_asin_f32),
        ("acos", "nat", "f32") => run_f32(kind, op_name, iters, nat_acos_f32),
        ("atan", "nat", "f32") => run_f32(kind, op_name, iters, nat_atan_f32),
        ("sinh", "nat", "f32") => run_f32(kind, op_name, iters, nat_sinh_f32),
        ("cosh", "nat", "f32") => run_f32(kind, op_name, iters, nat_cosh_f32),
        ("tanh", "nat", "f32") => run_f32(kind, op_name, iters, nat_tanh_f32),
        ("exp", "nat", "f32") => run_f32(kind, op_name, iters, nat_exp_f32),
        ("log", "nat", "f32") => run_f32(kind, op_name, iters, nat_log_f32),
        ("log1p", "nat", "f32") => run_f32(kind, op_name, iters, nat_log1p_f32),
        ("expm1", "nat", "f32") => run_f32(kind, op_name, iters, nat_expm1_f32),

        // ---------- native f64 ----------
        ("sin", "nat", "f64") => run_f64(kind, op_name, iters, nat_sin_f64),
        ("cos", "nat", "f64") => run_f64(kind, op_name, iters, nat_cos_f64),
        ("tan", "nat", "f64") => run_f64(kind, op_name, iters, nat_tan_f64),
        ("asin", "nat", "f64") => run_f64(kind, op_name, iters, nat_asin_f64),
        ("acos", "nat", "f64") => run_f64(kind, op_name, iters, nat_acos_f64),
        ("atan", "nat", "f64") => run_f64(kind, op_name, iters, nat_atan_f64),
        ("sinh", "nat", "f64") => run_f64(kind, op_name, iters, nat_sinh_f64),
        ("cosh", "nat", "f64") => run_f64(kind, op_name, iters, nat_cosh_f64),
        ("tanh", "nat", "f64") => run_f64(kind, op_name, iters, nat_tanh_f64),
        ("exp", "nat", "f64") => run_f64(kind, op_name, iters, nat_exp_f64),
        ("log", "nat", "f64") => run_f64(kind, op_name, iters, nat_log_f64),
        ("log1p", "nat", "f64") => run_f64(kind, op_name, iters, nat_log1p_f64),
        ("expm1", "nat", "f64") => run_f64(kind, op_name, iters, nat_expm1_f64),

        // ---------- native f16 ----------
        ("sin", "nat", "f16") => run_f16(kind, op_name, iters, nat_sin_f16),
        ("cos", "nat", "f16") => run_f16(kind, op_name, iters, nat_cos_f16),
        ("tan", "nat", "f16") => run_f16(kind, op_name, iters, nat_tan_f16),
        ("asin", "nat", "f16") => run_f16(kind, op_name, iters, nat_asin_f16),
        ("acos", "nat", "f16") => run_f16(kind, op_name, iters, nat_acos_f16),
        ("atan", "nat", "f16") => run_f16(kind, op_name, iters, nat_atan_f16),
        ("sinh", "nat", "f16") => run_f16(kind, op_name, iters, nat_sinh_f16),
        ("cosh", "nat", "f16") => run_f16(kind, op_name, iters, nat_cosh_f16),
        ("tanh", "nat", "f16") => run_f16(kind, op_name, iters, nat_tanh_f16),
        ("exp", "nat", "f16") => run_f16(kind, op_name, iters, nat_exp_f16),
        ("log", "nat", "f16") => run_f16(kind, op_name, iters, nat_log_f16),
        ("log1p", "nat", "f16") => run_f16(kind, op_name, iters, nat_log1p_f16),
        ("expm1", "nat", "f16") => run_f16(kind, op_name, iters, nat_expm1_f16),

        // ---------- native f128 ----------
        #[cfg(not(target_env = "musl"))]
        (_, "nat", "f128") => native_f128::run(kind, op_name, iters),

        // ---------- x87 f64-only ----------
        #[cfg(fpx_x87)]
        ("sin", "x87", "f64") => run_x87(kind, op_name, iters, fpx_sin_x87),
        #[cfg(fpx_x87)]
        ("cos", "x87", "f64") => run_x87(kind, op_name, iters, fpx_cos_x87),
        #[cfg(fpx_x87)]
        ("tan", "x87", "f64") => run_x87(kind, op_name, iters, fpx_tan_x87),
        #[cfg(fpx_x87)]
        ("sqrt", "x87", "f64") => run_x87(kind, op_name, iters, fpx_sqrt_x87),
        #[cfg(not(fpx_x87))]
        (_, "x87", _) => {
            eprintln!("x87 mode is unavailable: enable the x87 feature on a supported x86 target");
            std::process::exit(2);
        }

        _ => usage(),
    };

    consume(out);
}
