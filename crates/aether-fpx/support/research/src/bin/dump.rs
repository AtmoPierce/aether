#![no_main]

use std::hint::black_box;

/* ---------------- x87 externs ---------------- */

#[cfg(fpx_x87)]
extern "C" {
    fn fpx_sin_x87(x: f64) -> f64;
    fn fpx_cos_x87(x: f64) -> f64;
    fn fpx_sqrt_x87(x: f64) -> f64;
    fn fpx_tan_x87(x: f64) -> f64;
}

/* ---------------- FPX wrappers ---------------- */

#[no_mangle]
#[inline(never)]
pub extern "C" fn touch_fpx() {
    black_box(aether_fpx::sinf(0.5f32));
    black_box(aether_fpx::cosf(0.5f32));
    black_box(aether_fpx::tanf(0.5f32));

    black_box(aether_fpx::sin(0.5f64));
    black_box(aether_fpx::cos(0.5f64));
    black_box(aether_fpx::tan(0.5f64));
}

/* ---------------- direct-kernel assembly anchors ---------------- */

macro_rules! dump_kernel {
    ($name:ident, $function:path, $ty:ty) => {
        #[no_mangle]
        #[inline(never)]
        pub extern "C" fn $name(x: $ty) -> $ty {
            black_box($function(black_box(x)))
        }
    };
}

macro_rules! dump_kernel_set {
    (
        $horner_f32:ident => $horner_kernel_f32:ident,
        $estrin_f32:ident => $estrin_kernel_f32:ident,
        $horner_f64:ident => $horner_kernel_f64:ident,
        $estrin_f64:ident => $estrin_kernel_f64:ident $(,)?
    ) => {
        dump_kernel!($horner_f32, aether_fpx::bench::$horner_kernel_f32, f32);
        dump_kernel!($estrin_f32, aether_fpx::bench::$estrin_kernel_f32, f32);
        dump_kernel!($horner_f64, aether_fpx::bench::$horner_kernel_f64, f64);
        dump_kernel!($estrin_f64, aether_fpx::bench::$estrin_kernel_f64, f64);
    };
}

dump_kernel_set!(
    dump_sin_horner_f32 => sin_kernel_horner_f32,
    dump_sin_estrin_f32 => sin_kernel_estrin_f32,
    dump_sin_horner_f64 => sin_kernel_horner_f64,
    dump_sin_estrin_f64 => sin_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_cos_horner_f32 => cos_kernel_horner_f32,
    dump_cos_estrin_f32 => cos_kernel_estrin_f32,
    dump_cos_horner_f64 => cos_kernel_horner_f64,
    dump_cos_estrin_f64 => cos_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_tan_horner_f32 => tan_kernel_horner_f32,
    dump_tan_estrin_f32 => tan_kernel_estrin_f32,
    dump_tan_horner_f64 => tan_kernel_horner_f64,
    dump_tan_estrin_f64 => tan_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_asin_horner_f32 => asin_kernel_horner_f32,
    dump_asin_estrin_f32 => asin_kernel_estrin_f32,
    dump_asin_horner_f64 => asin_kernel_horner_f64,
    dump_asin_estrin_f64 => asin_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_acos_horner_f32 => acos_kernel_horner_f32,
    dump_acos_estrin_f32 => acos_kernel_estrin_f32,
    dump_acos_horner_f64 => acos_kernel_horner_f64,
    dump_acos_estrin_f64 => acos_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_atan_horner_f32 => atan_kernel_horner_f32,
    dump_atan_estrin_f32 => atan_kernel_estrin_f32,
    dump_atan_horner_f64 => atan_kernel_horner_f64,
    dump_atan_estrin_f64 => atan_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_sinh_horner_f32 => sinh_kernel_horner_f32,
    dump_sinh_estrin_f32 => sinh_kernel_estrin_f32,
    dump_sinh_horner_f64 => sinh_kernel_horner_f64,
    dump_sinh_estrin_f64 => sinh_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_cosh_horner_f32 => cosh_kernel_horner_f32,
    dump_cosh_estrin_f32 => cosh_kernel_estrin_f32,
    dump_cosh_horner_f64 => cosh_kernel_horner_f64,
    dump_cosh_estrin_f64 => cosh_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_tanh_horner_f32 => tanh_kernel_horner_f32,
    dump_tanh_estrin_f32 => tanh_kernel_estrin_f32,
    dump_tanh_horner_f64 => tanh_kernel_horner_f64,
    dump_tanh_estrin_f64 => tanh_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_exp_horner_f32 => exp_kernel_horner_f32,
    dump_exp_estrin_f32 => exp_kernel_estrin_f32,
    dump_exp_horner_f64 => exp_kernel_horner_f64,
    dump_exp_estrin_f64 => exp_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_expm1_horner_f32 => expm1_kernel_horner_f32,
    dump_expm1_estrin_f32 => expm1_kernel_estrin_f32,
    dump_expm1_horner_f64 => expm1_kernel_horner_f64,
    dump_expm1_estrin_f64 => expm1_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_log_horner_f32 => log_kernel_horner_f32,
    dump_log_estrin_f32 => log_kernel_estrin_f32,
    dump_log_horner_f64 => log_kernel_horner_f64,
    dump_log_estrin_f64 => log_kernel_estrin_f64,
);
dump_kernel_set!(
    dump_log1p_horner_f32 => log1p_kernel_horner_f32,
    dump_log1p_estrin_f32 => log1p_kernel_estrin_f32,
    dump_log1p_horner_f64 => log1p_kernel_horner_f64,
    dump_log1p_estrin_f64 => log1p_kernel_estrin_f64,
);

/* ---------------- measured FPX raw-polynomial path anchors ---------------- */

#[cfg(feature = "simd")]
macro_rules! dump_fpx_raw_path {
    ($name:ident, $function:path, $path:path, $ty:ty, $len:expr) => {
        #[no_mangle]
        #[inline(never)]
        pub unsafe extern "C" fn $name(values: *mut $ty) {
            let values = unsafe { std::slice::from_raw_parts_mut(values, $len) };
            $function($path, black_box(values));
            black_box(values);
        }
    };
}

#[cfg(feature = "simd")]
dump_fpx_raw_path!(
    dump_sin_raw_scalar_f32,
    aether_fpx::bench_paths::eval_sin_segment0_polynomial_f32,
    aether_fpx::bench_paths::ArithmeticPath::Scalar,
    f32,
    32
);
#[cfg(feature = "simd")]
dump_fpx_raw_path!(
    dump_sin_raw_fma_f32,
    aether_fpx::bench_paths::eval_sin_segment0_polynomial_f32,
    aether_fpx::bench_paths::ArithmeticPath::Fma,
    f32,
    32
);
#[cfg(feature = "simd")]
dump_fpx_raw_path!(
    dump_sin_raw_simd_horner_f32,
    aether_fpx::bench_paths::eval_sin_segment0_polynomial_f32,
    aether_fpx::bench_paths::ArithmeticPath::Simd,
    f32,
    32
);
#[cfg(feature = "simd")]
dump_fpx_raw_path!(
    dump_sin_raw_simd_estrin_f32,
    aether_fpx::bench_paths::eval_sin_segment0_polynomial_f32,
    aether_fpx::bench_paths::ArithmeticPath::SimdEstrin,
    f32,
    32
);
#[cfg(feature = "simd")]
dump_fpx_raw_path!(
    dump_sin_raw_scalar_f64,
    aether_fpx::bench_paths::eval_sin_segment0_polynomial_f64,
    aether_fpx::bench_paths::ArithmeticPath::Scalar,
    f64,
    16
);
#[cfg(feature = "simd")]
dump_fpx_raw_path!(
    dump_sin_raw_fma_f64,
    aether_fpx::bench_paths::eval_sin_segment0_polynomial_f64,
    aether_fpx::bench_paths::ArithmeticPath::Fma,
    f64,
    16
);
#[cfg(feature = "simd")]
dump_fpx_raw_path!(
    dump_sin_raw_simd_horner_f64,
    aether_fpx::bench_paths::eval_sin_segment0_polynomial_f64,
    aether_fpx::bench_paths::ArithmeticPath::Simd,
    f64,
    16
);
#[cfg(feature = "simd")]
dump_fpx_raw_path!(
    dump_sin_raw_simd_estrin_f64,
    aether_fpx::bench_paths::eval_sin_segment0_polynomial_f64,
    aether_fpx::bench_paths::ArithmeticPath::SimdEstrin,
    f64,
    16
);

/* ---------------- native wrappers ---------------- */

#[no_mangle]
#[inline(never)]
pub extern "C" fn touch_nat() {
    black_box(0.5_f32.sin());
    black_box(0.5_f32.cos());
    black_box(0.5_f32.tan());

    black_box(0.5_f64.sin());
    black_box(0.5_f64.cos());
    black_box(0.5_f64.tan());
}

/* ---------------- x87 wrappers ---------------- */

#[no_mangle]
#[inline(never)]
#[cfg(fpx_x87)]
pub extern "C" fn touch_x87() {
    unsafe {
        black_box(fpx_sin_x87(0.5));
        black_box(fpx_cos_x87(0.5));
        black_box(fpx_sqrt_x87(0.5));
        black_box(fpx_tan_x87(0.5));
    }
}

/* ---------------- entry ---------------- */

#[no_mangle]
pub extern "C" fn main() -> i32 {
    touch_fpx();
    touch_nat();
    #[cfg(fpx_x87)]
    touch_x87();
    0
}
