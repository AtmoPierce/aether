#![no_std]
#![no_main]
#![feature(f16, f128)]

use core::arch::asm;
use core::hint::{black_box, spin_loop};

use panic_halt as _;

macro_rules! wrap_f64 {
    ($name:ident, $target:path) => {
        #[unsafe(no_mangle)]
        #[inline(never)]
        pub fn $name(x: f64) -> f64 {
            $target(x)
        }
    };
}

macro_rules! wrap_f32 {
    ($name:ident, $target:path) => {
        #[unsafe(no_mangle)]
        #[inline(never)]
        pub fn $name(x: f32) -> f32 {
            $target(x)
        }
    };
}

macro_rules! wrap_f16 {
    ($name:ident, $target:path) => {
        #[unsafe(no_mangle)]
        #[inline(never)]
        pub fn $name(x: f16) -> f16 {
            $target(x)
        }
    };
}

macro_rules! wrap_f128 {
    ($name:ident, $target:path) => {
        #[unsafe(no_mangle)]
        #[inline(never)]
        pub fn $name(x: f128) -> f128 {
            $target(x)
        }
    };
}

wrap_f64!(call_sin_f64, aether_fpx::trig::sin);
wrap_f32!(call_sin_f32, aether_fpx::trig::sinf);
wrap_f16!(call_sin_f16, aether_fpx::trig::sins);
wrap_f128!(call_sin_f128, aether_fpx::trig::sinl);

wrap_f64!(call_cos_f64, aether_fpx::trig::cos);
wrap_f32!(call_cos_f32, aether_fpx::trig::cosf);
wrap_f16!(call_cos_f16, aether_fpx::trig::coss);
wrap_f128!(call_cos_f128, aether_fpx::trig::cosl);

wrap_f64!(call_tan_f64, aether_fpx::trig::tan);
wrap_f32!(call_tan_f32, aether_fpx::trig::tanf);
wrap_f16!(call_tan_f16, aether_fpx::trig::tans);
wrap_f128!(call_tan_f128, aether_fpx::trig::tanl);

wrap_f64!(call_asin_f64, aether_fpx::trig::asin);
wrap_f32!(call_asin_f32, aether_fpx::trig::asinf);
wrap_f16!(call_asin_f16, aether_fpx::trig::asins);
wrap_f128!(call_asin_f128, aether_fpx::trig::asinl);

wrap_f64!(call_acos_f64, aether_fpx::trig::acos);
wrap_f32!(call_acos_f32, aether_fpx::trig::acosf);
wrap_f16!(call_acos_f16, aether_fpx::trig::acoss);
wrap_f128!(call_acos_f128, aether_fpx::trig::acosl);

wrap_f64!(call_atan_f64, aether_fpx::trig::atan);
wrap_f32!(call_atan_f32, aether_fpx::trig::atanf);
wrap_f16!(call_atan_f16, aether_fpx::trig::atans);
wrap_f128!(call_atan_f128, aether_fpx::trig::atanl);

wrap_f64!(call_sinh_f64, aether_fpx::trig::sinh);
wrap_f32!(call_sinh_f32, aether_fpx::trig::sinhf);
wrap_f16!(call_sinh_f16, aether_fpx::trig::sinhs);
wrap_f128!(call_sinh_f128, aether_fpx::trig::sinhl);

wrap_f64!(call_cosh_f64, aether_fpx::trig::cosh);
wrap_f32!(call_cosh_f32, aether_fpx::trig::coshf);
wrap_f16!(call_cosh_f16, aether_fpx::trig::coshs);
wrap_f128!(call_cosh_f128, aether_fpx::trig::coshl);

wrap_f64!(call_tanh_f64, aether_fpx::trig::tanh);
wrap_f32!(call_tanh_f32, aether_fpx::trig::tanhf);
wrap_f16!(call_tanh_f16, aether_fpx::trig::tanhs);
wrap_f128!(call_tanh_f128, aether_fpx::trig::tanhl);

wrap_f64!(call_exp_f64, aether_fpx::trig::exp);
wrap_f32!(call_exp_f32, aether_fpx::trig::expf);
wrap_f16!(call_exp_f16, aether_fpx::trig::exps);
wrap_f128!(call_exp_f128, aether_fpx::trig::expl);

wrap_f64!(call_expm1_f64, aether_fpx::trig::expm1);
wrap_f32!(call_expm1_f32, aether_fpx::trig::expm1f);
wrap_f16!(call_expm1_f16, aether_fpx::trig::expm1s);
wrap_f128!(call_expm1_f128, aether_fpx::trig::expm1l);

wrap_f64!(call_log_f64, aether_fpx::trig::log);
wrap_f32!(call_log_f32, aether_fpx::trig::logf);
wrap_f16!(call_log_f16, aether_fpx::trig::logs);
wrap_f128!(call_log_f128, aether_fpx::trig::logl);

wrap_f64!(call_log1p_f64, aether_fpx::trig::log1p);
wrap_f32!(call_log1p_f32, aether_fpx::trig::log1pf);
wrap_f16!(call_log1p_f16, aether_fpx::trig::log1ps);
wrap_f128!(call_log1p_f128, aether_fpx::trig::log1pl);

#[cfg(target_arch = "arm")]
#[inline(never)]
unsafe fn probe_vcx3a_f<const IMM: u8>(sd: &mut f32, sn: f32, sm: f32) {
    asm!(
        "vcx3a.f p0, s0, s0, s1, s2, #{imm}",
        imm = const IMM,
        inout("s0") *sd,
        in("s1") sn,
        in("s2") sm,
        options(nostack)
    );
}

#[cfg(not(target_arch = "arm"))]
#[inline(never)]
unsafe fn probe_vcx3a_f<const IMM: u8>(sd: &mut f32, sn: f32, sm: f32) {
    let _ = IMM;
    let _ = sd;
    let _ = sn;
    let _ = sm;
}

#[unsafe(no_mangle)]
pub extern "C" fn main() -> ! {
    let x64 = 0.5_f64;
    let x32 = 0.5_f32;
    let x16 = 0.5_f16;
    let x128 = 0.5_f128;

    let mut vcx_acc = 0.25_f32;
    unsafe {
        probe_vcx3a_f::<0>(&mut vcx_acc, x32, x32);
    }
    black_box(vcx_acc);

    black_box(call_sin_f64(x64));
    black_box(call_sin_f32(x32));
    black_box(call_sin_f16(x16));
    black_box(call_sin_f128(x128));

    black_box(call_cos_f64(x64));
    black_box(call_cos_f32(x32));
    black_box(call_cos_f16(x16));
    black_box(call_cos_f128(x128));

    black_box(call_tan_f64(x64));
    black_box(call_tan_f32(x32));
    black_box(call_tan_f16(x16));
    black_box(call_tan_f128(x128));

    black_box(call_asin_f64(x64));
    black_box(call_asin_f32(x32));
    black_box(call_asin_f16(x16));
    black_box(call_asin_f128(x128));

    black_box(call_acos_f64(x64));
    black_box(call_acos_f32(x32));
    black_box(call_acos_f16(x16));
    black_box(call_acos_f128(x128));

    black_box(call_atan_f64(x64));
    black_box(call_atan_f32(x32));
    black_box(call_atan_f16(x16));
    black_box(call_atan_f128(x128));

    black_box(call_sinh_f64(x64));
    black_box(call_sinh_f32(x32));
    black_box(call_sinh_f16(x16));
    black_box(call_sinh_f128(x128));

    black_box(call_cosh_f64(x64));
    black_box(call_cosh_f32(x32));
    black_box(call_cosh_f16(x16));
    black_box(call_cosh_f128(x128));

    black_box(call_tanh_f64(x64));
    black_box(call_tanh_f32(x32));
    black_box(call_tanh_f16(x16));
    black_box(call_tanh_f128(x128));

    black_box(call_exp_f64(x64));
    black_box(call_exp_f32(x32));
    black_box(call_exp_f16(x16));
    black_box(call_exp_f128(x128));

    black_box(call_expm1_f64(x64));
    black_box(call_expm1_f32(x32));
    black_box(call_expm1_f16(x16));
    black_box(call_expm1_f128(x128));

    black_box(call_log_f64(1.5_f64));
    black_box(call_log_f32(1.5_f32));
    black_box(call_log_f16(1.5_f16));
    black_box(call_log_f128(1.5_f128));

    black_box(call_log1p_f64(x64));
    black_box(call_log1p_f32(x32));
    black_box(call_log1p_f16(x16));
    black_box(call_log1p_f128(x128));

    loop {
        spin_loop();
    }
}
