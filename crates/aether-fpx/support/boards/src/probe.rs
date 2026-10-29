use core::hint::{black_box, spin_loop};

#[inline(never)]
pub fn run_trig_probe() -> ! {
    let x64 = 0.5_f64;
    let x32 = 0.5_f32;
    let x16 = 0.5_f16;
    let x128 = 0.5_f128;

    black_box(aether_fpx::trig::sin(x64));
    black_box(aether_fpx::trig::sinf(x32));
    black_box(aether_fpx::trig::sins(x16));
    black_box(aether_fpx::trig::sinl(x128));

    black_box(aether_fpx::trig::cos(x64));
    black_box(aether_fpx::trig::cosf(x32));
    black_box(aether_fpx::trig::coss(x16));
    black_box(aether_fpx::trig::cosl(x128));

    black_box(aether_fpx::trig::tan(x64));
    black_box(aether_fpx::trig::tanf(x32));
    black_box(aether_fpx::trig::tans(x16));
    black_box(aether_fpx::trig::tanl(x128));

    black_box(aether_fpx::trig::asin(x64));
    black_box(aether_fpx::trig::asinf(x32));
    black_box(aether_fpx::trig::asins(x16));
    black_box(aether_fpx::trig::asinl(x128));

    black_box(aether_fpx::trig::acos(x64));
    black_box(aether_fpx::trig::acosf(x32));
    black_box(aether_fpx::trig::acoss(x16));
    black_box(aether_fpx::trig::acosl(x128));

    black_box(aether_fpx::trig::atan(x64));
    black_box(aether_fpx::trig::atanf(x32));
    black_box(aether_fpx::trig::atans(x16));
    black_box(aether_fpx::trig::atanl(x128));

    black_box(aether_fpx::trig::sinh(x64));
    black_box(aether_fpx::trig::sinhf(x32));
    black_box(aether_fpx::trig::sinhs(x16));
    black_box(aether_fpx::trig::sinhl(x128));

    black_box(aether_fpx::trig::cosh(x64));
    black_box(aether_fpx::trig::coshf(x32));
    black_box(aether_fpx::trig::coshs(x16));
    black_box(aether_fpx::trig::coshl(x128));

    black_box(aether_fpx::trig::tanh(x64));
    black_box(aether_fpx::trig::tanhf(x32));
    black_box(aether_fpx::trig::tanhs(x16));
    black_box(aether_fpx::trig::tanhl(x128));

    black_box(aether_fpx::trig::exp(x64));
    black_box(aether_fpx::trig::expf(x32));
    black_box(aether_fpx::trig::exps(x16));
    black_box(aether_fpx::trig::expl(x128));

    black_box(aether_fpx::trig::expm1(x64));
    black_box(aether_fpx::trig::expm1f(x32));
    black_box(aether_fpx::trig::expm1s(x16));
    black_box(aether_fpx::trig::expm1l(x128));

    black_box(aether_fpx::trig::log(1.5_f64));
    black_box(aether_fpx::trig::logf(1.5_f32));
    black_box(aether_fpx::trig::logs(1.5_f16));
    black_box(aether_fpx::trig::logl(1.5_f128));

    black_box(aether_fpx::trig::log1p(x64));
    black_box(aether_fpx::trig::log1pf(x32));
    black_box(aether_fpx::trig::log1ps(x16));
    black_box(aether_fpx::trig::log1pl(x128));

    loop {
        spin_loop();
    }
}
