use aether_rand::randomizers::XorShift64Star;
use core::hint::{black_box, spin_loop};
use rtt_target::{rprintln, rtt_init_print};

const RUNS_PER_FUNCTION: usize = 100;

#[cfg(target_arch = "arm")]
#[inline(always)]
fn init_cycle_counter() {
    const DEMCR: *mut u32 = 0xE000_EDFC as *mut u32;
    const DWT_CTRL: *mut u32 = 0xE000_1000 as *mut u32;
    const DWT_CYCCNT: *mut u32 = 0xE000_1004 as *mut u32;

    unsafe {
        let demcr = core::ptr::read_volatile(DEMCR);
        core::ptr::write_volatile(DEMCR, demcr | (1 << 24));
        core::ptr::write_volatile(DWT_CYCCNT, 0);
        let dwt_ctrl = core::ptr::read_volatile(DWT_CTRL);
        core::ptr::write_volatile(DWT_CTRL, dwt_ctrl | 1);
    }
}

#[cfg(target_arch = "arm")]
#[inline(always)]
fn read_cycle_counter() -> u64 {
    const DWT_CYCCNT: *const u32 = 0xE000_1004 as *const u32;
    unsafe { core::ptr::read_volatile(DWT_CYCCNT) as u64 }
}

#[cfg(target_arch = "riscv32")]
#[inline(always)]
fn init_cycle_counter() {}

#[cfg(target_arch = "riscv32")]
#[inline(always)]
fn read_cycle_counter() -> u64 {
    let cycles: u32;
    unsafe {
        core::arch::asm!("rdcycle {0}", out(reg) cycles, options(nomem, nostack, preserves_flags));
    }
    cycles as u64
}

#[cfg(not(any(target_arch = "arm", target_arch = "riscv32")))]
#[inline(always)]
fn init_cycle_counter() {}

#[cfg(not(any(target_arch = "arm", target_arch = "riscv32")))]
#[inline(always)]
fn read_cycle_counter() -> u64 {
    0
}

#[inline(always)]
fn sample_uniform(rng: &mut XorShift64Star, min: f32, max: f32) -> f32 {
    let t = rng.next_f64() as f32;
    min + (max - min) * t
}

#[inline(never)]
fn bench_function(rng: &mut XorShift64Star, name: &str, min: f32, max: f32, f: fn(f32) -> f32) {
    let batch_start = read_cycle_counter();

    for iteration in 0..RUNS_PER_FUNCTION {
        let input = sample_uniform(rng, min, max);
        let start = read_cycle_counter();
        let output = black_box(f(black_box(input)));
        let end = read_cycle_counter();
        let elapsed = end.wrapping_sub(start);

        rprintln!(
            "sample,{},{},{:.9e},{:.9e},{}",
            name,
            iteration,
            input,
            output,
            elapsed
        );
    }

    let batch_end = read_cycle_counter();
    let batch_elapsed = batch_end.wrapping_sub(batch_start);
    rprintln!("batch,{},{}", name, batch_elapsed);
}

#[inline(never)]
pub fn run_trig_bench_rtt() -> ! {
    rtt_init_print!();
    init_cycle_counter();

    let mut rng = XorShift64Star::new(0xF0A7_2350_A5A5_0001);

    rprintln!("meta,runs_per_function,{}", RUNS_PER_FUNCTION);
    rprintln!("header,sample,function,iter,input,output,cycles");
    rprintln!("header,batch,function,total_cycles");

    bench_function(
        &mut rng,
        "sinf",
        aether_fpx::trig::sin::SIN_INTERVAL_MIN_F32,
        aether_fpx::trig::sin::SIN_INTERVAL_MAX_F32,
        aether_fpx::trig::sinf,
    );
    bench_function(
        &mut rng,
        "cosf",
        aether_fpx::trig::cos::COS_INTERVAL_MIN_F32,
        aether_fpx::trig::cos::COS_INTERVAL_MAX_F32,
        aether_fpx::trig::cosf,
    );
    bench_function(
        &mut rng,
        "tanf",
        aether_fpx::trig::tan::TAN_INTERVAL_MIN_F32,
        aether_fpx::trig::tan::TAN_INTERVAL_MAX_F32,
        aether_fpx::trig::tanf,
    );
    bench_function(
        &mut rng,
        "asinf",
        aether_fpx::trig::asin::ASIN_INTERVAL_MIN_F32,
        aether_fpx::trig::asin::ASIN_INTERVAL_MAX_F32,
        aether_fpx::trig::asinf,
    );
    bench_function(
        &mut rng,
        "acosf",
        aether_fpx::trig::acos::ACOS_INTERVAL_MIN_F32,
        aether_fpx::trig::acos::ACOS_INTERVAL_MAX_F32,
        aether_fpx::trig::acosf,
    );
    bench_function(
        &mut rng,
        "atanf",
        aether_fpx::trig::atan::ATAN_INTERVAL_MIN_F32,
        aether_fpx::trig::atan::ATAN_INTERVAL_MAX_F32,
        aether_fpx::trig::atanf,
    );
    bench_function(
        &mut rng,
        "sinhf",
        aether_fpx::trig::sinh::SINH_INTERVAL_MIN_F32,
        aether_fpx::trig::sinh::SINH_INTERVAL_MAX_F32,
        aether_fpx::trig::sinhf,
    );
    bench_function(
        &mut rng,
        "coshf",
        aether_fpx::trig::cosh::COSH_INTERVAL_MIN_F32,
        aether_fpx::trig::cosh::COSH_INTERVAL_MAX_F32,
        aether_fpx::trig::coshf,
    );
    bench_function(
        &mut rng,
        "tanhf",
        aether_fpx::trig::tanh::TANH_INTERVAL_MIN_F32,
        aether_fpx::trig::tanh::TANH_INTERVAL_MAX_F32,
        aether_fpx::trig::tanhf,
    );
    bench_function(
        &mut rng,
        "expf",
        aether_fpx::trig::exp::EXP_INTERVAL_MIN_F32,
        aether_fpx::trig::exp::EXP_INTERVAL_MAX_F32,
        aether_fpx::trig::expf,
    );
    bench_function(
        &mut rng,
        "expm1f",
        aether_fpx::trig::expm1::EXPM1_INTERVAL_MIN_F32,
        aether_fpx::trig::expm1::EXPM1_INTERVAL_MAX_F32,
        aether_fpx::trig::expm1f,
    );
    bench_function(
        &mut rng,
        "logf",
        aether_fpx::trig::log::LOG_INTERVAL_MIN_F32,
        aether_fpx::trig::log::LOG_INTERVAL_MAX_F32,
        aether_fpx::trig::logf,
    );
    bench_function(
        &mut rng,
        "log1pf",
        aether_fpx::trig::log1p::LOG1P_INTERVAL_MIN_F32,
        aether_fpx::trig::log1p::LOG1P_INTERVAL_MAX_F32,
        aether_fpx::trig::log1pf,
    );

    rprintln!("meta,done,1");

    loop {
        spin_loop();
    }
}
