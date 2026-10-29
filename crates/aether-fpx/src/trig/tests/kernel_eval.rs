use crate::{bench, trig};

struct KernelF32 {
    name: &'static str,
    horner: fn(f32) -> f32,
    estrin: fn(f32) -> f32,
    lo: f32,
    hi: f32,
}

struct KernelF64 {
    name: &'static str,
    horner: fn(f64) -> f64,
    estrin: fn(f64) -> f64,
    lo: f64,
    hi: f64,
}

fn assert_f32_kernel_pair(case: &KernelF32) {
    const SAMPLES: usize = 33;
    for index in 0..SAMPLES {
        let t = index as f32 / (SAMPLES - 1) as f32;
        let x = case.lo + (case.hi - case.lo) * t;
        let horner = (case.horner)(x);
        let estrin = (case.estrin)(x);
        let scale = horner.abs().max(estrin.abs()).max(1.0);
        assert!(
            (horner - estrin).abs() <= 1.0e-4 * scale,
            "{} f32 Horner/Estrin mismatch at x={x:?}: horner={horner:?}, estrin={estrin:?}",
            case.name
        );
    }
}

fn assert_f64_kernel_pair(case: &KernelF64) {
    const SAMPLES: usize = 33;
    for index in 0..SAMPLES {
        let t = index as f64 / (SAMPLES - 1) as f64;
        let x = case.lo + (case.hi - case.lo) * t;
        let horner = (case.horner)(x);
        let estrin = (case.estrin)(x);
        let scale = horner.abs().max(estrin.abs()).max(1.0);
        assert!(
            (horner - estrin).abs() <= 2.0e-12 * scale,
            "{} f64 Horner/Estrin mismatch at x={x:?}: horner={horner:?}, estrin={estrin:?}",
            case.name
        );
    }
}

#[test]
fn all_direct_f32_kernels_have_horner_and_estrin_paths() {
    let cases = [
        KernelF32 {
            name: "sin",
            horner: bench::sin_kernel_horner_f32,
            estrin: bench::sin_kernel_estrin_f32,
            lo: trig::sin::SIN_INTERVAL_MIN_F32,
            hi: trig::sin::SIN_INTERVAL_MAX_F32,
        },
        KernelF32 {
            name: "cos",
            horner: bench::cos_kernel_horner_f32,
            estrin: bench::cos_kernel_estrin_f32,
            lo: trig::cos::COS_INTERVAL_MIN_F32,
            hi: trig::cos::COS_INTERVAL_MAX_F32,
        },
        KernelF32 {
            name: "tan",
            horner: bench::tan_kernel_horner_f32,
            estrin: bench::tan_kernel_estrin_f32,
            lo: trig::tan::TAN_INTERVAL_MIN_F32,
            hi: trig::tan::TAN_INTERVAL_MAX_F32,
        },
        KernelF32 {
            name: "asin",
            horner: bench::asin_kernel_horner_f32,
            estrin: bench::asin_kernel_estrin_f32,
            lo: trig::asin::ASIN_INTERVAL_MIN_F32,
            hi: trig::asin::ASIN_INTERVAL_MAX_F32,
        },
        KernelF32 {
            name: "acos",
            horner: bench::acos_kernel_horner_f32,
            estrin: bench::acos_kernel_estrin_f32,
            lo: trig::acos::ACOS_INTERVAL_MIN_F32,
            hi: trig::acos::ACOS_INTERVAL_MAX_F32,
        },
        KernelF32 {
            name: "atan",
            horner: bench::atan_kernel_horner_f32,
            estrin: bench::atan_kernel_estrin_f32,
            lo: trig::atan::ATAN_INTERVAL_MIN_F32,
            hi: trig::atan::ATAN_INTERVAL_MAX_F32,
        },
        KernelF32 {
            name: "sinh",
            horner: bench::sinh_kernel_horner_f32,
            estrin: bench::sinh_kernel_estrin_f32,
            lo: 0.0,
            hi: 0.5,
        },
        KernelF32 {
            name: "cosh",
            horner: bench::cosh_kernel_horner_f32,
            estrin: bench::cosh_kernel_estrin_f32,
            lo: 0.0,
            hi: 0.5,
        },
        KernelF32 {
            name: "tanh",
            horner: bench::tanh_kernel_horner_f32,
            estrin: bench::tanh_kernel_estrin_f32,
            lo: 0.0,
            hi: 0.25,
        },
        KernelF32 {
            name: "exp",
            horner: bench::exp_kernel_horner_f32,
            estrin: bench::exp_kernel_estrin_f32,
            lo: trig::exp::EXP_INTERVAL_MIN_F32,
            hi: trig::exp::EXP_INTERVAL_MAX_F32,
        },
        KernelF32 {
            name: "expm1",
            horner: bench::expm1_kernel_horner_f32,
            estrin: bench::expm1_kernel_estrin_f32,
            lo: trig::expm1::EXPM1_INTERVAL_MIN_F32,
            hi: trig::expm1::EXPM1_INTERVAL_MAX_F32,
        },
        KernelF32 {
            name: "log",
            horner: bench::log_kernel_horner_f32,
            estrin: bench::log_kernel_estrin_f32,
            lo: trig::log::LOG_INTERVAL_MIN_F32,
            hi: trig::log::LOG_INTERVAL_MAX_F32,
        },
        KernelF32 {
            name: "log1p",
            horner: bench::log1p_kernel_horner_f32,
            estrin: bench::log1p_kernel_estrin_f32,
            lo: trig::log1p::LOG1P_INTERVAL_MIN_F32,
            hi: trig::log1p::LOG1P_INTERVAL_MAX_F32,
        },
    ];

    for case in &cases {
        assert_f32_kernel_pair(case);
    }
}

#[test]
fn all_direct_f64_kernels_have_horner_and_estrin_paths() {
    let cases = [
        KernelF64 {
            name: "sin",
            horner: bench::sin_kernel_horner_f64,
            estrin: bench::sin_kernel_estrin_f64,
            lo: trig::sin::SIN_INTERVAL_MIN_F64,
            hi: trig::sin::SIN_INTERVAL_MAX_F64,
        },
        KernelF64 {
            name: "cos",
            horner: bench::cos_kernel_horner_f64,
            estrin: bench::cos_kernel_estrin_f64,
            lo: trig::cos::COS_INTERVAL_MIN_F64,
            hi: trig::cos::COS_INTERVAL_MAX_F64,
        },
        KernelF64 {
            name: "tan",
            horner: bench::tan_kernel_horner_f64,
            estrin: bench::tan_kernel_estrin_f64,
            lo: trig::tan::TAN_INTERVAL_MIN_F64,
            hi: trig::tan::TAN_INTERVAL_MAX_F64,
        },
        KernelF64 {
            name: "asin",
            horner: bench::asin_kernel_horner_f64,
            estrin: bench::asin_kernel_estrin_f64,
            lo: trig::asin::ASIN_INTERVAL_MIN_F64,
            hi: trig::asin::ASIN_INTERVAL_MAX_F64,
        },
        KernelF64 {
            name: "acos",
            horner: bench::acos_kernel_horner_f64,
            estrin: bench::acos_kernel_estrin_f64,
            lo: trig::acos::ACOS_INTERVAL_MIN_F64,
            hi: trig::acos::ACOS_INTERVAL_MAX_F64,
        },
        KernelF64 {
            name: "atan",
            horner: bench::atan_kernel_horner_f64,
            estrin: bench::atan_kernel_estrin_f64,
            lo: trig::atan::ATAN_INTERVAL_MIN_F64,
            hi: trig::atan::ATAN_INTERVAL_MAX_F64,
        },
        KernelF64 {
            name: "sinh",
            horner: bench::sinh_kernel_horner_f64,
            estrin: bench::sinh_kernel_estrin_f64,
            lo: 0.0,
            hi: 0.125,
        },
        KernelF64 {
            name: "cosh",
            horner: bench::cosh_kernel_horner_f64,
            estrin: bench::cosh_kernel_estrin_f64,
            lo: 0.0,
            hi: 0.125,
        },
        KernelF64 {
            name: "tanh",
            horner: bench::tanh_kernel_horner_f64,
            estrin: bench::tanh_kernel_estrin_f64,
            lo: 0.0,
            hi: 0.0625,
        },
        KernelF64 {
            name: "exp",
            horner: bench::exp_kernel_horner_f64,
            estrin: bench::exp_kernel_estrin_f64,
            lo: trig::exp::EXP_INTERVAL_MIN_F64,
            hi: trig::exp::EXP_INTERVAL_MAX_F64,
        },
        KernelF64 {
            name: "expm1",
            horner: bench::expm1_kernel_horner_f64,
            estrin: bench::expm1_kernel_estrin_f64,
            lo: trig::expm1::EXPM1_INTERVAL_MIN_F64,
            hi: trig::expm1::EXPM1_INTERVAL_MAX_F64,
        },
        KernelF64 {
            name: "log",
            horner: bench::log_kernel_horner_f64,
            estrin: bench::log_kernel_estrin_f64,
            lo: trig::log::LOG_INTERVAL_MIN_F64,
            hi: trig::log::LOG_INTERVAL_MAX_F64,
        },
        KernelF64 {
            name: "log1p",
            horner: bench::log1p_kernel_horner_f64,
            estrin: bench::log1p_kernel_estrin_f64,
            lo: trig::log1p::LOG1P_INTERVAL_MIN_F64,
            hi: trig::log1p::LOG1P_INTERVAL_MAX_F64,
        },
    ];

    for case in &cases {
        assert_f64_kernel_pair(case);
    }
}
