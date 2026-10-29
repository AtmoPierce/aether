#![feature(f16, f128)]

use core::f64::consts::PI;
use core::{f128, f16};

use aether_fpx_bench::precision::{
    test_f128_stratified_bits, test_f16_exhaustive, test_f32_stratified_bits, test_f64_boundaries,
    test_f64_pi_refinement, test_f64_stratified_bits, write_precision_csv, PrecisionReport,
    PrecisionSummary,
};
use aether_fpx_bench::{valid_ulp_dist_f128, valid_ulp_dist_f64};
use rayon::prelude::*;

/* ============================================================
 * Function descriptors
 * ============================================================ */

#[derive(Clone, Copy)]
struct FuncF64 {
    name: &'static str,
    reference: fn(f64) -> f64,
    fpx: fn(f64) -> f64,
    lo: f64,
    hi: f64,
}

#[derive(Clone, Copy)]
struct FuncF32 {
    name: &'static str,
    reference: fn(f32) -> f32,
    fpx: fn(f32) -> f32,
    lo: f32,
    hi: f32,
}

#[derive(Clone, Copy)]
struct FuncF16 {
    name: &'static str,
    reference: fn(f16) -> f16,
    fpx: fn(f16) -> f16,
    lo: f16,
    hi: f16,
}

#[derive(Clone, Copy)]
struct FuncF128 {
    name: &'static str,
    reference: fn(f128) -> f128,
    fpx: fn(f128) -> f128,
    lo: f128,
    hi: f128,
}

fn publication_interval(name: &str) -> (f64, f64) {
    match name {
        "sin" | "cos" => (-4.0, 4.0),
        "tan" => (-1.2, 1.2),
        "asin" | "acos" => (-1.0, 1.0),
        "atan" => (-8.0, 8.0),
        "sinh" | "cosh" | "tanh" | "exp" | "expm1" => (-4.0, 4.0),
        "log" => (1.0e-6, 10.0),
        "log1p" => (-0.999_999, 10.0),
        _ => unreachable!("missing publication interval for {name}"),
    }
}

fn test_f64_publication_interval(function: FuncF64) -> PrecisionReport {
    const SAMPLES: usize = 2001;
    let (lo, hi) = publication_interval(function.name);
    let mut max_abs = -1.0_f64;
    let mut max_abs_values = (0.0_f64, 0.0_f64, 0.0_f64);
    let mut max_ulp = 0_u64;
    let mut max_ulp_values = (0.0_f64, 0.0_f64, 0.0_f64);

    for index in 0..SAMPLES {
        let x = lo + (hi - lo) * index as f64 / (SAMPLES - 1) as f64;
        let reference = (function.reference)(x);
        let fpx = (function.fpx)(x);
        if !reference.is_finite() || !fpx.is_finite() {
            continue;
        }

        let abs = (fpx - reference).abs();
        if abs > max_abs {
            max_abs = abs;
            max_abs_values = (x, reference, fpx);
        }

        if let Some(ulp) = valid_ulp_dist_f64(fpx, reference) {
            if ulp > max_ulp {
                max_ulp = ulp;
                max_ulp_values = (x, reference, fpx);
            }
        }
    }

    PrecisionReport {
        report: format!(
            "{} f64 publication_interval_2001: max_abs={:.17e}, max_ulp={}",
            function.name, max_abs, max_ulp
        ),
        summary: PrecisionSummary {
            name: function.name.to_string(),
            precision: "f64".to_string(),
            test: "publication_interval_2001".to_string(),
            max_abs: format!("{max_abs:.17e}"),
            max_ulp: max_ulp.to_string(),
            max_abs_x: format!("{:.17e}", max_abs_values.0),
            max_abs_ref: format!("{:.17e}", max_abs_values.1),
            max_abs_fpx: format!("{:.17e}", max_abs_values.2),
            max_ulp_x: format!("{:.17e}", max_ulp_values.0),
            max_ulp_ref: format!("{:.17e}", max_ulp_values.1),
            max_ulp_fpx: format!("{:.17e}", max_ulp_values.2),
        },
    }
}

fn test_f128_publication_interval(function: FuncF128) -> PrecisionReport {
    const SAMPLES: usize = 2001;
    let (lo, hi) = publication_interval(function.name);
    let lo = lo as f128;
    let hi = hi as f128;
    let mut max_abs = -1.0_f128;
    let mut max_abs_values = (0.0_f128, 0.0_f128, 0.0_f128);
    let mut max_ulp = 0_u128;
    let mut max_ulp_values = (0.0_f128, 0.0_f128, 0.0_f128);

    for index in 0..SAMPLES {
        let x = lo + (hi - lo) * (index as f128) / ((SAMPLES - 1) as f128);
        let reference = (function.reference)(x);
        let fpx = (function.fpx)(x);
        if !reference.is_finite() || !fpx.is_finite() {
            continue;
        }

        let abs = (fpx - reference).abs();
        if abs > max_abs {
            max_abs = abs;
            max_abs_values = (x, reference, fpx);
        }

        if let Some(ulp) = valid_ulp_dist_f128(fpx, reference) {
            if ulp > max_ulp {
                max_ulp = ulp;
                max_ulp_values = (x, reference, fpx);
            }
        }
    }

    PrecisionReport {
        report: format!(
            "{} f128 publication_interval_2001: max_abs={:?}, max_ulp={}",
            function.name, max_abs, max_ulp
        ),
        summary: PrecisionSummary {
            name: function.name.to_string(),
            precision: "f128".to_string(),
            test: "publication_interval_2001".to_string(),
            max_abs: format!("{max_abs:?}"),
            max_ulp: max_ulp.to_string(),
            max_abs_x: format!("{:?}", max_abs_values.0),
            max_abs_ref: format!("{:?}", max_abs_values.1),
            max_abs_fpx: format!("{:?}", max_abs_values.2),
            max_ulp_x: format!("{:?}", max_ulp_values.0),
            max_ulp_ref: format!("{:?}", max_ulp_values.1),
            max_ulp_fpx: format!("{:?}", max_ulp_values.2),
        },
    }
}

/* ============================================================
 * References
 * ============================================================ */

#[inline]
fn ref_sin_f64(x: f64) -> f64 {
    x.sin()
}
#[inline]
fn ref_cos_f64(x: f64) -> f64 {
    x.cos()
}
#[inline]
fn ref_tan_f64(x: f64) -> f64 {
    x.tan()
}
#[inline]
fn ref_sinh_f64(x: f64) -> f64 {
    x.sinh()
}
#[inline]
fn ref_cosh_f64(x: f64) -> f64 {
    x.cosh()
}
#[inline]
fn ref_tanh_f64(x: f64) -> f64 {
    x.tanh()
}
#[inline]
fn ref_asin_f64(x: f64) -> f64 {
    x.asin()
}
#[inline]
fn ref_acos_f64(x: f64) -> f64 {
    x.acos()
}
#[inline]
fn ref_atan_f64(x: f64) -> f64 {
    x.atan()
}
#[inline]
fn ref_exp_f64(x: f64) -> f64 {
    x.exp()
}
#[inline]
fn ref_log_f64(x: f64) -> f64 {
    x.ln()
}
#[inline]
fn ref_log1p_f64(x: f64) -> f64 {
    x.ln_1p()
}
#[inline]
fn ref_expm1_f64(x: f64) -> f64 {
    x.exp_m1()
}

#[inline]
fn ref_sin_f32(x: f32) -> f32 {
    x.sin()
}
#[inline]
fn ref_cos_f32(x: f32) -> f32 {
    x.cos()
}
#[inline]
fn ref_tan_f32(x: f32) -> f32 {
    x.tan()
}
#[inline]
fn ref_sinh_f32(x: f32) -> f32 {
    x.sinh()
}
#[inline]
fn ref_cosh_f32(x: f32) -> f32 {
    x.cosh()
}
#[inline]
fn ref_tanh_f32(x: f32) -> f32 {
    x.tanh()
}
#[inline]
fn ref_asin_f32(x: f32) -> f32 {
    x.asin()
}
#[inline]
fn ref_acos_f32(x: f32) -> f32 {
    x.acos()
}
#[inline]
fn ref_atan_f32(x: f32) -> f32 {
    x.atan()
}
#[inline]
fn ref_exp_f32(x: f32) -> f32 {
    x.exp()
}
#[inline]
fn ref_log_f32(x: f32) -> f32 {
    x.ln()
}
#[inline]
fn ref_log1p_f32(x: f32) -> f32 {
    x.ln_1p()
}
#[inline]
fn ref_expm1_f32(x: f32) -> f32 {
    x.exp_m1()
}

#[inline]
fn ref_sin_f16(x: f16) -> f16 {
    (x as f32).sin() as f16
}
#[inline]
fn ref_cos_f16(x: f16) -> f16 {
    (x as f32).cos() as f16
}
#[inline]
fn ref_tan_f16(x: f16) -> f16 {
    (x as f32).tan() as f16
}
#[inline]
fn ref_sinh_f16(x: f16) -> f16 {
    (x as f32).sinh() as f16
}
#[inline]
fn ref_cosh_f16(x: f16) -> f16 {
    (x as f32).cosh() as f16
}
#[inline]
fn ref_tanh_f16(x: f16) -> f16 {
    (x as f32).tanh() as f16
}
#[inline]
fn ref_asin_f16(x: f16) -> f16 {
    (x as f32).asin() as f16
}
#[inline]
fn ref_acos_f16(x: f16) -> f16 {
    (x as f32).acos() as f16
}
#[inline]
fn ref_atan_f16(x: f16) -> f16 {
    (x as f32).atan() as f16
}
#[inline]
fn ref_exp_f16(x: f16) -> f16 {
    (x as f32).exp() as f16
}
#[inline]
fn ref_log_f16(x: f16) -> f16 {
    (x as f32).ln() as f16
}
#[inline]
fn ref_log1p_f16(x: f16) -> f16 {
    (x as f32).ln_1p() as f16
}
#[inline]
fn ref_expm1_f16(x: f16) -> f16 {
    (x as f32).exp_m1() as f16
}

#[inline]
fn ref_sin_f128(x: f128) -> f128 {
    x.sin()
}
#[inline]
fn ref_cos_f128(x: f128) -> f128 {
    x.cos()
}
#[inline]
fn ref_tan_f128(x: f128) -> f128 {
    x.tan()
}
#[inline]
fn ref_sinh_f128(x: f128) -> f128 {
    x.sinh()
}
#[inline]
fn ref_cosh_f128(x: f128) -> f128 {
    x.cosh()
}
#[inline]
fn ref_tanh_f128(x: f128) -> f128 {
    x.tanh()
}
#[inline]
fn ref_asin_f128(x: f128) -> f128 {
    x.asin()
}
#[inline]
fn ref_acos_f128(x: f128) -> f128 {
    x.acos()
}
#[inline]
fn ref_atan_f128(x: f128) -> f128 {
    x.atan()
}
#[inline]
fn ref_exp_f128(x: f128) -> f128 {
    x.exp()
}
#[inline]
fn ref_log_f128(x: f128) -> f128 {
    x.ln()
}
#[inline]
fn ref_log1p_f128(x: f128) -> f128 {
    x.ln_1p()
}
#[inline]
fn ref_expm1_f128(x: f128) -> f128 {
    x.exp_m1()
}

/* ============================================================
 * FPX
 * ============================================================ */

#[inline]
fn fpx_sin_f64(x: f64) -> f64 {
    aether_fpx::trig::sin(x)
}
#[inline]
fn fpx_cos_f64(x: f64) -> f64 {
    aether_fpx::trig::cos(x)
}
#[inline]
fn fpx_tan_f64(x: f64) -> f64 {
    aether_fpx::trig::tan(x)
}
#[inline]
fn fpx_sinh_f64(x: f64) -> f64 {
    aether_fpx::trig::sinh(x)
}
#[inline]
fn fpx_cosh_f64(x: f64) -> f64 {
    aether_fpx::trig::cosh(x)
}
#[inline]
fn fpx_tanh_f64(x: f64) -> f64 {
    aether_fpx::trig::tanh(x)
}
#[inline]
fn fpx_asin_f64(x: f64) -> f64 {
    aether_fpx::trig::asin(x)
}
#[inline]
fn fpx_acos_f64(x: f64) -> f64 {
    aether_fpx::trig::acos(x)
}
#[inline]
fn fpx_atan_f64(x: f64) -> f64 {
    aether_fpx::trig::atan(x)
}
#[inline]
fn fpx_exp_f64(x: f64) -> f64 {
    aether_fpx::trig::exp(x)
}
#[inline]
fn fpx_log_f64(x: f64) -> f64 {
    aether_fpx::trig::log(x)
}
#[inline]
fn fpx_log1p_f64(x: f64) -> f64 {
    aether_fpx::trig::log1p(x)
}
#[inline]
fn fpx_expm1_f64(x: f64) -> f64 {
    aether_fpx::trig::expm1(x)
}

#[inline]
fn fpx_sin_f32(x: f32) -> f32 {
    aether_fpx::trig::sinf(x)
}
#[inline]
fn fpx_cos_f32(x: f32) -> f32 {
    aether_fpx::trig::cosf(x)
}
#[inline]
fn fpx_tan_f32(x: f32) -> f32 {
    aether_fpx::trig::tanf(x)
}
#[inline]
fn fpx_sinh_f32(x: f32) -> f32 {
    aether_fpx::trig::sinhf(x)
}
#[inline]
fn fpx_cosh_f32(x: f32) -> f32 {
    aether_fpx::trig::coshf(x)
}
#[inline]
fn fpx_tanh_f32(x: f32) -> f32 {
    aether_fpx::trig::tanhf(x)
}
#[inline]
fn fpx_asin_f32(x: f32) -> f32 {
    aether_fpx::trig::asinf(x)
}
#[inline]
fn fpx_acos_f32(x: f32) -> f32 {
    aether_fpx::trig::acosf(x)
}
#[inline]
fn fpx_atan_f32(x: f32) -> f32 {
    aether_fpx::trig::atanf(x)
}
#[inline]
fn fpx_exp_f32(x: f32) -> f32 {
    aether_fpx::trig::expf(x)
}
#[inline]
fn fpx_log_f32(x: f32) -> f32 {
    aether_fpx::trig::logf(x)
}
#[inline]
fn fpx_log1p_f32(x: f32) -> f32 {
    aether_fpx::trig::log1pf(x)
}
#[inline]
fn fpx_expm1_f32(x: f32) -> f32 {
    aether_fpx::trig::expm1f(x)
}

#[inline]
fn fpx_sin_f16(x: f16) -> f16 {
    aether_fpx::trig::sins(x)
}
#[inline]
fn fpx_cos_f16(x: f16) -> f16 {
    aether_fpx::trig::coss(x)
}
#[inline]
fn fpx_tan_f16(x: f16) -> f16 {
    aether_fpx::trig::tans(x)
}
#[inline]
fn fpx_sinh_f16(x: f16) -> f16 {
    aether_fpx::trig::sinhs(x)
}
#[inline]
fn fpx_cosh_f16(x: f16) -> f16 {
    aether_fpx::trig::coshs(x)
}
#[inline]
fn fpx_tanh_f16(x: f16) -> f16 {
    aether_fpx::trig::tanhs(x)
}
#[inline]
fn fpx_asin_f16(x: f16) -> f16 {
    aether_fpx::trig::asins(x)
}
#[inline]
fn fpx_acos_f16(x: f16) -> f16 {
    aether_fpx::trig::acoss(x)
}
#[inline]
fn fpx_atan_f16(x: f16) -> f16 {
    aether_fpx::trig::atans(x)
}
#[inline]
fn fpx_exp_f16(x: f16) -> f16 {
    aether_fpx::trig::exps(x)
}
#[inline]
fn fpx_log_f16(x: f16) -> f16 {
    aether_fpx::trig::logs(x)
}
#[inline]
fn fpx_log1p_f16(x: f16) -> f16 {
    aether_fpx::trig::log1ps(x)
}
#[inline]
fn fpx_expm1_f16(x: f16) -> f16 {
    aether_fpx::trig::expm1s(x)
}

#[inline]
fn fpx_sin_f128(x: f128) -> f128 {
    aether_fpx::trig::sinl(x)
}
#[inline]
fn fpx_cos_f128(x: f128) -> f128 {
    aether_fpx::trig::cosl(x)
}
#[inline]
fn fpx_tan_f128(x: f128) -> f128 {
    aether_fpx::trig::tanl(x)
}
#[inline]
fn fpx_sinh_f128(x: f128) -> f128 {
    aether_fpx::trig::sinhl(x)
}
#[inline]
fn fpx_cosh_f128(x: f128) -> f128 {
    aether_fpx::trig::coshl(x)
}
#[inline]
fn fpx_tanh_f128(x: f128) -> f128 {
    aether_fpx::trig::tanhl(x)
}
#[inline]
fn fpx_asin_f128(x: f128) -> f128 {
    aether_fpx::trig::asinl(x)
}
#[inline]
fn fpx_acos_f128(x: f128) -> f128 {
    aether_fpx::trig::acosl(x)
}
#[inline]
fn fpx_atan_f128(x: f128) -> f128 {
    aether_fpx::trig::atanl(x)
}
#[inline]
fn fpx_exp_f128(x: f128) -> f128 {
    aether_fpx::trig::expl(x)
}
#[inline]
fn fpx_log_f128(x: f128) -> f128 {
    aether_fpx::trig::logl(x)
}
#[inline]
fn fpx_log1p_f128(x: f128) -> f128 {
    aether_fpx::trig::log1pl(x)
}
#[inline]
fn fpx_expm1_f128(x: f128) -> f128 {
    aether_fpx::trig::expm1l(x)
}

/* ============================================================
 * main
 * ============================================================ */

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut only_arg: Option<String> = None;
    let mut output_path = std::path::PathBuf::from("target/precision_summary.csv");

    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "--only" => {
                if i + 1 >= args.len() {
                    eprintln!("missing value for --only (example: --only tanh,sinh,cosh)");
                    std::process::exit(2);
                }
                only_arg = Some(args[i + 1].clone());
                i += 2;
            }
            "--out" => {
                if i + 1 >= args.len() {
                    eprintln!("missing value for --out");
                    std::process::exit(2);
                }
                output_path = std::path::PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--help" | "-h" => {
                println!(
                    "Usage: precision [--only <func1,func2,...>] [--out <summary.csv>]\n\nIf --only is omitted, all precision tests run."
                );
                return;
            }
            other => {
                eprintln!("unknown argument: {other}");
                eprintln!("Use --help for usage.");
                std::process::exit(2);
            }
        }
    }

    let only_funcs: Option<std::collections::BTreeSet<String>> = only_arg.as_ref().map(|s| {
        s.split(',')
            .map(|v| v.trim().to_lowercase())
            .filter(|v| !v.is_empty())
            .collect()
    });

    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent).unwrap();
    }

    let mut summaries: Vec<PrecisionSummary> = Vec::new();

    fn append_reports(summaries: &mut Vec<PrecisionSummary>, reports: Vec<PrecisionReport>) {
        for report in reports {
            report.print();
            summaries.push(report.summary);
            println!();
        }
    }

    let funcs_f64 = [
        FuncF64 {
            name: "sin",
            reference: ref_sin_f64,
            fpx: fpx_sin_f64,
            lo: f64::MIN,
            hi: f64::MAX,
        },
        FuncF64 {
            name: "cos",
            reference: ref_cos_f64,
            fpx: fpx_cos_f64,
            lo: f64::MIN,
            hi: f64::MAX,
        },
        FuncF64 {
            name: "tan",
            reference: ref_tan_f64,
            fpx: fpx_tan_f64,
            lo: f64::MIN,
            hi: f64::MAX,
        },
        FuncF64 {
            name: "sinh",
            reference: ref_sinh_f64,
            fpx: fpx_sinh_f64,
            lo: f64::MIN,
            hi: f64::MAX,
        },
        FuncF64 {
            name: "cosh",
            reference: ref_cosh_f64,
            fpx: fpx_cosh_f64,
            lo: f64::MIN,
            hi: f64::MAX,
        },
        FuncF64 {
            name: "tanh",
            reference: ref_tanh_f64,
            fpx: fpx_tanh_f64,
            lo: f64::MIN,
            hi: f64::MAX,
        },
        FuncF64 {
            name: "asin",
            reference: ref_asin_f64,
            fpx: fpx_asin_f64,
            lo: -1.0,
            hi: 1.0,
        },
        FuncF64 {
            name: "acos",
            reference: ref_acos_f64,
            fpx: fpx_acos_f64,
            lo: -1.0,
            hi: 1.0,
        },
        FuncF64 {
            name: "atan",
            reference: ref_atan_f64,
            fpx: fpx_atan_f64,
            lo: f64::MIN,
            hi: f64::MAX,
        },
        FuncF64 {
            name: "exp",
            reference: ref_exp_f64,
            fpx: fpx_exp_f64,
            lo: f64::MIN,
            hi: f64::MAX,
        },
        FuncF64 {
            name: "log",
            reference: ref_log_f64,
            fpx: fpx_log_f64,
            lo: f64::from_bits(1),
            hi: f64::MAX,
        },
        FuncF64 {
            name: "log1p",
            reference: ref_log1p_f64,
            fpx: fpx_log1p_f64,
            lo: -1.0,
            hi: f64::MAX,
        },
        FuncF64 {
            name: "expm1",
            reference: ref_expm1_f64,
            fpx: fpx_expm1_f64,
            lo: f64::MIN,
            hi: f64::MAX,
        },
    ];
    let funcs_f64: Vec<FuncF64> = funcs_f64
        .into_iter()
        .filter(|f| {
            only_funcs
                .as_ref()
                .map(|set| set.contains(f.name))
                .unwrap_or(true)
        })
        .collect();

    let funcs_f32 = [
        FuncF32 {
            name: "sin",
            reference: ref_sin_f32,
            fpx: fpx_sin_f32,
            lo: f32::MIN,
            hi: f32::MAX,
        },
        FuncF32 {
            name: "cos",
            reference: ref_cos_f32,
            fpx: fpx_cos_f32,
            lo: f32::MIN,
            hi: f32::MAX,
        },
        FuncF32 {
            name: "tan",
            reference: ref_tan_f32,
            fpx: fpx_tan_f32,
            lo: f32::MIN,
            hi: f32::MAX,
        },
        FuncF32 {
            name: "sinh",
            reference: ref_sinh_f32,
            fpx: fpx_sinh_f32,
            lo: f32::MIN,
            hi: f32::MAX,
        },
        FuncF32 {
            name: "cosh",
            reference: ref_cosh_f32,
            fpx: fpx_cosh_f32,
            lo: f32::MIN,
            hi: f32::MAX,
        },
        FuncF32 {
            name: "tanh",
            reference: ref_tanh_f32,
            fpx: fpx_tanh_f32,
            lo: f32::MIN,
            hi: f32::MAX,
        },
        FuncF32 {
            name: "asin",
            reference: ref_asin_f32,
            fpx: fpx_asin_f32,
            lo: -1.0,
            hi: 1.0,
        },
        FuncF32 {
            name: "acos",
            reference: ref_acos_f32,
            fpx: fpx_acos_f32,
            lo: -1.0,
            hi: 1.0,
        },
        FuncF32 {
            name: "atan",
            reference: ref_atan_f32,
            fpx: fpx_atan_f32,
            lo: f32::MIN,
            hi: f32::MAX,
        },
        FuncF32 {
            name: "exp",
            reference: ref_exp_f32,
            fpx: fpx_exp_f32,
            lo: f32::MIN,
            hi: f32::MAX,
        },
        FuncF32 {
            name: "log",
            reference: ref_log_f32,
            fpx: fpx_log_f32,
            lo: f32::from_bits(1),
            hi: f32::MAX,
        },
        FuncF32 {
            name: "log1p",
            reference: ref_log1p_f32,
            fpx: fpx_log1p_f32,
            lo: -1.0,
            hi: f32::MAX,
        },
        FuncF32 {
            name: "expm1",
            reference: ref_expm1_f32,
            fpx: fpx_expm1_f32,
            lo: f32::MIN,
            hi: f32::MAX,
        },
    ];
    let funcs_f32: Vec<FuncF32> = funcs_f32
        .into_iter()
        .filter(|f| {
            only_funcs
                .as_ref()
                .map(|set| set.contains(f.name))
                .unwrap_or(true)
        })
        .collect();

    let funcs_f16 = [
        FuncF16 {
            name: "sin",
            reference: ref_sin_f16,
            fpx: fpx_sin_f16,
            lo: f16::MIN,
            hi: f16::MAX,
        },
        FuncF16 {
            name: "cos",
            reference: ref_cos_f16,
            fpx: fpx_cos_f16,
            lo: f16::MIN,
            hi: f16::MAX,
        },
        FuncF16 {
            name: "tan",
            reference: ref_tan_f16,
            fpx: fpx_tan_f16,
            lo: f16::MIN,
            hi: f16::MAX,
        },
        FuncF16 {
            name: "sinh",
            reference: ref_sinh_f16,
            fpx: fpx_sinh_f16,
            lo: f16::MIN,
            hi: f16::MAX,
        },
        FuncF16 {
            name: "cosh",
            reference: ref_cosh_f16,
            fpx: fpx_cosh_f16,
            lo: f16::MIN,
            hi: f16::MAX,
        },
        FuncF16 {
            name: "tanh",
            reference: ref_tanh_f16,
            fpx: fpx_tanh_f16,
            lo: f16::MIN,
            hi: f16::MAX,
        },
        FuncF16 {
            name: "asin",
            reference: ref_asin_f16,
            fpx: fpx_asin_f16,
            lo: -1.0,
            hi: 1.0,
        },
        FuncF16 {
            name: "acos",
            reference: ref_acos_f16,
            fpx: fpx_acos_f16,
            lo: -1.0,
            hi: 1.0,
        },
        FuncF16 {
            name: "atan",
            reference: ref_atan_f16,
            fpx: fpx_atan_f16,
            lo: f16::MIN,
            hi: f16::MAX,
        },
        FuncF16 {
            name: "exp",
            reference: ref_exp_f16,
            fpx: fpx_exp_f16,
            lo: f16::MIN,
            hi: f16::MAX,
        },
        FuncF16 {
            name: "log",
            reference: ref_log_f16,
            fpx: fpx_log_f16,
            lo: f16::from_bits(1),
            hi: f16::MAX,
        },
        FuncF16 {
            name: "log1p",
            reference: ref_log1p_f16,
            fpx: fpx_log1p_f16,
            lo: -1.0,
            hi: f16::MAX,
        },
        FuncF16 {
            name: "expm1",
            reference: ref_expm1_f16,
            fpx: fpx_expm1_f16,
            lo: f16::MIN,
            hi: f16::MAX,
        },
    ];
    let funcs_f16: Vec<FuncF16> = funcs_f16
        .into_iter()
        .filter(|f| {
            only_funcs
                .as_ref()
                .map(|set| set.contains(f.name))
                .unwrap_or(true)
        })
        .collect();

    let funcs_f128 = [
        FuncF128 {
            name: "sin",
            reference: ref_sin_f128,
            fpx: fpx_sin_f128,
            lo: f128::MIN,
            hi: f128::MAX,
        },
        FuncF128 {
            name: "cos",
            reference: ref_cos_f128,
            fpx: fpx_cos_f128,
            lo: f128::MIN,
            hi: f128::MAX,
        },
        FuncF128 {
            name: "tan",
            reference: ref_tan_f128,
            fpx: fpx_tan_f128,
            lo: f128::MIN,
            hi: f128::MAX,
        },
        FuncF128 {
            name: "sinh",
            reference: ref_sinh_f128,
            fpx: fpx_sinh_f128,
            lo: f128::MIN,
            hi: f128::MAX,
        },
        FuncF128 {
            name: "cosh",
            reference: ref_cosh_f128,
            fpx: fpx_cosh_f128,
            lo: f128::MIN,
            hi: f128::MAX,
        },
        FuncF128 {
            name: "tanh",
            reference: ref_tanh_f128,
            fpx: fpx_tanh_f128,
            lo: f128::MIN,
            hi: f128::MAX,
        },
        FuncF128 {
            name: "asin",
            reference: ref_asin_f128,
            fpx: fpx_asin_f128,
            lo: -1.0,
            hi: 1.0,
        },
        FuncF128 {
            name: "acos",
            reference: ref_acos_f128,
            fpx: fpx_acos_f128,
            lo: -1.0,
            hi: 1.0,
        },
        FuncF128 {
            name: "atan",
            reference: ref_atan_f128,
            fpx: fpx_atan_f128,
            lo: f128::MIN,
            hi: f128::MAX,
        },
        FuncF128 {
            name: "exp",
            reference: ref_exp_f128,
            fpx: fpx_exp_f128,
            lo: f128::MIN,
            hi: f128::MAX,
        },
        FuncF128 {
            name: "log",
            reference: ref_log_f128,
            fpx: fpx_log_f128,
            lo: f128::from_bits(1),
            hi: f128::MAX,
        },
        FuncF128 {
            name: "log1p",
            reference: ref_log1p_f128,
            fpx: fpx_log1p_f128,
            lo: -1.0,
            hi: f128::MAX,
        },
        FuncF128 {
            name: "expm1",
            reference: ref_expm1_f128,
            fpx: fpx_expm1_f128,
            lo: f128::MIN,
            hi: f128::MAX,
        },
    ];
    let funcs_f128: Vec<FuncF128> = funcs_f128
        .into_iter()
        .filter(|f| {
            only_funcs
                .as_ref()
                .map(|set| set.contains(f.name))
                .unwrap_or(true)
        })
        .collect();

    /* ================= f32 full-domain stratified bits ================= */

    println!("==== f32 full-domain stratified bits ====");
    let exponents32: Vec<u8> = (0..=254).collect();
    let mantissas32: [u32; 8] = [
        0,
        1,
        1_u32 << 22,
        (1_u32 << 23) - 1,
        0x0055_5555,
        0x002a_aaaa,
        0x007f_ffff,
        0x0000_0001,
    ];
    let reports: Vec<PrecisionReport> = funcs_f32
        .par_iter()
        .map(|f| {
            test_f32_stratified_bits(
                f.name,
                f.reference,
                f.fpx,
                f.lo,
                f.hi,
                &exponents32,
                &mantissas32,
            )
        })
        .collect();
    append_reports(&mut summaries, reports);

    /* ================= publication-interval regression sweeps ================= */

    println!("==== f64 publication intervals ====");
    let reports: Vec<PrecisionReport> = funcs_f64
        .par_iter()
        .copied()
        .map(test_f64_publication_interval)
        .collect();
    append_reports(&mut summaries, reports);

    println!("==== f128 publication intervals ====");
    let reports: Vec<PrecisionReport> = funcs_f128
        .par_iter()
        .copied()
        .map(test_f128_publication_interval)
        .collect();
    append_reports(&mut summaries, reports);

    /* ================= f16 exhaustive ================= */

    println!("==== f16 exhaustive ====");
    let reports: Vec<PrecisionReport> = funcs_f16
        .par_iter()
        .map(|f| test_f16_exhaustive(f.name, f.reference, f.fpx, f.lo, f.hi))
        .collect();
    append_reports(&mut summaries, reports);

    /* ================= f64 stratified bits ================= */

    println!("==== f64 full-domain stratified bits ====");
    let exponents: Vec<u16> = (0..=2046).collect();
    let mantissas: [u64; 8] = [
        0,
        1,
        (1u64 << 51),
        (1u64 << 52) - 1,
        0xAAAAAAAAAAAAA,
        0x5555555555555,
        0x000FFFFFFFFFFFFF,
        0x0000000000000001,
    ];

    let reports: Vec<PrecisionReport> = funcs_f64
        .par_iter()
        .map(|f| {
            test_f64_stratified_bits(
                f.name,
                f.reference,
                f.fpx,
                f.lo,
                f.hi,
                &exponents,
                &mantissas,
            )
        })
        .collect();
    append_reports(&mut summaries, reports);

    /* ================= f128 stratified bits ================= */

    println!("==== f128 full-domain stratified bits ====");
    let mut exponents128: Vec<u16> = (0..=32766).step_by(257).collect();
    exponents128.extend([1, 16382, 16383, 16384, 32765, 32766]);
    exponents128.sort_unstable();
    exponents128.dedup();
    let mantissas128: [u128; 8] = [
        0,
        1,
        1u128 << 111,
        (1u128 << 112) - 1,
        0xAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAu128 & ((1u128 << 112) - 1),
        0x55555555555555555555555555555555u128 & ((1u128 << 112) - 1),
        0x0000FFFFFFFFFFFFFFFFFFFFFFFFFFFFu128 & ((1u128 << 112) - 1),
        0x00000000000000000000000000000001u128,
    ];

    let reports: Vec<PrecisionReport> = funcs_f128
        .par_iter()
        .map(|f| {
            test_f128_stratified_bits(
                f.name,
                f.reference,
                f.fpx,
                f.lo,
                f.hi,
                &exponents128,
                &mantissas128,
            )
        })
        .collect();
    append_reports(&mut summaries, reports);

    /* ================= f64 PI refinement ================= */

    println!("==== f64 PI refinement ====");
    let reports: Vec<PrecisionReport> = funcs_f64
        .par_iter()
        .copied()
        .filter(|f| matches!(f.name, "sin" | "cos"))
        .map(|f| test_f64_pi_refinement(f.name, f.reference, f.fpx, 4096))
        .collect();
    append_reports(&mut summaries, reports);

    /* ================= f64 boundaries ================= */

    println!("==== f64 boundaries ====");
    let reports: Vec<PrecisionReport> = funcs_f64
        .par_iter()
        .copied()
        .map(|f| {
            let points = match f.name {
                "sin" | "cos" | "tan" => {
                    let mut pts = vec![f64::MIN, -PI, 0.0, PI / 6.0, PI / 4.0, PI / 3.0];
                    if f.name != "tan" {
                        pts.push(PI / 2.0);
                    }
                    pts.extend([PI, f64::MAX]);
                    pts
                }
                "asin" | "acos" => vec![-1.0, -0.5, 0.0, 0.5, 1.0],
                "atan" => vec![f64::MIN, -1.0, 0.0, 1.0, f64::MAX],
                "sinh" | "cosh" | "tanh" => {
                    vec![f64::MIN, -20.0, -1.0, 0.0, 1.0, 20.0, f64::MAX]
                }
                "exp" | "expm1" => {
                    vec![f64::MIN, -745.0, -1.0, 0.0, 1.0, 709.0, f64::MAX]
                }
                "log" => vec![
                    f64::from_bits(1),
                    f64::MIN_POSITIVE,
                    0.5,
                    1.0,
                    2.0,
                    f64::MAX,
                ],
                "log1p" => vec![
                    -1.0,
                    f64::from_bits((-1.0_f64).to_bits() - 1),
                    -0.5,
                    0.0,
                    1.0,
                    f64::MAX,
                ],
                _ => unreachable!(),
            };

            test_f64_boundaries(f.name, f.reference, f.fpx, &points)
        })
        .collect();
    append_reports(&mut summaries, reports);

    write_precision_csv(&output_path, &summaries).unwrap();
    let violations: Vec<_> = summaries
        .iter()
        .filter_map(|summary| {
            summary
                .max_ulp
                .parse::<u128>()
                .ok()
                .filter(|&ulp| ulp > 5)
                .map(|ulp| {
                    format!(
                        "{} {} {}: {} ULP",
                        summary.name, summary.precision, summary.test, ulp
                    )
                })
        })
        .collect();
    if !violations.is_empty() {
        panic!("5 ULP contract violated:\n{}", violations.join("\n"));
    }
    println!("\nWrote {}", output_path.display());
}
