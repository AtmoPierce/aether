// aether-fpx-gen/src/design/remez/create/generators.rs

use aether_core::real::real::RealCast;
use core::{f128, f16};
use aether_fpx_bench::ulp::{ulp_dist_f128, ulp_dist_f32, ulp_dist_f64};

use super::config::{
    candidate_meets_acceptance, effective_config, resolve_acceptance, FunctionSpec, SplitStrategy,
    SweepSpec,
};

use super::coefficients::{
    eval_coeff128_log_atanh, eval_coeff128_shape, eval_coeff32_log_atanh, eval_coeff32_shape,
    eval_coeff64_log_atanh, eval_coeff64_shape, poly_domain, trig_kernel, CoeffSet128, CoeffSet16,
    CoeffSet32, CoeffSet64, FunctionPlan, Nk, PolyDomain, SegmentPlan, SegmentShape, TrigKernel,
};

use crate::export::{
    write_coeffs_to_c, write_function_coeffs_to_json, write_function_gpu_kernels,
    write_function_to_rs, FunctionRemezExport,
};

use crate::remez::analysis::*;
use crate::remez::{max_error_on_grid, RemezApprox};

/* ================================================================================================
 * Entry points
 * ================================================================================================ */

pub fn generate_all_from_config() {
    generate_runtime_support();
    let mut cfg = effective_config();

    // Keep the manifest declarative while allowing surgical regeneration of
    // a changed generated module. Example: FPX_GEN_ONLY=sinh,cosh,tanh.
    if let Ok(only) = std::env::var("FPX_GEN_ONLY") {
        let selected: std::collections::BTreeSet<_> = only
            .split(',')
            .map(str::trim)
            .filter(|name| !name.is_empty())
            .collect();
        cfg.functions
            .retain(|name, _| selected.contains(name.as_str()));
        if cfg.functions.is_empty() {
            panic!("[aether-fpx-gen] FPX_GEN_ONLY selected no enabled functions");
        }
    }

    eprintln!("FPX config loaded: {} functions", cfg.functions.len());

    for (name, func) in &cfg.functions {
        println!("Generating Function: {}", name);

        eprintln!("[aether-fpx-gen] {}: building plan", name);
        let plan = build_function_plan(name, func, &cfg);

        // Rust
        eprintln!("[aether-fpx-gen] {}: writing Rust", name);
        write_function_to_rs(name, func, &plan);

        // JSON
        {
            eprintln!("[aether-fpx-gen] {}: writing JSON/GPU", name);
            let dom = plan.dom;
            let mut export = FunctionRemezExport::new(name.clone(), dom.table_domain_str().into());

            for seg in &plan.segments {
                export.add_segment_from_plan(seg);
            }

            write_function_gpu_kernels(&export);
            write_function_coeffs_to_json(&export);
        }

        // C
        eprintln!("[aether-fpx-gen] {}: writing C", name);
        generate_c_from_plan(&plan, &cfg);

        // Plots
        eprintln!("[aether-fpx-gen] {}: writing plots", name);
        generate_plots_from_plan(&plan, &cfg);
        eprintln!("[aether-fpx-gen] {}: done", name);
    }
    println!("Finished Generating!");
}

fn generate_runtime_support() {
    generate_payne_hanek_support();
    generate_bench_support();
}

fn generate_payne_hanek_support() {
    let manifest = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let source = manifest.join("templates/payne_hanek.rs");
    let destination = manifest.join("../../src/trig/payne_hanek.rs");
    let template = std::fs::read_to_string(&source).unwrap_or_else(|error| {
        panic!(
            "failed to read Payne-Hanek template {}: {error}",
            source.display()
        )
    });
    let generated = super::payne_hanek::render_source(&template);
    std::fs::write(&destination, generated)
        .unwrap_or_else(|error| panic!("failed to write {}: {error}", destination.display()));
    eprintln!("[aether-fpx-gen] wrote {}", destination.display());
}

fn generate_bench_support() {
    let manifest = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let source = manifest.join("templates/bench.rs");
    let destination = manifest.join("../../src/bench.rs");
    let generated = std::fs::read_to_string(&source).unwrap_or_else(|error| {
        panic!(
            "failed to read bench template {}: {error}",
            source.display()
        )
    });
    std::fs::write(&destination, generated)
        .unwrap_or_else(|error| panic!("failed to write {}: {error}", destination.display()));
    eprintln!("[aether-fpx-gen] wrote {}", destination.display());
}

/* ================================================================================================
 * Function plan construction
 * ================================================================================================ */

#[derive(Clone, Copy)]
struct IntervalX {
    a: f64,
    b: f64,
}

#[derive(Clone)]
struct FitResult {
    plan: Option<SegmentPlan>,
    worst_score: f64,
    worst_x: f64,
}

fn segment_shape(dom: PolyDomain, kernel: Option<TrigKernel>) -> SegmentShape {
    match dom {
        PolyDomain::X => SegmentShape::X,
        PolyDomain::Z => SegmentShape::Z(kernel.expect("Z-domain segment missing kernel")),
    }
}

fn build_segment_plan(
    dom: PolyDomain,
    kernel: Option<TrigKernel>,
    a: f128,
    b: f128,
    nk: Nk,
    coeffs: (CoeffSet128, CoeffSet64, CoeffSet32, CoeffSet16),
) -> SegmentPlan {
    let (c128, c64, c32, c16) = coeffs;
    SegmentPlan {
        a128: a,
        b128: b,
        a64: a as f64,
        b64: b as f64,
        a32: a as f32,
        b32: b as f32,
        c128,
        c64,
        c32,
        c16,
        nk,
        shape: segment_shape(dom, kernel),
    }
}

fn normalized_segment_score(
    name: &str,
    dom: PolyDomain,
    seg: &SegmentPlan,
    interval: IntervalX,
    samples: usize,
    acc: &super::config::AcceptanceResolved,
) -> (f64, f64) {
    let n = samples.max(2);
    let step = (interval.b - interval.a) / ((n - 1) as f64);
    let (f_runtime, _) = runtime_f_df(name);

    let mut worst_score = 0.0f64;
    let mut worst_x = interval.a;

    for i in 0..n {
        let x = interval.a + step * (i as f64);
        let y_ref128 = f_runtime(x as f128);
        let y_poly128 = if name == "log" {
            eval_coeff128_log_atanh(x as f128, &seg.c128)
        } else {
            eval_coeff128_shape(x as f128, &seg.c128, seg.shape)
        };

        let y_ref64 = y_ref128 as f64;
        let y_ref32 = y_ref128 as f32;
        let y_poly64 = if name == "log" {
            eval_coeff64_log_atanh(x, &seg.c64)
        } else {
            eval_coeff64_shape(x, &seg.c64, seg.shape)
        };
        let y_poly32 = if name == "log" {
            eval_coeff32_log_atanh(x as f32, &seg.c32)
        } else {
            eval_coeff32_shape(x as f32, &seg.c32, seg.shape)
        };

        let abs128 = (y_ref128 - y_poly128).abs();
        let abs64 = (y_ref64 - y_poly64).abs();
        let abs32 = (y_ref32 - y_poly32).abs();

        let ulp128 = ulp_dist_f128(y_ref128, y_poly128) as f64;
        let ulp64 = ulp_dist_f64(y_ref64, y_poly64) as f64;
        let ulp32 = ulp_dist_f32(y_ref32, y_poly32) as f64;

        let score = [
            (abs128 as f64) / ((acc.abs_tol_f128 as f64).max(f64::MIN_POSITIVE)),
            ulp128 / (acc.max_ulp_f128.max(1) as f64),
            abs64 / acc.abs_tol_f64.max(f64::MIN_POSITIVE),
            ulp64 / (acc.max_ulp_f64.max(1) as f64),
            (abs32 as f64) / (acc.abs_tol_f32.max(f32::MIN_POSITIVE) as f64),
            ulp32 / (acc.max_ulp_f32.max(1) as f64),
        ]
        .into_iter()
        .fold(0.0f64, f64::max);

        if score > worst_score {
            worst_score = score;
            worst_x = x;
        }
    }

    let _ = dom;
    (worst_score, worst_x)
}

fn try_fit_interval(
    name: &str,
    func: &FunctionSpec,
    cfg: &super::config::FpxRemezConfig,
    dom: PolyDomain,
    kernel: Option<TrigKernel>,
    sweep: &SweepSpec,
    interval: IntervalX,
) -> FitResult {
    let acc = resolve_acceptance(&cfg.defaults, func, Some(sweep), None);
    let samples = if sweep.samples != 0 {
        sweep.samples
    } else {
        cfg.defaults.samples
    };

    let mut best_fail: Option<(f64, f64)> = None;
    let trace_candidates = std::env::var_os("FPX_GEN_TRACE_CANDIDATES").is_some();

    for cand in &sweep.candidates {
        let nk = Nk {
            n: cand.n,
            k: cand.k,
        };
        if nk.k != nk.n + 1 {
            continue;
        }
        if trace_candidates {
            eprintln!(
                "[aether-fpx-gen] {name}: trying n={} k={} over [{:.17e}, {:.17e}]",
                nk.n, nk.k, interval.a, interval.b,
            );
        }

        let coeffs = match catch_unwind_quiet(|| {
            solve_segment_coeffs(
                name,
                func,
                cfg,
                dom,
                (interval.a as f128, interval.b as f128),
                nk,
            )
        }) {
            Ok(coeffs) => coeffs,
            Err(payload) => {
                eprintln!(
                    "[aether-fpx-gen] {name}: skipping singular n={} k={} candidate: {}",
                    nk.n,
                    nk.k,
                    panic_payload_to_string(payload),
                );
                continue;
            }
        };
        if trace_candidates {
            eprintln!("[aether-fpx-gen] {name}: solved n={} k={}", nk.n, nk.k);
        }

        if !coeffs_are_sane(name, &coeffs.0, &coeffs.1, &coeffs.2, &coeffs.3) {
            continue;
        }

        let seg = build_segment_plan(
            dom,
            kernel,
            interval.a as f128,
            interval.b as f128,
            nk,
            coeffs,
        );
        let (score, worst_x) = normalized_segment_score(name, dom, &seg, interval, samples, &acc);

        let passes = {
            candidate_meets_acceptance(
                name,
                dom,
                (interval.a, interval.b),
                &seg.c128,
                &seg.c64,
                &seg.c32,
                &acc,
                samples,
                cfg.defaults.samples,
            )
        };
        if trace_candidates {
            eprintln!(
                "[aether-fpx-gen] {name}: n={} k={} score={score:.6e} passes={passes}",
                nk.n, nk.k,
            );
        }

        if passes {
            // Candidate order is the configured cost/degree preference.  Stop
            // once that contract is met instead of solving larger candidates
            // merely to reduce an already-accepted error score.
            return FitResult {
                plan: Some(seg),
                worst_score: score,
                worst_x,
            };
        } else if best_fail.as_ref().map(|(s, _)| score < *s).unwrap_or(true) {
            best_fail = Some((score, worst_x));
        }
    }

    let (worst_score, worst_x) =
        best_fail.unwrap_or((f64::INFINITY, 0.5 * (interval.a + interval.b)));

    FitResult {
        plan: None,
        worst_score,
        worst_x,
    }
}

fn choose_split_point(interval: IntervalX, fit: &FitResult, sweep: &SweepSpec) -> f64 {
    let mid = 0.5 * (interval.a + interval.b);
    let lo = interval.a + sweep.min_width;
    let hi = interval.b - sweep.min_width;

    if !(lo < hi) {
        return mid;
    }

    let raw = match sweep.split_strategy {
        SplitStrategy::Mid => mid,
        SplitStrategy::Worst => fit.worst_x,
    };

    if raw.is_finite() {
        raw.clamp(lo, hi)
    } else {
        mid
    }
}

fn fail_sweep_acceptance(name: &str, sweep: &SweepSpec, fits: &[FitResult], limit: &str) -> ! {
    let best_overrun = fits
        .iter()
        .filter(|fit| fit.plan.is_none())
        .map(|fit| fit.worst_score)
        .fold(0.0_f64, f64::max);
    panic!(
        "[aether-fpx-gen] `{name}` failed configured acceptance after {limit} over [{:.17e}, {:.17e}]; best failing candidate was {:.3}x the configured error budget; refusing to export an over-budget approximation",
        sweep.interval[0].as_f64(),
        sweep.interval[1].as_f64(),
        best_overrun,
    );
}

fn build_segments_from_explicit(
    name: &str,
    func: &FunctionSpec,
    cfg: &super::config::FpxRemezConfig,
    dom: PolyDomain,
    kernel: Option<TrigKernel>,
) -> Vec<SegmentPlan> {
    let mut out = Vec::with_capacity(func.segments.len());
    for seg in &func.segments {
        let nk = Nk { n: seg.n, k: seg.k };
        let (a, b) = seg_interval_f128(seg);
        let coeffs = solve_segment_coeffs(name, func, cfg, dom, (a, b), nk);

        if !coeffs_are_sane(name, &coeffs.0, &coeffs.1, &coeffs.2, &coeffs.3) {
            panic!(
                "[aether-fpx-gen] explicit segment for `{}` over [{:.17e}, {:.17e}] produced non-finite or unstable coefficients; lower degree, widen the segment, or use sweep",
                name,
                a as f64,
                b as f64,
            );
        }

        let acc = resolve_acceptance(&cfg.defaults, func, None, Some(seg));
        let samples = seg
            .samples
            .filter(|&samples| samples != 0)
            .unwrap_or(cfg.defaults.samples)
            .max(2);
        if !candidate_meets_acceptance(
            name,
            dom,
            (a as f64, b as f64),
            &coeffs.0,
            &coeffs.1,
            &coeffs.2,
            &acc,
            samples,
            cfg.defaults.samples,
        ) {
            panic!(
                "[aether-fpx-gen] explicit segment for `{name}` over [{:.17e}, {:.17e}] failed configured acceptance after {samples} samples; refusing to export an over-budget approximation",
                a as f64,
                b as f64,
            );
        }

        out.push(build_segment_plan(dom, kernel, a, b, nk, coeffs));
    }
    out.sort_by(|lhs, rhs| lhs.a64.partial_cmp(&rhs.a64).unwrap());
    out
}

fn build_segments_from_sweep(
    name: &str,
    func: &FunctionSpec,
    cfg: &super::config::FpxRemezConfig,
    dom: PolyDomain,
    kernel: Option<TrigKernel>,
    sweep: &SweepSpec,
) -> Vec<SegmentPlan> {
    let mut intervals = vec![IntervalX {
        a: sweep.interval[0].as_f64(),
        b: sweep.interval[1].as_f64(),
    }];
    let mut fit_cache = std::collections::BTreeMap::<(u64, u64), FitResult>::new();

    loop {
        let fits: Vec<FitResult> = intervals
            .iter()
            .map(|iv| {
                let key = (iv.a.to_bits(), iv.b.to_bits());
                if let Some(fit) = fit_cache.get(&key) {
                    return fit.clone();
                }

                let fit = try_fit_interval(name, func, cfg, dom, kernel, sweep, *iv);
                fit_cache.insert(key, fit.clone());
                fit
            })
            .collect();

        if fits.iter().all(|fit| fit.plan.is_some()) {
            let mut segments: Vec<SegmentPlan> =
                fits.into_iter().map(|fit| fit.plan.unwrap()).collect();

            if sweep.merge_pass {
                loop {
                    let mut merged = false;
                    let mut next = Vec::new();
                    let mut i = 0usize;

                    while i < segments.len() {
                        if i + 1 < segments.len() {
                            let merged_iv = IntervalX {
                                a: segments[i].a64,
                                b: segments[i + 1].b64,
                            };
                            let fit =
                                try_fit_interval(name, func, cfg, dom, kernel, sweep, merged_iv);
                            if let Some(plan) = fit.plan {
                                next.push(plan);
                                merged = true;
                                i += 2;
                                continue;
                            }
                        }

                        next.push(segments[i].clone());
                        i += 1;
                    }

                    segments = next;
                    if !merged {
                        break;
                    }
                }
            }

            segments.sort_by(|lhs, rhs| lhs.a64.partial_cmp(&rhs.a64).unwrap());
            return segments;
        }

        if intervals.len() >= sweep.max_segments {
            fail_sweep_acceptance(
                name,
                sweep,
                &fits,
                &format!("reaching max_segments={}", sweep.max_segments),
            );
        }

        let mut split_idx = None;
        let mut split_score = f64::NEG_INFINITY;
        for (idx, (iv, fit)) in intervals.iter().zip(fits.iter()).enumerate() {
            if fit.plan.is_some() {
                continue;
            }
            if (iv.b - iv.a) < (2.0 * sweep.min_width) {
                continue;
            }
            if fit.worst_score > split_score {
                split_score = fit.worst_score;
                split_idx = Some(idx);
            }
        }

        let Some(idx) = split_idx else {
            fail_sweep_acceptance(
                name,
                sweep,
                &fits,
                &format!("reaching min_width={}", sweep.min_width),
            );
        };

        let iv = intervals[idx];
        let fit = &fits[idx];
        let x_split = choose_split_point(iv, fit, sweep);

        intervals.splice(
            idx..=idx,
            [
                IntervalX {
                    a: iv.a,
                    b: x_split,
                },
                IntervalX {
                    a: x_split,
                    b: iv.b,
                },
            ],
        );
    }
}

fn build_function_plan(
    name: &str,
    func: &FunctionSpec,
    cfg: &super::config::FpxRemezConfig,
) -> FunctionPlan {
    if matches!(name, "sinh" | "cosh" | "tanh") {
        return build_identity_function_plan(name, func);
    }

    let dom = poly_domain(name);
    let kernel = if dom == PolyDomain::Z {
        Some(trig_kernel(name))
    } else {
        None
    };

    let mut segments = if !func.segments.is_empty() {
        build_segments_from_explicit(name, func, cfg, dom, kernel)
    } else {
        let sweep = func.sweep.as_ref().unwrap_or_else(|| {
            panic!(
                "function `{}` has neither explicit segments nor sweep",
                name
            )
        });
        build_segments_from_sweep(name, func, cfg, dom, kernel, sweep)
    };

    segments.sort_by(|lhs, rhs| lhs.a64.partial_cmp(&rhs.a64).unwrap());
    validate_segment_transitions(name, dom, &segments);

    let first = segments.first().expect("empty segment plan");
    let last = segments.last().expect("empty segment plan");
    let a = first.a64;
    let b = last.b64;

    eprintln!("[aether-fpx-gen] {}: {} segment(s)", name, segments.len());
    for (idx, seg) in segments.iter().enumerate() {
        eprintln!(
            "  seg{}: [{:.17e}, {:.17e}] n={} k={}",
            idx, seg.a64, seg.b64, seg.nk.n, seg.nk.k
        );
    }

    FunctionPlan {
        name: name.to_string(),
        dom,
        kernel,
        interval_min_f128: first.a128,
        interval_max_f128: last.b128,
        interval_min_f64: a,
        interval_max_f64: b,
        interval_min_f32: a as f32,
        interval_max_f32: b as f32,
        segments,
    }
}

/// The generated public sinh/cosh/tanh implementations are identity based
/// (`expm1` away from zero and parity-constrained series near zero). They do
/// not evaluate a raw-x minimax table. Retain a minimal table only for legacy
/// low-precision backend plumbing, avoiding ill-conditioned degree-24 global
/// monomials and an unnecessary Remez solve during scalar regeneration.
fn build_identity_function_plan(name: &str, func: &FunctionSpec) -> FunctionPlan {
    let sweep = func
        .sweep
        .as_ref()
        .unwrap_or_else(|| panic!("identity-based function `{name}` has no sweep interval"));
    let a128 = sweep.interval[0].as_f64() as f128;
    let b128 = sweep.interval[1].as_f64() as f128;
    let segment = SegmentPlan {
        a128,
        b128,
        a64: a128 as f64,
        b64: b128 as f64,
        a32: a128 as f32,
        b32: b128 as f32,
        c128: CoeffSet128::N2([0.0_f128; 2]),
        c64: CoeffSet64::N2([0.0_f64; 2]),
        c32: CoeffSet32::N2([0.0_f32; 2]),
        c16: CoeffSet16::N2([0.0_f16; 2]),
        nk: Nk { n: 2, k: 3 },
        shape: SegmentShape::X,
    };

    eprintln!("[aether-fpx-gen] {name}: identity-based scalar plan");
    FunctionPlan {
        name: name.to_string(),
        dom: PolyDomain::X,
        kernel: None,
        segments: vec![segment],
        interval_min_f128: a128,
        interval_max_f128: b128,
        interval_min_f64: a128 as f64,
        interval_max_f64: b128 as f64,
        interval_min_f32: a128 as f32,
        interval_max_f32: b128 as f32,
    }
}

/// Diagnose a plan whose adjacent kernels disagree at their selector.
/// The checks use exactly the precision and evaluation order used by the
/// generated polynomial helpers.
fn validate_segment_transitions(name: &str, dom: PolyDomain, segments: &[SegmentPlan]) {
    if segments.len() < 2 {
        return;
    }

    for pair in segments.windows(2) {
        let left = &pair[0];
        let right = &pair[1];
        if left.b128 != right.a128 {
            panic!("[aether-fpx-gen] `{name}` has non-contiguous segment geometry");
        }

        let x128 = left.b128;
        let y128_left = if name == "log" {
            eval_coeff128_log_atanh(x128, &left.c128)
        } else {
            eval_coeff128_shape(x128, &left.c128, left.shape)
        };
        let y128_right = if name == "log" {
            eval_coeff128_log_atanh(x128, &right.c128)
        } else {
            eval_coeff128_shape(x128, &right.c128, right.shape)
        };

        let x64 = left.b64;
        let y64_left = if name == "log" {
            eval_coeff64_log_atanh(x64, &left.c64)
        } else {
            eval_coeff64_shape(x64, &left.c64, left.shape)
        };
        let y64_right = if name == "log" {
            eval_coeff64_log_atanh(x64, &right.c64)
        } else {
            eval_coeff64_shape(x64, &right.c64, right.shape)
        };

        let x32 = left.b32;
        let y32_left = if name == "log" {
            eval_coeff32_log_atanh(x32, &left.c32)
        } else {
            eval_coeff32_shape(x32, &left.c32, left.shape)
        };
        let y32_right = if name == "log" {
            eval_coeff32_log_atanh(x32, &right.c32)
        } else {
            eval_coeff32_shape(x32, &right.c32, right.shape)
        };

        let scale128 = y128_left.abs().max(y128_right.abs()).max(1.0_f128);
        let scale64 = y64_left.abs().max(y64_right.abs()).max(1.0_f64);
        let scale32 = y32_left.abs().max(y32_right.abs()).max(1.0_f32);
        let scalar_uses_f128_segment = !matches!(name, "asin" | "atan");
        let scalar_uses_f64_segment = name != "asin" || x64 > 1.25e-1_f64;
        let scalar_uses_f32_segment = name != "asin" || x32 > 1.25e-1_f32;
        if scalar_uses_f128_segment
            && (y128_left - y128_right).abs() > 64.0_f128 * f128::EPSILON * scale128
            || scalar_uses_f64_segment
                && (y64_left - y64_right).abs() > 64.0_f64 * f64::EPSILON * scale64
            || scalar_uses_f32_segment
                && (y32_left - y32_right).abs() > 64.0_f32 * f32::EPSILON * scale32
        {
            eprintln!(
                "[aether-fpx-gen] warning: experimental `{name}` kernel fails segment continuity validation at {:?}",
                x128
            );
        }
    }

    let _ = dom;
}

/* ================================================================================================
 * Coefficient solving (single source of truth)
 * ================================================================================================ */

pub fn solve_segment_coeffs(
    name: &str,
    func: &FunctionSpec,
    cfg: &super::config::FpxRemezConfig,
    dom: PolyDomain,
    interval: (f128, f128),
    nk: Nk,
) -> (CoeffSet128, CoeffSet64, CoeffSet32, CoeffSet16) {
    let _ = resolve_acceptance(&cfg.defaults, func, func.sweep.as_ref(), None);
    if nk.k != nk.n + 1 {
        panic!("solver requires k = n + 1 for {}", name);
    }

    let interval_x = interval;
    let interval = solve_interval_from_x_for_name(name, dom, interval.0, interval.1);
    let dom_name = match dom {
        PolyDomain::X => "x",
        PolyDomain::Z => "z",
    };

    match catch_unwind_quiet(|| match nk.n {
        2 => coeffs_from_remez::<2, 3>(dom, name, interval),
        3 => coeffs_from_remez::<3, 4>(dom, name, interval),
        4 => coeffs_from_remez::<4, 5>(dom, name, interval),
        5 => coeffs_from_remez::<5, 6>(dom, name, interval),
        6 => coeffs_from_remez::<6, 7>(dom, name, interval),
        7 => coeffs_from_remez::<7, 8>(dom, name, interval),
        8 => coeffs_from_remez::<8, 9>(dom, name, interval),
        9 => coeffs_from_remez::<9, 10>(dom, name, interval),
        10 => coeffs_from_remez::<10, 11>(dom, name, interval),
        11 => coeffs_from_remez::<11, 12>(dom, name, interval),
        12 => coeffs_from_remez::<12, 13>(dom, name, interval),
        13 => coeffs_from_remez::<13, 14>(dom, name, interval),
        14 => coeffs_from_remez::<14, 15>(dom, name, interval),
        15 => coeffs_from_remez::<15, 16>(dom, name, interval),
        16 => coeffs_from_remez::<16, 17>(dom, name, interval),
        17 => coeffs_from_remez::<17, 18>(dom, name, interval),
        18 => coeffs_from_remez::<18, 19>(dom, name, interval),
        19 => coeffs_from_remez::<19, 20>(dom, name, interval),
        20 => coeffs_from_remez::<20, 21>(dom, name, interval),
        21 => coeffs_from_remez::<21, 22>(dom, name, interval),
        22 => coeffs_from_remez::<22, 23>(dom, name, interval),
        23 => coeffs_from_remez::<23, 24>(dom, name, interval),
        24 => coeffs_from_remez::<24, 25>(dom, name, interval),
        25 => coeffs_from_remez::<25, 26>(dom, name, interval),
        n => panic!("unsupported n={n}"),
    }) {
        Ok(coeffs) => coeffs,
        Err(payload) => {
            let cause = panic_payload_to_string(payload);
            panic!(
                "[aether-fpx-gen] `{}` remez solve failed for n={} k={} on {}-domain x-interval [{:.17e}, {:.17e}] mapped to [{:.17e}, {:.17e}]: {}",
                name,
                nk.n,
                nk.k,
                dom_name,
                interval_x.0 as f64,
                interval_x.1 as f64,
                interval.0 as f64,
                interval.1 as f64,
                cause,
            );
        }
    }
}

fn catch_unwind_quiet<F, R>(f: F) -> Result<R, Box<dyn std::any::Any + Send>>
where
    F: FnOnce() -> R,
{
    let hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(f));
    std::panic::set_hook(hook);
    result
}

fn panic_payload_to_string(payload: Box<dyn std::any::Any + Send>) -> String {
    match payload.downcast::<String>() {
        Ok(message) => *message,
        Err(payload) => match payload.downcast::<&'static str>() {
            Ok(message) => (*message).to_string(),
            Err(_) => "non-string panic payload".to_string(),
        },
    }
}

#[inline]
fn checked_log_f16_coeff(name: &str, i: usize, c: f128) -> f16 {
    let out = c as f16;

    if !out.is_finite() {
        panic!(
            "[aether-fpx-gen] f16 coefficient overflow for `{}`: c[{}] = {:?}; \
             log f16 cannot store this monomial coefficient. Use lower degree, split the interval, \
             or evaluate the f32 polynomial for the f16 public function.",
            name, i, c,
        );
    }

    out
}

fn coeffs_from_remez<const N: usize, const K: usize>(
    dom: PolyDomain,
    name: &str,
    interval: (f128, f128),
) -> (CoeffSet128, CoeffSet64, CoeffSet32, CoeffSet16) {
    let rz128 = solve_remez::<f128, N, K>(dom, name, (interval.0, interval.1));
    // The log tail lives on z < 0.03 and binary128 accuracy requires a much
    // higher degree than binary64.  Solving the same high-degree monomial
    // Vandermonde again in f64 is ill-conditioned; the resulting coefficients
    // can overflow even though the binary128 Remez solution is well behaved.
    // Cast the accepted high-precision design for log, while retaining the
    // existing precision-specific solve for every other function.
    let rz64 = (name != "log")
        .then(|| solve_remez::<f64, N, K>(dom, name, (interval.0 as f64, interval.1 as f64)));

    let mut a128 = [0.0f128; N];
    let mut a64 = [0.0f64; N];
    let mut a32 = [0.0f32; N];
    let mut a16 = [0.0f32 as f16; N];

    for i in 0..N {
        a128[i] = rz128.coefficients[i];
        let lower = rz64
            .as_ref()
            .map(|rz| rz.coefficients[i])
            .unwrap_or(rz128.coefficients[i] as f64);
        a64[i] = lower;
        a32[i] = lower as f32;
        a16[i] = lower as f16;
    }
    (
        CoeffSet128::from_slice(&a128),
        CoeffSet64::from_slice(&a64),
        CoeffSet32::from_slice(&a32),
        CoeffSet16::from_slice(&a16),
    )
}

/* ================================================================================================
 * C emission (FROM PLAN)
 * ================================================================================================ */

fn generate_c_from_plan(plan: &FunctionPlan, cfg: &super::config::FpxRemezConfig) {
    let name = plan.name.as_str();
    let dom = plan.dom;
    let samples = cfg.defaults.samples;

    let (f_solve, _) = solve_f_df(dom, name);

    for (i, seg) in plan.segments.iter().enumerate() {
        let interval_x = (seg.a128, seg.b128);
        let solve_interval = solve_interval_from_x_for_name(name, dom, interval_x.0, interval_x.1);

        macro_rules! go {
            ($N:expr, $K:expr) => {{
                let mut rz = RemezApprox::<f128, $N, $K>::new(solve_interval.0, solve_interval.1);

                for j in 0..$N {
                    rz.coefficients[j] = f128::from_f64(seg.c64.as_slice()[j]);
                }

                let max_err =
                    max_error_on_grid(&rz, f_solve, solve_interval.0, solve_interval.1, samples);

                write_coeffs_to_c(
                    &rz,
                    name,
                    f_solve,
                    interval_x,
                    max_err,
                    &format!("{}_seg{}", name, i),
                );
            }};
        }

        match seg.nk.n {
            2 => go!(2, 3),
            3 => go!(3, 4),
            4 => go!(4, 5),
            5 => go!(5, 6),
            6 => go!(6, 7),
            7 => go!(7, 8),
            8 => go!(8, 9),
            9 => go!(9, 10),
            10 => go!(10, 11),
            11 => go!(11, 12),
            12 => go!(12, 13),
            13 => go!(13, 14),
            14 => go!(14, 15),
            15 => go!(15, 16),
            16 => go!(16, 17),
            17 => go!(17, 18),
            18 => go!(18, 19),
            19 => go!(19, 20),
            20 => go!(20, 21),
            21 => go!(21, 22),
            22 => go!(22, 23),
            23 => go!(23, 24),
            24 => go!(24, 25),
            25 => go!(25, 26),
            n => panic!("unsupported n={n}"),
        }
    }
}

/* ================================================================================================
 * Plots (FROM PLAN)
 * ================================================================================================ */

fn generate_plots_from_plan(plan: &FunctionPlan, cfg: &super::config::FpxRemezConfig) {
    let samples = cfg.defaults.samples;
    let name = plan.name.as_str();

    let (f_runtime, _) = runtime_f_df(name);
    let (f_solve, _) = solve_f_df(plan.dom, name);

    for (i, seg) in plan.segments.iter().enumerate() {
        let file_prefix = format!("{}_seg{}", name, i);

        let kernel = match seg.shape {
            SegmentShape::X => None,
            SegmentShape::Z(k) => Some(k),
        };

        let interval_x = (seg.a128, seg.b128);
        let solve_interval =
            solve_interval_from_x_for_name(name, plan.dom, interval_x.0, interval_x.1);

        macro_rules! go {
            ($N:expr, $K:expr) => {{
                let mut rz = RemezApprox::<f128, $N, $K>::new(solve_interval.0, solve_interval.1);

                for j in 0..$N {
                    rz.coefficients[j] = f128::from_f64(seg.c64.as_slice()[j]);
                }

                crate::plots::write_segment_plots(
                    name,
                    &file_prefix,
                    plan.dom,
                    kernel,
                    &rz,
                    interval_x,
                    solve_interval,
                    samples,
                    f_runtime,
                    f_solve,
                );
            }};
        }

        match seg.nk.n {
            2 => go!(2, 3),
            3 => go!(3, 4),
            4 => go!(4, 5),
            5 => go!(5, 6),
            6 => go!(6, 7),
            7 => go!(7, 8),
            8 => go!(8, 9),
            9 => go!(9, 10),
            10 => go!(10, 11),
            11 => go!(11, 12),
            12 => go!(12, 13),
            13 => go!(13, 14),
            14 => go!(14, 15),
            15 => go!(15, 16),
            16 => go!(16, 17),
            17 => go!(17, 18),
            18 => go!(18, 19),
            19 => go!(19, 20),
            20 => go!(20, 21),
            21 => go!(21, 22),
            22 => go!(22, 23),
            23 => go!(23, 24),
            24 => go!(24, 25),
            25 => go!(25, 26),
            n => panic!("unsupported n={n}"),
        }
    }
}

fn coeffs_are_sane(
    name: &str,
    c128: &CoeffSet128,
    c64: &CoeffSet64,
    c32: &CoeffSet32,
    c16: &CoeffSet16,
) -> bool {
    for &c in c128.as_slice() {
        if !c.is_finite() {
            return false;
        }
    }

    for &c in c64.as_slice() {
        if !c.is_finite() {
            return false;
        }

        // Tan residual monomial fits with huge coefficients are usually bad.
        if name == "tan" && c.abs() > 1.0e12 {
            return false;
        }
    }

    for &c in c32.as_slice() {
        if !c.is_finite() {
            return false;
        }

        if name == "tan" && c.abs() > 1.0e8 {
            return false;
        }
    }

    for &c in c16.as_slice() {
        if !c.is_finite() {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::super::config::{AcceptanceSpec, NkCandidate, SegmentSpec, TomlNumber};
    use super::*;

    fn permissive_acceptance() -> AcceptanceSpec {
        AcceptanceSpec {
            max_ulp_f128: Some(u128::MAX),
            max_ulp_f64: Some(u64::MAX),
            max_ulp_f32: Some(u32::MAX),
            abs_tol_f128: Some(1.0),
            abs_tol_f64: Some(1.0),
            abs_tol_f32: Some(1.0),
        }
    }

    fn exact_acceptance() -> AcceptanceSpec {
        AcceptanceSpec {
            max_ulp_f128: Some(0),
            max_ulp_f64: Some(0),
            max_ulp_f32: Some(0),
            abs_tol_f128: Some(0.0),
            abs_tol_f64: Some(0.0),
            abs_tol_f32: Some(0.0),
        }
    }

    fn explicit_exp_function(
        function_acceptance: AcceptanceSpec,
        segment_acceptance: AcceptanceSpec,
        samples: usize,
    ) -> FunctionSpec {
        FunctionSpec {
            file_prefix: "exp".into(),
            quadrants: Vec::new(),
            mirror: false,
            segments: vec![SegmentSpec {
                seg: Some("test".into()),
                n: 2,
                k: 3,
                interval: [TomlNumber::Num(0.0), TomlNumber::Num(0.01)],
                samples: Some(samples),
                acceptance: Some(segment_acceptance),
            }],
            sweep: None,
            acceptance: Some(function_acceptance),
        }
    }

    #[test]
    fn sweep_selects_the_first_candidate_that_meets_acceptance() {
        let sweep = SweepSpec {
            enabled: true,
            interval: [TomlNumber::Num(0.0), TomlNumber::Num(0.01)],
            min_width: 0.001,
            max_segments: 1,
            split_strategy: SplitStrategy::Worst,
            merge_pass: false,
            samples: 17,
            candidates: vec![NkCandidate { n: 2, k: 3 }, NkCandidate { n: 5, k: 6 }],
            w_madds: 1.0,
            w_segments: 3.0,
            w_selectors: 1.0,
            acceptance: None,
        };
        let function = FunctionSpec {
            file_prefix: "exp".into(),
            quadrants: Vec::new(),
            mirror: false,
            segments: Vec::new(),
            sweep: Some(sweep.clone()),
            acceptance: Some(AcceptanceSpec {
                max_ulp_f128: Some(u128::MAX),
                max_ulp_f64: Some(u64::MAX),
                max_ulp_f32: Some(u32::MAX),
                abs_tol_f128: Some(1.0),
                abs_tol_f64: Some(1.0),
                abs_tol_f32: Some(1.0),
            }),
        };

        let config = super::super::config::FpxRemezConfig::default();
        let segments =
            build_segments_from_sweep("exp", &function, &config, PolyDomain::X, None, &sweep);

        assert_eq!(segments.len(), 1);
        assert_eq!(segments[0].nk.n, 2);
    }

    #[test]
    #[should_panic(expected = "refusing to export an over-budget approximation")]
    fn exhausted_sweep_never_returns_a_failing_fallback() {
        let sweep = SweepSpec {
            enabled: true,
            interval: [TomlNumber::Num(0.0), TomlNumber::Num(1.0)],
            min_width: 0.1,
            max_segments: 1,
            split_strategy: SplitStrategy::Worst,
            merge_pass: false,
            samples: 17,
            candidates: vec![NkCandidate { n: 2, k: 3 }],
            w_madds: 1.0,
            w_segments: 3.0,
            w_selectors: 1.0,
            acceptance: None,
        };
        let function = FunctionSpec {
            file_prefix: "exp".into(),
            quadrants: Vec::new(),
            mirror: false,
            segments: Vec::new(),
            sweep: Some(sweep.clone()),
            acceptance: Some(AcceptanceSpec {
                max_ulp_f128: Some(1),
                max_ulp_f64: Some(1),
                max_ulp_f32: Some(1),
                abs_tol_f128: Some(1.0),
                abs_tol_f64: Some(1.0),
                abs_tol_f32: Some(1.0),
            }),
        };

        let config = super::super::config::FpxRemezConfig::default();
        let _ = build_segments_from_sweep("exp", &function, &config, PolyDomain::X, None, &sweep);
    }

    #[test]
    fn explicit_segment_uses_its_acceptance_override() {
        let function = explicit_exp_function(exact_acceptance(), permissive_acceptance(), 7);
        let config = super::super::config::FpxRemezConfig::default();

        let segments = build_segments_from_explicit("exp", &function, &config, PolyDomain::X, None);

        assert_eq!(segments.len(), 1);
        assert_eq!(segments[0].nk.n, 2);
    }

    #[test]
    #[should_panic(expected = "failed configured acceptance after 7 samples")]
    fn explicit_segment_refuses_an_over_budget_approximation() {
        let function = explicit_exp_function(permissive_acceptance(), exact_acceptance(), 7);
        let config = super::super::config::FpxRemezConfig::default();

        let _ = build_segments_from_explicit("exp", &function, &config, PolyDomain::X, None);
    }

    #[test]
    fn degree_dispatch_preserves_n11_and_n13_coefficient_counts() {
        let function = FunctionSpec {
            file_prefix: "exp".into(),
            quadrants: Vec::new(),
            mirror: false,
            segments: Vec::new(),
            sweep: None,
            acceptance: None,
        };
        let config = super::super::config::FpxRemezConfig::default();

        for n in [11, 13] {
            let (c128, c64, c32, c16) = solve_segment_coeffs(
                "exp",
                &function,
                &config,
                PolyDomain::X,
                (-1.0_f128, 1.0_f128),
                Nk { n, k: n + 1 },
            );

            assert_eq!(c128.as_slice().len(), n);
            assert_eq!(c64.as_slice().len(), n);
            assert_eq!(c32.as_slice().len(), n);
            assert_eq!(c16.as_slice().len(), n);
        }
    }
}
