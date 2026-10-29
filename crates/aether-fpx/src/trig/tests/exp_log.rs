use crate::trig::{exp, expf, expm1, expm1f, log, log1pf, log1pl, logf, logl, logs};

use super::support::{ulp_diff_f128, ulp_diff_f32, ulp_diff_f64};

#[test]
fn expf_is_stable_near_zero_and_ln2_boundaries() {
    let points = [
        -0.25_f32,
        -1.0e-6_f32,
        0.0_f32,
        1.0e-6_f32,
        0.25_f32,
        core::f32::consts::LN_2 - 1.0e-5_f32,
        core::f32::consts::LN_2,
        core::f32::consts::LN_2 + 1.0e-5_f32,
    ];

    for x in points {
        let got = expf(x);
        let want = x.exp();
        assert!(
            ulp_diff_f32(got, want) <= 4,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn exp_is_stable_near_zero_and_ln2_boundaries() {
    let points = [
        -0.25_f64,
        -1.0e-12_f64,
        0.0_f64,
        1.0e-12_f64,
        0.25_f64,
        core::f64::consts::LN_2 - 1.0e-12_f64,
        core::f64::consts::LN_2,
        core::f64::consts::LN_2 + 1.0e-12_f64,
    ];

    for x in points {
        let got = exp(x);
        let want = x.exp();
        assert!(
            ulp_diff_f64(got, want) <= 5,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn logf_is_stable_near_one() {
    let points = [
        f32::from_bits(1.0_f32.to_bits() - 8),
        f32::from_bits(1.0_f32.to_bits() - 1),
        1.0_f32,
        f32::from_bits(1.0_f32.to_bits() + 1),
        f32::from_bits(1.0_f32.to_bits() + 8),
        0.75_f32,
        1.25_f32,
    ];

    for x in points {
        let got = logf(x);
        let want = x.ln();
        assert!(
            ulp_diff_f32(got, want) <= 4,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn logf_generated_kernel_is_within_five_ulp_for_small_offsets() {
    let points = [
        1.0_f32 - 9.0e-7_f32,
        1.0_f32 + 9.0e-7_f32,
        1.0_f32 - 5.0e-5_f32,
        1.0_f32 + 5.0e-5_f32,
    ];

    for x in points {
        let got = logf(x);
        let want = x.ln();
        let ulp = ulp_diff_f32(got, want);
        assert!(ulp <= 5, "x={x:?} got={got:?} want={want:?} ulp={ulp}");
    }
}

#[test]
fn log_is_stable_near_one() {
    let points = [
        f64::from_bits(1.0_f64.to_bits() - 16),
        f64::from_bits(1.0_f64.to_bits() - 1),
        1.0_f64,
        f64::from_bits(1.0_f64.to_bits() + 1),
        f64::from_bits(1.0_f64.to_bits() + 16),
        0.75_f64,
        1.25_f64,
    ];

    for x in points {
        let got = log(x);
        let want = x.ln();
        assert!(
            ulp_diff_f64(got, want) <= 5,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn generated_log_remez_path_stays_within_five_ulp() {
    const SAMPLES: usize = 4097;

    for index in 0..SAMPLES {
        let t64 = index as f64 / (SAMPLES - 1) as f64;
        let x64 = core::f64::consts::FRAC_1_SQRT_2
            + (core::f64::consts::SQRT_2 - core::f64::consts::FRAC_1_SQRT_2) * t64;
        let got64 = log(x64);
        let want64 = x64.ln();
        let ulp64 = ulp_diff_f64(got64, want64);
        assert!(
            ulp64 <= 5,
            "f64 generated log kernel exceeded 5 ULP at x={x64:?}: got={got64:?}, want={want64:?}, ulp={ulp64}"
        );

        let t128 = (index as f128) / ((SAMPLES - 1) as f128);
        let x128 = core::f128::consts::FRAC_1_SQRT_2
            + (core::f128::consts::SQRT_2 - core::f128::consts::FRAC_1_SQRT_2) * t128;
        let got128 = logl(x128);
        let want128 = x128.ln();
        let ulp128 = ulp_diff_f128(got128, want128);
        assert!(
            ulp128 <= 5,
            "f128 generated log kernel exceeded 5 ULP at x={x128:?}: got={got128:?}, want={want128:?}, ulp={ulp128}"
        );
    }

    // This was the worst observed input for the former unconstrained fit.
    let former_failure = 1.000_976_093_018_660_6_f64;
    assert!(ulp_diff_f64(log(former_failure), former_failure.ln()) <= 5);

    // Exercise the first values outside logl's near-one log1p bypass.
    let transition = 0.000_244_140_625_f128;
    for boundary in [1.0_f128 - transition, 1.0_f128 + transition] {
        for x in [
            f128::from_bits(boundary.to_bits() - 1),
            boundary,
            f128::from_bits(boundary.to_bits() + 1),
        ] {
            let got = logl(x);
            let want = x.ln();
            let ulp = ulp_diff_f128(got, want);
            assert!(
                ulp <= 5,
                "f128 log transition exceeded 5 ULP at x={x:?}: got={got:?}, want={want:?}, ulp={ulp}"
            );
        }
    }

    for x in [
        f128::from_bits(1.0_f128.to_bits() - 1),
        f128::from_bits(1.0_f128.to_bits() + 1),
    ] {
        let got = logl(x);
        let want = x.ln();
        let ulp = ulp_diff_f128(got, want);
        assert!(
            ulp <= 5,
            "f128 log neighbor exceeded 5 ULP at x={x:?}: got={got:?}, want={want:?}, ulp={ulp}"
        );
    }
}

#[test]
fn logs_uses_the_generated_kernel_around_one() {
    for x in [
        f16::from_bits(1.0_f16.to_bits() - 1),
        f16::from_bits(1.0_f16.to_bits() + 1),
    ] {
        assert_eq!(logs(x), (x as f32).ln() as f16, "x={x:?}");
    }
}

#[test]
fn logl_is_exact_at_one() {
    assert_eq!(logl(1.0_f128), 0.0_f128);
}

#[test]
fn logl_matches_core_near_one_for_small_f128_offsets() {
    let points = [
        1.0_f128 - 9.0e-7_f128,
        1.0_f128 + 9.0e-7_f128,
        1.0_f128 - 5.0e-5_f128,
        1.0_f128 + 5.0e-5_f128,
    ];

    for x in points {
        let got = logl(x);
        let want = x.ln();
        assert!(
            ulp_diff_f128(got, want) <= 2,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn log1pl_is_within_five_ulp_at_f128_epsilon() {
    assert!(ulp_diff_f128(log1pl(f128::EPSILON), f128::EPSILON) <= 5);
}

#[test]
fn log1pl_matches_core_near_zero_for_small_f128_offsets() {
    let points = [-9.0e-7_f128, 9.0e-7_f128, -5.0e-5_f128, 5.0e-5_f128];

    for x in points {
        let got = log1pl(x);
        let want = x.ln_1p();
        assert!(
            ulp_diff_f128(got, want) <= 2,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn log1pf_matches_core_near_zero_for_small_f32_offsets() {
    let points = [-9.0e-7_f32, 9.0e-7_f32, -5.0e-5_f32, 5.0e-5_f32];

    for x in points {
        assert_eq!(log1pf(x), x.ln_1p(), "x={x:?}");
    }
}

#[test]
fn expm1f_matches_core_near_zero_for_small_f32_offsets() {
    let points = [-9.0e-7_f32, 9.0e-7_f32, -5.0e-5_f32, 5.0e-5_f32];

    for x in points {
        let got = expm1f(x);
        let want = x.exp_m1();
        assert!(
            ulp_diff_f32(got, want) <= 2,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn expm1_matches_core_near_zero_for_small_f64_offsets() {
    let points = [
        -1.0e-12_f64,
        1.0e-12_f64,
        -5.0e-5_f64,
        5.0e-5_f64,
        -3.5e-3_f64,
        3.5e-3_f64,
    ];

    for x in points {
        let got = expm1(x);
        let want = x.exp_m1();
        assert!(
            ulp_diff_f64(got, want) <= 2,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn expm1f_matches_core_through_transition_band() {
    let points = [-3.5e-3_f32, -2.5e-3_f32, 2.5e-3_f32, 3.5e-3_f32];

    for x in points {
        let got = expm1f(x);
        let want = x.exp_m1();
        assert!(
            ulp_diff_f32(got, want) <= 2,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}
