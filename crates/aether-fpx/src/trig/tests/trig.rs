use crate::trig::{
    acos, acosf, asin, asinf, cos, cosf, coshl, sin, sinf, sinh, sinhf, tan, tanf, tanh, tanhf,
};

use super::support::{ulp_diff_f128, ulp_diff_f32, ulp_diff_f64};

fn assert_abs_close_f32(got: f32, want: f32, abs_tol: f32, context: &str) {
    assert!(
        (got - want).abs() <= abs_tol,
        "{context} got={got:?} want={want:?}"
    );
}

fn assert_abs_close_f64(got: f64, want: f64, abs_tol: f64, context: &str) {
    assert!(
        (got - want).abs() <= abs_tol,
        "{context} got={got:?} want={want:?}"
    );
}

#[test]
fn sinf_is_stable_across_quadrants_and_reduction_points() {
    let points = [
        -1.0e-6_f32,
        0.0_f32,
        1.0e-6_f32,
        core::f32::consts::FRAC_PI_6,
        core::f32::consts::FRAC_PI_2 - 1.0e-5_f32,
        core::f32::consts::FRAC_PI_2,
        core::f32::consts::FRAC_PI_2 + 1.0e-5_f32,
        core::f32::consts::PI,
        10.0_f32 * core::f32::consts::PI + 0.125_f32,
    ];

    for x in points {
        let got = sinf(x);
        let want = x.sin();
        assert_abs_close_f32(got, want, 2.0e-6_f32, &format!("x={x:?}"));
    }
}

#[test]
fn sin_is_stable_across_quadrants_and_reduction_points() {
    let points = [
        -1.0e-12_f64,
        0.0_f64,
        1.0e-12_f64,
        core::f64::consts::FRAC_PI_6,
        core::f64::consts::FRAC_PI_2 - 1.0e-12_f64,
        core::f64::consts::FRAC_PI_2,
        core::f64::consts::FRAC_PI_2 + 1.0e-12_f64,
        core::f64::consts::PI,
        10.0_f64 * core::f64::consts::PI + 0.125_f64,
    ];

    for x in points {
        let got = sin(x);
        let want = x.sin();
        assert_abs_close_f64(got, want, 2.0e-15_f64, &format!("x={x:?}"));
    }
}

#[test]
fn cosf_is_stable_across_quadrants_and_reduction_points() {
    let points = [
        -1.0e-6_f32,
        0.0_f32,
        1.0e-6_f32,
        core::f32::consts::FRAC_PI_3,
        core::f32::consts::FRAC_PI_2 - 1.0e-5_f32,
        core::f32::consts::FRAC_PI_2,
        core::f32::consts::FRAC_PI_2 + 1.0e-5_f32,
        core::f32::consts::PI,
        10.0_f32 * core::f32::consts::PI + 0.125_f32,
    ];

    for x in points {
        let got = cosf(x);
        let want = x.cos();
        assert_abs_close_f32(got, want, 2.0e-6_f32, &format!("x={x:?}"));
    }
}

#[test]
fn cos_is_stable_across_quadrants_and_reduction_points() {
    let points = [
        -1.0e-12_f64,
        0.0_f64,
        1.0e-12_f64,
        core::f64::consts::FRAC_PI_3,
        core::f64::consts::FRAC_PI_2 - 1.0e-12_f64,
        core::f64::consts::FRAC_PI_2,
        core::f64::consts::FRAC_PI_2 + 1.0e-12_f64,
        core::f64::consts::PI,
        10.0_f64 * core::f64::consts::PI + 0.125_f64,
    ];

    for x in points {
        let got = cos(x);
        let want = x.cos();
        assert_abs_close_f64(got, want, 2.0e-15_f64, &format!("x={x:?}"));
    }
}

#[test]
fn asinf_matches_core_near_domain_edges() {
    let points = [
        1.0_f32 - 5.0e-5_f32,
        1.0_f32 - 5.0e-7_f32,
        -1.0_f32 + 5.0e-5_f32,
        -1.0_f32 + 5.0e-7_f32,
    ];

    for x in points {
        let got = asinf(x);
        let want = x.asin();
        assert!(
            ulp_diff_f32(got, want) <= 2,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn asin_matches_core_near_domain_edges() {
    let points = [
        1.0_f64 - 5.0e-5_f64,
        1.0_f64 - 5.0e-12_f64,
        -1.0_f64 + 5.0e-5_f64,
        -1.0_f64 + 5.0e-12_f64,
    ];

    for x in points {
        assert_eq!(asin(x), x.asin(), "x={x:?}");
    }
}

#[test]
fn acosf_matches_core_near_domain_edges() {
    let points = [
        1.0_f32 - 5.0e-5_f32,
        1.0_f32 - 5.0e-7_f32,
        -1.0_f32 + 5.0e-5_f32,
        -1.0_f32 + 5.0e-7_f32,
    ];

    for x in points {
        assert_eq!(acosf(x), x.acos(), "x={x:?}");
    }
}

#[test]
fn acos_matches_core_near_domain_edges() {
    let points = [
        1.0_f64 - 5.0e-5_f64,
        1.0_f64 - 5.0e-12_f64,
        -1.0_f64 + 5.0e-5_f64,
        -1.0_f64 + 5.0e-12_f64,
    ];

    for x in points {
        let got = acos(x);
        let want = x.acos();
        assert!(
            ulp_diff_f64(got, want) <= 64,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn sinhf_matches_core_near_zero() {
    let points = [-9.0e-7_f32, 9.0e-7_f32, -5.0e-5_f32, 5.0e-5_f32];

    for x in points {
        assert_eq!(sinhf(x), x.sinh(), "x={x:?}");
    }
}

#[test]
fn sinh_matches_core_near_zero() {
    let points = [
        -1.0e-12_f64,
        1.0e-12_f64,
        -5.0e-5_f64,
        5.0e-5_f64,
        -2.7e-3_f64,
        2.7e-3_f64,
    ];

    for x in points {
        let got = sinh(x);
        let want = x.sinh();
        assert!(
            ulp_diff_f64(got, want) <= 2,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn tanhf_matches_core_near_zero() {
    let points = [-9.0e-7_f32, 9.0e-7_f32, -5.0e-5_f32, 5.0e-5_f32];

    for x in points {
        assert_eq!(tanhf(x), x.tanh(), "x={x:?}");
    }
}

#[test]
fn tanh_matches_core_near_zero() {
    let points = [
        -1.0e-12_f64,
        1.0e-12_f64,
        -5.0e-5_f64,
        5.0e-5_f64,
        -3.95e-3_f64,
        3.95e-3_f64,
    ];

    for x in points {
        let got = tanh(x);
        let want = x.tanh();
        assert!(
            ulp_diff_f64(got, want) <= 2,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}

#[test]
fn sin_matches_core_near_zero_crossings() {
    let points = [
        core::f64::consts::PI - 1.0e-12_f64,
        core::f64::consts::PI,
        core::f64::consts::PI + 1.0e-12_f64,
    ];

    for x in points {
        let got = sin(x);
        let want = x.sin();
        assert_abs_close_f64(got, want, 2.0e-15_f64, &format!("x={x:?}"));
    }
}

#[test]
fn cos_matches_core_near_zero_crossings() {
    let points = [
        core::f64::consts::FRAC_PI_2 - 1.0e-12_f64,
        core::f64::consts::FRAC_PI_2,
        core::f64::consts::FRAC_PI_2 + 1.0e-12_f64,
    ];

    for x in points {
        let got = cos(x);
        let want = x.cos();
        assert_abs_close_f64(got, want, 2.0e-15_f64, &format!("x={x:?}"));
    }
}

#[test]
fn tan_matches_core_near_zero() {
    let points = [
        -1.7e-3_f64,
        -9.0e-4_f64,
        -2.5e-4_f64,
        2.5e-4_f64,
        9.0e-4_f64,
        1.7e-3_f64,
    ];

    for x in points {
        assert_eq!(tan(x), x.tan(), "x={x:?}");
    }
}

#[test]
fn tanf_matches_core_near_zero() {
    let points = [
        -1.7e-3_f32,
        -9.0e-4_f32,
        -2.5e-4_f32,
        2.5e-4_f32,
        9.0e-4_f32,
        1.7e-3_f32,
    ];

    for x in points {
        assert_eq!(tanf(x), x.tan(), "x={x:?}");
    }
}

#[test]
fn coshl_matches_core_for_small_f128_arguments() {
    let points = [
        3.333333333333333333333333333333333e-1_f128,
        5.0e-1_f128,
        1.0_f128,
        f128::from_bits(0x40000000000000000000000000000001u128),
        4.0_f128,
    ];

    for x in points {
        let got = coshl(x);
        let want = x.cosh();
        assert!(
            ulp_diff_f128(got, want) <= 2,
            "x={x:?} got={got:?} want={want:?}"
        );
    }
}
