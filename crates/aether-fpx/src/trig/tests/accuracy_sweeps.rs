use crate::trig::{asinl, cosh, expm1, log, log1p, sinh, sinhl, tanh, tanl};

use super::support::{ulp_diff_f128, ulp_diff_f64};

fn assert_f64_sweep(name: &str, fpx: fn(f64) -> f64, reference: fn(f64) -> f64, lo: f64, hi: f64) {
    const SAMPLES: usize = 2001;
    for index in 0..SAMPLES {
        let x = lo + (hi - lo) * index as f64 / (SAMPLES - 1) as f64;
        let got = fpx(x);
        let want = reference(x);
        let ulp = ulp_diff_f64(got, want);
        assert!(
            ulp <= 5,
            "{name} exceeded 5 ULP at x={x:?}: got={got:?}, want={want:?}, ulp={ulp}"
        );
    }
}

fn assert_f128_sweep(
    name: &str,
    fpx: fn(f128) -> f128,
    reference: fn(f128) -> f128,
    lo: f128,
    hi: f128,
) {
    const SAMPLES: usize = 2001;
    for index in 0..SAMPLES {
        let x = lo + (hi - lo) * (index as f128) / ((SAMPLES - 1) as f128);
        let got = fpx(x);
        let want = reference(x);
        let ulp = ulp_diff_f128(got, want);
        assert!(
            ulp <= 5,
            "{name} exceeded 5 ULP at x={x:?}: got={got:?}, want={want:?}, ulp={ulp}"
        );
    }
}

#[test]
fn publication_intervals_stay_within_five_ulp_for_f64() {
    assert_f64_sweep("expm1", expm1, f64::exp_m1, -4.0, 4.0);
    assert_f64_sweep("sinh", sinh, f64::sinh, -4.0, 4.0);
    assert_f64_sweep("cosh", cosh, f64::cosh, -4.0, 4.0);
    assert_f64_sweep("tanh", tanh, f64::tanh, -4.0, 4.0);
    assert_f64_sweep("log", log, f64::ln, 1.0e-6, 10.0);
    assert_f64_sweep("log1p", log1p, f64::ln_1p, -0.999_999, 10.0);
}

#[test]
fn publication_intervals_stay_within_five_ulp_for_f128_regressions() {
    assert_f128_sweep("asin", asinl, f128::asin, -1.0, 1.0);
    assert_f128_sweep("sinh", sinhl, f128::sinh, -4.0, 4.0);
    assert_f128_sweep("tan", tanl, f128::tan, -1.2, 1.2);
}

#[test]
fn asin_f128_range_reduction_regression() {
    let x = f128::from_bits(0xbffe6353f7ced916872b020c49ba5e35);
    let got = asinl(x);
    let reference = x.asin();
    assert!(
        ulp_diff_f128(got, reference) <= 5,
        "asin exceeded 5 ULP at the former range-reduction failure: got={got:?}, reference={reference:?}"
    );
}
