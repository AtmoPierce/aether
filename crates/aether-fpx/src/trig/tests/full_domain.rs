use crate::trig::{acos, asin, atan, cos, cosh, exp, expm1, log, log1p, sin, sinh, tan, tanh};

fn assert_same_specials(fpx: fn(f64) -> f64, reference: fn(f64) -> f64, points: &[f64]) {
    for &x in points {
        let got = fpx(x);
        let want = reference(x);
        if want.is_nan() {
            assert!(got.is_nan(), "x={x:?} got={got:?} want=NaN");
        } else {
            assert_eq!(got.to_bits(), want.to_bits(), "x={x:?}");
        }
    }
}

#[test]
fn generated_functions_match_full_domain_endpoint_semantics() {
    let all_reals = [f64::NEG_INFINITY, -0.0, 0.0, f64::INFINITY, f64::NAN];
    assert_same_specials(sin, f64::sin, &all_reals);
    assert_same_specials(cos, f64::cos, &all_reals);
    assert_same_specials(tan, f64::tan, &all_reals);
    assert_same_specials(sinh, f64::sinh, &all_reals);
    assert_same_specials(cosh, f64::cosh, &all_reals);
    assert_same_specials(tanh, f64::tanh, &all_reals);
    assert_same_specials(atan, f64::atan, &all_reals);
    assert_same_specials(exp, f64::exp, &all_reals);
    assert_same_specials(expm1, f64::exp_m1, &all_reals);

    let unit_domain = [
        f64::NEG_INFINITY,
        f64::from_bits((-1.0_f64).to_bits() + 1),
        -1.0,
        -0.0,
        0.0,
        1.0,
        f64::from_bits(1.0_f64.to_bits() + 1),
        f64::INFINITY,
        f64::NAN,
    ];
    assert_same_specials(asin, f64::asin, &unit_domain);
    assert_same_specials(acos, f64::acos, &unit_domain);

    let log_domain = [
        f64::NEG_INFINITY,
        -1.0,
        -0.0,
        0.0,
        f64::from_bits(1),
        1.0,
        f64::INFINITY,
        f64::NAN,
    ];
    assert_same_specials(log, f64::ln, &log_domain);

    let log1p_domain = [
        f64::NEG_INFINITY,
        f64::from_bits((-1.0_f64).to_bits() + 1),
        -1.0,
        -0.0,
        0.0,
        f64::INFINITY,
        f64::NAN,
    ];
    assert_same_specials(log1p, f64::ln_1p, &log1p_domain);
}
