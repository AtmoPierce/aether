#![cfg(feature = "f128")]

use crate::remez::RemezApprox;
use aether_core::real::Real;
use core::f128;

pub struct RemezQuantResult<F: Real, const N: usize, const K: usize> {
    pub coeffs: [F; N],
    pub max_err_vs_poly: F,
    pub max_err_vs_target: F,
}

#[inline]
fn eval_poly_generic<F: Real, const N: usize>(coeffs: &[F; N], x: F) -> F {
    // p(x) = c0 + c1 x + c2 x^2 + ...
    let mut acc = F::ZERO;
    let mut power = F::ONE;
    for i in 0..N {
        acc = acc + coeffs[i] * power;
        power = power * x;
    }
    acc
}

fn max_error_vs_poly<F: Real, const N: usize, const K: usize>(
    rz128: &RemezApprox<f128, N, K>,
    coeffs: &[F; N],
    a: f128,
    b: f128,
    samples: usize,
) -> F {
    let n = samples.max(1);
    let mut max_err = F::ZERO;

    for j in 0..n {
        let t64 = j as f64 / (n as f64 - 1.0).max(1.0);
        let t128 = t64 as f128;
        let x128 = a + (b - a) * t128;
        let xF = F::from_f128(x128);

        let p_ref_128 = rz128.eval_poly(x128);
        let p_ref = F::from_f128(p_ref_128);

        let p_q = eval_poly_generic::<F, N>(coeffs, xF);
        let e = (p_q - p_ref).abs();
        if e > max_err {
            max_err = e;
        }
    }

    max_err
}

fn max_error_vs_target<F, G, const N: usize, const K: usize>(
    coeffs: &[F; N],
    f: G,
    a: f128,
    b: f128,
    samples: usize,
) -> F
where
    F: Real,
    G: Fn(f128) -> f128,
{
    let n = samples.max(1);
    let mut max_err = F::ZERO;

    for j in 0..n {
        let t64 = j as f64 / (n as f64 - 1.0).max(1.0);
        let t128 = t64 as f128;
        let x128 = a + (b - a) * t128;
        let xF = F::from_f128(x128);

        let y_true_128 = f(x128);
        let y_true = F::from_f128(y_true_128);

        let y_q = eval_poly_generic::<F, N>(coeffs, xF);
        let e = (y_q - y_true).abs();
        if e > max_err {
            max_err = e;
        }
    }

    max_err
}

pub fn optimize_quantized_coeffs<F, G, const N: usize, const K: usize>(
    rz128: &RemezApprox<f128, N, K>,
    f: G,
    a: f128,
    b: f128,
    samples: usize,
    iters: usize,
    step: F,
) -> RemezQuantResult<F, N, K>
where
    F: Real,
    G: Fn(f128) -> f128,
{
    let n = samples.max(1);

    let mut coeffs = [F::ZERO; N];
    for i in 0..N {
        coeffs[i] = F::from_f128(rz128.coefficients[i]);
    }

    let inv_n = F::ONE / F::from_usize(n);
    let two = F::from_f32(2.0);

    for _iter in 0..iters {
        let mut grad = [F::ZERO; N];
        for j in 0..n {
            let t64 = j as f64 / (n as f64 - 1.0).max(1.0);
            let t128 = t64 as f128;
            let x128 = a + (b - a) * t128;
            let xF = F::from_f128(x128);

            let y_true_128 = f(x128);
            let y_true = F::from_f128(y_true_128);

            let y_q = eval_poly_generic::<F, N>(&coeffs, xF);
            let err = y_q - y_true;

            let mut basis = F::ONE;
            for i in 0..N {
                grad[i] = grad[i] + two * err * basis;
                basis = basis * xF;
            }
        }

        for i in 0..N {
            coeffs[i] = coeffs[i] - step * inv_n * grad[i];
        }
    }

    let max_err_vs_poly = max_error_vs_poly::<F, N, K>(rz128, &coeffs, a, b, n);
    let max_err_vs_target = max_error_vs_target::<F, G, N, K>(&coeffs, f, a, b, n);

    RemezQuantResult {
        coeffs,
        max_err_vs_poly,
        max_err_vs_target,
    }
}
