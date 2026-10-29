// aether-fpx-gen/src/design/remez/create/accuracy.rs
//! Accuracy analysis: ULP and error metrics per segment

use core::f128;

use aether_core::real::real::{Real, RealCast};
use aether_fpx_bench::ulp::ulp_dist_f64;

use super::coefficients::TrigKernel;
use crate::remez::RemezApprox;

/// Accuracy metrics for one segment
#[derive(Clone, Debug)]
pub struct SegmentAccuracy {
    pub segment_idx: usize,
    pub n_coeffs: usize, // polynomial degree N
    pub max_ulp: f64,
    pub max_abs_error: f64,
    pub mean_abs_error: f64,
    pub max_rel_error: f64,
}

impl SegmentAccuracy {
    pub fn report(&self) -> String {
        format!(
            "  Segment {}: N={:2}, max_ULP={:7.2}, max_err={:.2e}, mean_err={:.2e}, max_rel_err={:.2e}",
            self.segment_idx,
            self.n_coeffs,
            self.max_ulp,
            self.max_abs_error,
            self.mean_abs_error,
            self.max_rel_error,
        )
    }
}

/// Analyze accuracy of a segment polynomial
///
/// IMPORTANT:
/// This MUST mirror runtime kernel semantics exactly.
/// Any mismatch here will invalidate ULP guarantees.
pub fn analyze_segment_accuracy<const N: usize, const K: usize>(
    rz: &RemezApprox<f128, N, K>,
    kernel: TrigKernel,
    segment_idx: usize,
    interval_x: (f128, f128),
    f_runtime: fn(f128) -> f128,
    samples: usize,
) -> SegmentAccuracy {
    let (ax, bx) = interval_x;

    // ULP is meaningless near zero due to dense spacing
    const NEAR_ZERO_THRESHOLD: f64 = 1.0e-10;

    let mut ulps = Vec::with_capacity(samples);
    let mut abs_errors = Vec::with_capacity(samples);
    let mut rel_errors = Vec::with_capacity(samples);

    for i in 0..=samples {
        let t = i as f64 / samples as f64;
        let x_f64 = ax.to_f64() + (bx.to_f64() - ax.to_f64()) * t;
        let x: f128 = x_f64.into();

        // ---- Approximation (mirrors kernel evaluation)
        let approx_f128 = match kernel {
            TrigKernel::Sin => {
                if x == f128::ZERO {
                    f128::ZERO
                } else {
                    let z = x * x;
                    x * rz.eval_poly(z)
                }
            }
            TrigKernel::Cos => {
                let z = x * x;
                rz.eval_poly(z)
            }
            TrigKernel::Tan => {
                let z = x * x;
                x + x * z * rz.eval_poly(z)
            }
        };

        let approx = approx_f128.to_f64();
        let true_val = f_runtime(x).to_f64();
        let err = (approx - true_val).abs();

        let ulp = if true_val.abs() < NEAR_ZERO_THRESHOLD {
            0.0
        } else {
            ulp_dist_f64(approx, true_val) as f64
        };

        abs_errors.push(err);
        ulps.push(ulp);

        if true_val != 0.0 {
            rel_errors.push(err / true_val.abs());
        }
    }

    let max_ulp = ulps.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    let max_abs_error = abs_errors.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    let mean_abs_error = if abs_errors.is_empty() {
        0.0
    } else {
        abs_errors.iter().sum::<f64>() / abs_errors.len() as f64
    };

    let max_rel_error = if rel_errors.is_empty() {
        0.0
    } else {
        rel_errors.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    };

    SegmentAccuracy {
        segment_idx,
        n_coeffs: N,
        max_ulp,
        max_abs_error,
        mean_abs_error,
        max_rel_error,
    }
}

/// Check if a segment meets ULP acceptance criteria
#[inline]
pub fn meets_ulp_acceptance(max_ulp: f64, target_max_ulp: f64) -> bool {
    max_ulp <= target_max_ulp
}

/// Function-wide accuracy report
pub struct FunctionAccuracyReport {
    pub func_name: String,
    pub segments: Vec<SegmentAccuracy>,
}

impl FunctionAccuracyReport {
    pub fn new(func_name: String) -> Self {
        Self {
            func_name,
            segments: Vec::new(),
        }
    }

    pub fn add_segment(&mut self, acc: SegmentAccuracy) {
        self.segments.push(acc);
    }

    pub fn report(&self) -> String {
        let mut out = format!("=== {} Accuracy Report ===\n", self.func_name);

        if self.segments.is_empty() {
            out.push_str("  (no segments)\n");
            return out;
        }

        let max_ulp_overall = self
            .segments
            .iter()
            .map(|s| s.max_ulp)
            .fold(f64::NEG_INFINITY, f64::max);

        out.push_str(&format!("Overall: max_ULP={:.2}\n\n", max_ulp_overall));

        for seg in &self.segments {
            out.push_str(&seg.report());
            out.push('\n');
        }

        out
    }
}
