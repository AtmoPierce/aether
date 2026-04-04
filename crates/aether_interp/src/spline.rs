use aether_core::math::Vector;

/// A cubic spline over scalar samples.
#[derive(Debug, Clone, Copy)]
pub struct CubicSpline<const N: usize> {
    xs: Vector<f64, N>,
    ys: Vector<f64, N>,
    second_derivatives: Vector<f64, N>,
}

impl<const N: usize> CubicSpline<N> {
    /// Builds a cubic spline through the given samples.
    pub fn new(xs: Vector<f64, N>, ys: Vector<f64, N>) -> Result<Self, &'static str> {
        if N < 2 {
            return Err("at least two samples are required");
        }
        for i in 0..N - 1 {
            if xs[i + 1] <= xs[i] {
                return Err("x samples must be strictly increasing");
            }
        }

        let mut second_derivatives = Vector::new([0.0_f64; N]);
        let mut u = Vector::new([0.0_f64; N]);

        for i in 1..N - 1 {
            let sig = (xs[i] - xs[i - 1]) / (xs[i + 1] - xs[i - 1]);
            let p = sig * second_derivatives[i - 1] + 2.0;
            second_derivatives[i] = (sig - 1.0) / p;

            let ddydx = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])
                - (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]);
            u[i] = (6.0 * ddydx / (xs[i + 1] - xs[i - 1]) - sig * u[i - 1]) / p;
        }

        for k in (0..N - 1).rev() {
            second_derivatives[k] = second_derivatives[k] * second_derivatives[k + 1] + u[k];
        }

        Ok(Self {
            xs,
            ys,
            second_derivatives,
        })
    }

    /// Evaluates the spline at `x`. Values outside the sample range are clamped
    /// to the nearest segment for stable endpoint extrapolation.
    pub fn evaluate(&self, x: f64) -> f64 {
        let idx = self.segment_index(x);
        let x0 = self.xs[idx];
        let x1 = self.xs[idx + 1];
        let y0 = self.ys[idx];
        let y1 = self.ys[idx + 1];
        let y2_0 = self.second_derivatives[idx];
        let y2_1 = self.second_derivatives[idx + 1];
        let h = x1 - x0;

        let a = (x1 - x) / h;
        let b = (x - x0) / h;
        a * y0
            + b * y1
            + ((a * a * a - a) * y2_0 + (b * b * b - b) * y2_1) * (h * h / 6.0)
    }

    pub const fn sample_count(&self) -> usize {
        N
    }

    fn segment_index(&self, x: f64) -> usize {
        if x <= self.xs[0] {
            return 0;
        }
        if x >= self.xs[N - 1] {
            return N - 2;
        }

        let mut lo = 0usize;
        let mut hi = N - 1;
        while lo + 1 < hi {
            let mid = (lo + hi) / 2;
            if self.xs[mid] <= x {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        lo
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_spline_reproduces_linear_samples() {
        let spline = CubicSpline::new(Vector::new([0.0, 1.0, 2.0]), Vector::new([1.0, 3.0, 5.0])).unwrap();
        assert!((spline.evaluate(0.5) - 2.0).abs() < 1.0e-12);
        assert!((spline.evaluate(1.5) - 4.0).abs() < 1.0e-12);
    }

    #[test]
    fn cubic_spline_hits_knots() {
        let spline = CubicSpline::new(Vector::new([0.0, 1.0, 3.0]), Vector::new([0.0, 2.0, 1.0])).unwrap();
        assert!((spline.evaluate(0.0) - 0.0).abs() < 1.0e-12);
        assert!((spline.evaluate(1.0) - 2.0).abs() < 1.0e-12);
        assert!((spline.evaluate(3.0) - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn cubic_spline_rejects_unsorted_samples() {
        assert!(CubicSpline::new(Vector::new([0.0, 1.0, 1.0]), Vector::new([0.0, 1.0, 2.0])).is_err());
    }
}
