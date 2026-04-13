use aether_core::math::Vector;

/// Barycentric Lagrange polynomial interpolation over scalar samples.
#[derive(Debug, Clone, Copy)]
pub struct PolynomialInterpolator<const N: usize> {
    xs: Vector<f64, N>,
    ys: Vector<f64, N>,
    weights: Vector<f64, N>,
}

impl<const N: usize> PolynomialInterpolator<N> {
    /// Builds a polynomial interpolator that passes through all sample points.
    pub fn new(xs: Vector<f64, N>, ys: Vector<f64, N>) -> Result<Self, &'static str> {
        if N == 0 {
            return Err("at least one sample is required");
        }

        let mut weights = Vector::new([1.0_f64; N]);
        for j in 0..N {
            for k in 0..N {
                if j == k {
                    continue;
                }
                let diff = xs[j] - xs[k];
                if diff.abs() <= f64::EPSILON {
                    return Err("x samples must be unique");
                }
                weights[j] /= diff;
            }
        }

        Ok(Self { xs, ys, weights })
    }

    /// Evaluates the interpolating polynomial at `x`.
    pub fn evaluate(&self, x: f64) -> f64 {
        for i in 0..N {
            if (x - self.xs[i]).abs() <= f64::EPSILON {
                return self.ys[i];
            }
        }

        let mut numerator = 0.0_f64;
        let mut denominator = 0.0_f64;
        for i in 0..N {
            let term = self.weights[i] / (x - self.xs[i]);
            numerator += term * self.ys[i];
            denominator += term;
        }
        numerator / denominator
    }

    pub const fn degree(&self) -> usize {
        N.saturating_sub(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polynomial_interpolates_quadratic() {
        let xs = Vector::new([-1.0, 0.0, 2.0]);
        let ys = Vector::new([0.0, 1.0, 9.0]);
        let interp = PolynomialInterpolator::new(xs, ys).unwrap();
        assert!((interp.evaluate(1.5) - 6.25).abs() < 1.0e-12);
    }

    #[test]
    fn polynomial_returns_exact_knot_value() {
        let interp = PolynomialInterpolator::new(Vector::new([0.0, 2.0]), Vector::new([1.0, 5.0])).unwrap();
        assert_eq!(interp.evaluate(2.0), 5.0);
    }

    #[test]
    fn polynomial_rejects_duplicate_x_values() {
        assert!(PolynomialInterpolator::new(Vector::new([0.0, 0.0]), Vector::new([1.0, 2.0])).is_err());
    }
}
