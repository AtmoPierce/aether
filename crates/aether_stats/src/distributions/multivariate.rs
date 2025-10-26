#![cfg(feature = "std")]
use aether_core::math::{Vector, Matrix};
use num_traits::{Float, cast::cast };

/// Multivariate Normal distribution with mean vector and covariance matrix.
#[derive(Clone, Debug)]
pub struct MultivariateNormal<F: Float + , const N: usize> {
    pub mean: Vector<F, N>,
    pub cov: Matrix<F, N, N>,
}

impl<F: Float + , const N: usize> MultivariateNormal<F, N> {
    pub fn new(mean: Vector<F, N>, cov: Matrix<F, N, N>) -> Self { Self { mean, cov } }

    /// Evaluate the density at x using the formula
    /// (2π)^{-N/2} |Σ|^{-1/2} exp(-0.5 (x-μ)^T Σ^{-1} (x-μ))
    pub fn pdf(&self, x: &Vector<F, N>) -> F {
        let mut diff = x.clone();
        for i in 0..N { diff[i] = diff[i] - self.mean[i]; }
        let inv = self.cov.clone().inverse_gauss_jordan().expect("Error: Covariance matrix is singular");
        let quad = {
            let tmp = inv * diff.clone();
            diff.dot(&tmp)
        };
        let det = self.cov.determinant();
        let two = F::one() + F::one();
        let pi = cast(std::f64::consts::PI).unwrap();
        let denom = ( (two * pi).powf(cast(N as u32).unwrap()) * det ).sqrt();
        ( (-F::from(0.5).unwrap()) * quad ).exp() / denom
    }
}

pub type MultivariateNormalF64<const N: usize> = MultivariateNormal<f64, N>;

#[cfg(test)]
mod tests {
    use super::*;
    use aether_core::math::{Vector, Matrix};

    #[test]
    fn mvn_pdf_1d_matches_univariate() {
        let mean = Vector::new([0.0]);
        let cov = Matrix::identity();
        let mvn = MultivariateNormal::new(mean, cov);
        let x = Vector::new([0.0]);
        let p = mvn.pdf(&x);
        let un = crate::distributions::gaussian::Gaussian::new(0.0, 1.0);
        assert!((p - un.pdf(0.0)).abs() < 1e-12);
    }
}
