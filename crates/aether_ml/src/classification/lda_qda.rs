#![cfg(feature = "std")]
use aether_core::math::{Matrix, Vector};

/// Linear Discriminant Analysis (Gaussian with shared covariance)
pub struct Lda<const N: usize> {
    pub means: Vec<Vector<f64, N>>,
    pub priors: Vec<f64>,
    pub pooled_cov: Matrix<f64, N, N>,
}

/// Quadratic Discriminant Analysis (class-specific covariances)
pub struct Qda<const N: usize> {
    pub means: Vec<Vector<f64, N>>,
    pub priors: Vec<f64>,
    pub covs: Vec<Matrix<f64, N, N>>,
}

impl<const N: usize> Lda<N> {
    pub fn fit(x_data: &[Vector<f64, N>], y: &[usize], num_classes: usize) -> Self {
        let mut means: Vec<Vector<f64, N>> = Vec::with_capacity(num_classes);
        let mut priors: Vec<f64> = Vec::with_capacity(num_classes);
        let mut pooled = Matrix::<f64, N, N>::zeros();

        for c in 0..num_classes {
            // collect class samples
            let mut sum = Vector::new([0.0; N]);
            let mut count = 0usize;
            let mut cov = Matrix::<f64, N, N>::zeros();
            for (x, &label) in x_data.iter().zip(y.iter()) {
                if label == c {
                    count += 1;
                    for i in 0..N {
                        sum[i] += x[i];
                    }
                }
            }
            let prior = (count as f64) / (x_data.len() as f64);
            priors.push(prior);
            if count == 0 {
                means.push(Vector::new([0.0; N]));
                continue;
            }
            for i in 0..N {
                sum[i] /= (count as f64);
            }
            means.push(sum);

            // compute class covariance
            for (x, &label) in x_data.iter().zip(y.iter()) {
                if label == c {
                    let mut diff = Vector::new([0.0; N]);
                    for i in 0..N {
                        diff[i] = x[i] - means[c][i];
                    }
                    // outer product diff * diff^T
                    for i in 0..N {
                        for j in 0..N {
                            cov[(i, j)] += diff[i] * diff[j];
                        }
                    }
                }
            }
            // normalize
            if count > 1 {
                cov = cov / (count as f64 - 1.0);
            }
            pooled = pooled + cov * prior;
        }

        Lda {
            means,
            priors,
            pooled_cov: pooled,
        }
    }

    /// Classify a single sample, returns class index with max discriminant
    pub fn predict(&self, x: &Vector<f64, N>) -> usize {
        let mut best = 0usize;
        let mut best_score = f64::NEG_INFINITY;
        for (i, mu) in self.means.iter().enumerate() {
            // discriminant: x^T Sigma^{-1} mu - 0.5 mu^T Sigma^{-1} mu + log prior
            if let Some(sinv) = self.pooled_cov.inverse_gauss_jordan() {
                let mut xm = Vector::new([0.0; N]);
                for k in 0..N {
                    xm[k] = x[k];
                }
                // compute x^T S^{-1} mu
                let tmp = (sinv * mu.clone()); // Matrix * Vector -> Vector
                let mut term1 = 0.0;
                for k in 0..N {
                    term1 += x[k] * tmp[k];
                }

                // mu^T S^{-1} mu
                let tmp2 = (sinv * mu.clone());
                let mut term2 = 0.0;
                for k in 0..N {
                    term2 += mu[k] * tmp2[k];
                }

                let score = term1 - 0.5 * term2 + self.priors[i].ln();
                if score > best_score {
                    best_score = score;
                    best = i;
                }
            }
        }
        best
    }
}

impl<const N: usize> Qda<N> {
    pub fn fit(x_data: &[Vector<f64, N>], y: &[usize], num_classes: usize) -> Self {
        let mut means: Vec<Vector<f64, N>> = Vec::with_capacity(num_classes);
        let mut priors: Vec<f64> = Vec::with_capacity(num_classes);
        let mut covs: Vec<Matrix<f64, N, N>> = Vec::with_capacity(num_classes);

        for c in 0..num_classes {
            let mut sum = Vector::new([0.0; N]);
            let mut count = 0usize;
            for (x, &label) in x_data.iter().zip(y.iter()) {
                if label == c {
                    count += 1;
                    for i in 0..N {
                        sum[i] += x[i];
                    }
                }
            }
            if count == 0 {
                means.push(Vector::new([0.0; N]));
                priors.push(0.0);
                covs.push(Matrix::zeros());
                continue;
            }
            for i in 0..N {
                sum[i] /= (count as f64);
            }
            means.push(sum);
            priors.push((count as f64) / (x_data.len() as f64));

            let mut cov = Matrix::<f64, N, N>::zeros();
            for (x, &label) in x_data.iter().zip(y.iter()) {
                if label == c {
                    let mut diff = Vector::new([0.0; N]);
                    for i in 0..N {
                        diff[i] = x[i] - means[c][i];
                    }
                    for i in 0..N {
                        for j in 0..N {
                            cov[(i, j)] += diff[i] * diff[j];
                        }
                    }
                }
            }
            if count > 1 {
                cov = cov / (count as f64 - 1.0);
            }
            covs.push(cov);
        }
        Qda {
            means,
            priors,
            covs,
        }
    }

    pub fn predict(&self, x: &Vector<f64, N>) -> usize {
        let mut best = 0usize;
        let mut best_score = f64::NEG_INFINITY;
        for (i, mu) in self.means.iter().enumerate() {
            let cov = &self.covs[i];
            if let Some(cinv) = cov.inverse_gauss_jordan() {
                // score = -0.5 ln|Sigma| - 0.5 (x-mu)^T Sigma^{-1} (x-mu) + ln prior
                // determinant for supported sizes (generic fallback)
                let det = cov.determinant();
                let mut diff = Vector::new([0.0; N]);
                for k in 0..N {
                    diff[k] = x[k] - mu[k];
                }
                // tmp = Sigma^{-1} * diff
                let tmp = cinv * diff.clone();
                let mut quad = 0.0;
                for k in 0..N {
                    quad += diff[k] * tmp[k];
                }
                let score = -0.5 * det.ln() - 0.5 * quad + self.priors[i].ln();
                if score > best_score {
                    best_score = score;
                    best = i;
                }
            }
        }
        best
    }
}
