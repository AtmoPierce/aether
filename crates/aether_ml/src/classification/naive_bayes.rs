#![cfg(feature = "std")]
use aether_core::math::Vector;

/// Gaussian Naive Bayes for fixed-size feature vectors.
pub struct GaussianNB<const N: usize> {
    pub class_count: Vec<usize>,
    pub class_prior: Vec<f64>,
    pub class_mean: Vec<Vector<f64, N>>,
    pub class_var: Vec<Vector<f64, N>>, // variance (not covariance)
}

impl<const N: usize> GaussianNB<N> {
    pub fn fit(x_data: &[Vector<f64, N>], y: &[usize], num_classes: usize) -> Self {
        let mut class_count = vec![0usize; num_classes];
        let mut class_mean: Vec<Vector<f64, N>> = vec![Vector::new([0.0; N]); num_classes];
        let mut class_var: Vec<Vector<f64, N>> = vec![Vector::new([0.0; N]); num_classes];

        for (x, &label) in x_data.iter().zip(y.iter()) {
            class_count[label] += 1;
            for i in 0..N {
                class_mean[label][i] += x[i];
            }
        }
        for c in 0..num_classes {
            if class_count[c] > 0 {
                for i in 0..N {
                    class_mean[c][i] /= class_count[c] as f64;
                }
            }
        }
        // compute variance (MLE: divide by n)
        for (x, &label) in x_data.iter().zip(y.iter()) {
            for i in 0..N {
                let d = x[i] - class_mean[label][i];
                class_var[label][i] += d * d;
            }
        }
        for c in 0..num_classes {
            if class_count[c] > 0 {
                for i in 0..N {
                    class_var[c][i] /= class_count[c] as f64;
                }
            }
        }

        let total: usize = class_count.iter().sum();
        let class_prior: Vec<f64> = class_count
            .iter()
            .map(|&n| (n as f64) / (total as f64))
            .collect();

        GaussianNB {
            class_count,
            class_prior,
            class_mean,
            class_var,
        }
    }

    pub fn predict(&self, x: &Vector<f64, N>) -> usize {
        // compute log-likelihood for each class (sum log Gaussian probabilities)
        let mut best = 0usize;
        let mut best_score = f64::NEG_INFINITY;
        for (c, prior) in self.class_prior.iter().enumerate() {
            if self.class_count[c] == 0 {
                continue;
            }
            let mut score = prior.ln();
            for i in 0..N {
                let var = self.class_var[c][i].max(1e-12);
                let diff = x[i] - self.class_mean[c][i];
                // log Gaussian (univariate)
                score += -0.5 * ((2.0 * std::f64::consts::PI * var).ln() + (diff * diff) / var);
            }
            if score > best_score {
                best_score = score;
                best = c;
            }
        }
        best
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use aether_core::math::Vector;

    #[test]
    fn gnb_simple() {
        // 1D data, class 0 centered at 0, class1 centered at 5
        let xs: [Vector<f64, 1>; 6] = [
            Vector::new([0.0]),
            Vector::new([0.2]),
            Vector::new([-0.1]),
            Vector::new([5.0]),
            Vector::new([5.2]),
            Vector::new([4.8]),
        ];
        let ys: [usize; 6] = [0, 0, 0, 1, 1, 1];
        let model = GaussianNB::<1>::fit(&xs, &ys, 2);
        assert_eq!(model.predict(&xs[0]), 0);
        assert_eq!(model.predict(&xs[3]), 1);
    }
}
