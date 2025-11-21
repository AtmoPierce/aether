#![cfg(feature = "std")]
use aether_core::math::Vector;
use aether_opt::gradient_descent::GradientDescentGeneric;
use aether_core::real::{Real};
use std::vec::Vec;

/// Simple logistic regression for small fixed-size feature vectors.
/// Generic over N (number of features) and Real type F.
#[derive(Clone, Copy, Debug)]
pub struct LogisticRegression<F: Real + Copy, const N: usize> {
    pub weights: Vector<F, N>,
    pub bias: F,
    pub lr: F,
}

impl<F: Real + Copy, const N: usize> LogisticRegression<F, N> {
    pub fn new() -> Self {
        Self {
            weights: Vector::new([F::ZERO; N]),
            bias: F::ZERO,
            // neutral default lr
            lr: F::ONE,
        }
    }

    #[inline]
    fn sigmoid(x: F) -> F {
        F::ONE / (F::ONE + (-x).exp())
    }

    /// Predict probability for a single example
    pub fn predict_probability(&self, x: &Vector<F, N>) -> F {
        let mut s = self.bias;
        for i in 0..N {
            s = s + self.weights[i] * x[i];
        }
        Self::sigmoid(s)
    }

    /// Predict 0/1 label
    pub fn predict(&self, x: &Vector<F, N>) -> u8 {
        if self.predict_probability(x) > (F::ONE / (F::ONE + F::ONE)) {
            1
        } else {
            0
        }
    }

    /// One epoch of batch gradient descent over m examples given as slices
    /// x_data: &[Vector<F, N>]  y: &[u8]
    pub fn train_epoch(&mut self, x_data: &[Vector<F, N>], y: &[u8]) {
        let m_usize = x_data.len();
        if m_usize == 0 {
            return;
        }
        let m: F = F::from_usize(m_usize);

        let mut grad_w = Vector::new([F::ZERO; N]);
        let mut grad_b = F::ZERO;

        for (x, &label) in x_data.iter().zip(y.iter()) {
            let p = self.predict_probability(x);
            let y_f: F = if label == 0 { F::ZERO } else { F::ONE };
            let err = p - y_f;
            for i in 0..N {
                grad_w[i] = grad_w[i] + err * x[i];
            }
            grad_b = grad_b + err;
        }

        for i in 0..N {
            self.weights[i] = self.weights[i] - self.lr * (grad_w[i] / m);
        }
        self.bias = self.bias - self.lr * (grad_b / m);
    }

    /// Fit using existing GradientDescent optimizer over augmented params [w..., bias]
    pub fn fit_with_optimizer<const M: usize>(
        &mut self,
        x_data: &[Vector<F, N>],
        y: &[u8],
        mut opt: GradientDescentGeneric<F, M>,
    ) {
        assert_eq!(M, N + 1, "fit_with_optimizer: expected M == N+1");
        let m_usize = x_data.len();
        if m_usize == 0 {
            return;
        }
        let m: F = F::from_usize(m_usize);
        let mut p_arr: [F; M] = [F::ZERO; M];
        for i in 0..N {
            p_arr[i] = self.weights[i];
        }
        p_arr[N] = self.bias;
        let x0 = Vector::new(p_arr);

        let tiny = F::ONE / (F::ONE + (F::ONE + F::ONE));

        let f = |pv: &Vector<F, M>| -> F {
            let mut loss = F::ZERO;
            for (x, &yi) in x_data.iter().zip(y.iter()) {
                let mut s = pv[N];
                for i in 0..N {
                    s = s + pv[i] * x[i];
                }
                let p = F::ONE / (F::ONE + (-s).exp());
                let yv: F = if yi == 0 { F::ZERO } else { F::ONE };
                // binary cross-entropy: -(y ln p + (1-y) ln (1-p))
                loss = loss - (yv * (p + tiny).ln() + (F::ONE - yv) * (F::ONE - p + tiny).ln());
            }
            loss / m
        };

        let g = |pv: &Vector<F, M>| -> Vector<F, M> {
            let mut grads: [F; M] = [F::ZERO; M];
            for (x, &yi) in x_data.iter().zip(y.iter()) {
                let mut s = pv[N];
                for i in 0..N {
                    s = s + pv[i] * x[i];
                }
                let p = F::ONE / (F::ONE + (-s).exp());
                let e = p - if yi == 0 { F::ZERO } else { F::ONE };
                for i in 0..N {
                    grads[i] = grads[i] + e * x[i];
                }
                grads[N] = grads[N] + e;
            }
            for i in 0..M {
                grads[i] = grads[i] / m;
            }
            Vector::new(grads)
        };

        let (p_opt, _fval, _iters, _conv) = opt.minimize(x0, f, g);
        for i in 0..N {
            self.weights[i] = p_opt[i];
        }
        self.bias = p_opt[N];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use aether_core::math::Vector;

    #[test]
    fn train_linearly_separable() {
        // dataset: AND-like gate (two features)
        let x: [Vector<f64, 2>; 4] = [
            Vector::new([0.0, 0.0]),
            Vector::new([0.0, 1.0]),
            Vector::new([1.0, 0.0]),
            Vector::new([1.0, 1.0]),
        ];
        let y: [u8; 4] = [0, 0, 0, 1];

        let mut model = LogisticRegression::<f64, 2> {
            weights: Vector::new([0.0, 0.0]),
            bias: 0.0,
            lr: 0.5,
        };

        for _ in 0..2000 {
            model.train_epoch(&x, &y);
        }

        // Check predictions
        let preds: Vec<u8> = x.iter().map(|xi| model.predict(xi)).collect();
        assert_eq!(preds, vec![0u8, 0, 0, 1]);
    }
}
