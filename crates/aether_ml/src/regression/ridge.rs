#![cfg(feature = "std")]
use aether_core::math::Vector;
use aether_opt::gradient_descent::GradientDescentGeneric;
use aether_core::real::Real;

/// Ridge regression (L2) using batch gradient descent
#[derive(Clone, Copy, Debug)]
pub struct RidgeRegression<F: Real + Copy, const N: usize> {
    pub weights: Vector<F, N>,
    pub bias: F,
    pub lr: F,
    pub alpha: F,
}

impl<F: Real + Copy, const N: usize> RidgeRegression<F, N> {
    pub fn new() -> Self {
        // Neutral defaults: lr = 1, alpha = 1 (caller should set typical small lr/alpha)
        Self {
            weights: Vector::new([F::ZERO; N]),
            bias: F::ZERO,
            lr: F::ONE,
            alpha: F::ONE,
        }
    }

    pub fn predict(&self, x: &Vector<F, N>) -> F {
        let mut s = self.bias;
        for i in 0..N {
            s = s + self.weights[i] * x[i];
        }
        s
    }

    pub fn train_epoch(&mut self, x_data: &[Vector<F, N>], y: &[F]) {
        let m_usize = x_data.len();
        if m_usize == 0 {
            return;
        }
        let m: F = F::from_usize(m_usize);

        let mut grad_w = Vector::new([F::ZERO; N]);
        let mut grad_b = F::ZERO;

        for (x, &yi) in x_data.iter().zip(y.iter()) {
            let pred = self.predict(x);
            let err = pred - yi;
            for i in 0..N {
                grad_w[i] = grad_w[i] + err * x[i];
            }
            grad_b = grad_b + err;
        }

        let two = F::ONE + F::ONE;
        for i in 0..N {
            // L2 term gradient: 2*alpha*w
            self.weights[i] = self.weights[i]
                - self.lr * ((two * grad_w[i] / m) + two * self.alpha * self.weights[i]);
        }
        self.bias = self.bias - self.lr * (two * grad_b / m);
    }

    /// Fit using existing GradientDescent optimizer; packs weights and bias into size N+1 vector.
    pub fn fit_with_optimizer<const M: usize>(
        &mut self,
        x_data: &[Vector<F, N>],
        y: &[F],
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
        let alpha = self.alpha;
        let two = F::ONE + F::ONE;
        let f = move |pv: &Vector<F, M>| -> F {
            let mut loss = F::ZERO;
            for (x, &yi) in x_data.iter().zip(y.iter()) {
                let mut pred = pv[N];
                for i in 0..N {
                    pred = pred + pv[i] * x[i];
                }
                let e = pred - yi;
                loss = loss + e * e;
            }
            // L2 regularizer (alpha * ||w||^2)
            let mut rw = F::ZERO;
            for i in 0..N {
                rw = rw + pv[i] * pv[i];
            }
            (loss / m) + alpha * rw
        };

        let g = move |pv: &Vector<F, M>| -> Vector<F, M> {
            let mut grads: [F; M] = [F::ZERO; M];
            for (x, &yi) in x_data.iter().zip(y.iter()) {
                let mut pred = pv[N];
                for i in 0..N {
                    pred = pred + pv[i] * x[i];
                }
                let e = pred - yi;
                for i in 0..N {
                    grads[i] = grads[i] + two * e * x[i];
                }
                grads[N] = grads[N] + two * e;
            }
            for i in 0..N {
                grads[i] = grads[i] / m + two * alpha * pv[i];
            }
            grads[N] = grads[N] / m;
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
    fn ridge_handles_overfit() {
        // small dataset with noise
        let xs: [Vector<f64, 1>; 6] = [
            Vector::new([0.0]),
            Vector::new([1.0]),
            Vector::new([2.0]),
            Vector::new([3.0]),
            Vector::new([4.0]),
            Vector::new([5.0]),
        ];
        let ys: [f64; 6] = [1.0, 3.1, 4.9, 6.95, 9.05, 11.1];

        let mut model = RidgeRegression::<f64, 1> {
            weights: Vector::new([0.0]),
            bias: 0.0,
            lr: 0.05,
            alpha: 0.1,
        };
        for _ in 0..3000 {
            model.train_epoch(&xs, &ys);
        }

        let p = model.predict(&xs[3]);
        assert!((p - 6.95).abs() < 0.2, "pred={}", p);
    }
}
