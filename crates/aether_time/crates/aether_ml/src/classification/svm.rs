#![cfg(feature = "std")]
use aether::math::Vector;
use aether_opt::gradient_descent::GradientDescentGeneric;
use num_traits::{cast::cast, Float};

/// Linear SVM using hinge loss and L2 regularization trained with batch GD.
#[derive(Clone, Copy, Debug)]
pub struct LinearSvm<F: Float + Copy, const N: usize> {
    pub weights: Vector<F, N>,
    pub bias: F,
    pub lr: F,
    pub c: F, // regularization weight (higher = less regularization)
}

impl<F: Float + Copy, const N: usize> LinearSvm<F, N> {
    pub fn new() -> Self {
        Self {
            weights: Vector::new([F::zero(); N]),
            bias: F::zero(),
            lr: F::one(),
            c: F::one(),
        }
    }

    pub fn decision(&self, x: &Vector<F, N>) -> F {
        let mut s = self.bias;
        for i in 0..N {
            s = s + self.weights[i] * x[i];
        }
        s
    }

    /// One epoch of batch gradient descent on hinge loss
    /// y in {-1, +1}
    pub fn train_epoch(&mut self, x_data: &[Vector<F, N>], y: &[i8]) {
        let m_usize = x_data.len();
        if m_usize == 0 {
            return;
        }
        let m: F = cast(m_usize).unwrap();

        let mut grad_w = Vector::new([F::zero(); N]);
        let mut grad_b = F::zero();

        for (x, &yi) in x_data.iter().zip(y.iter()) {
            let y_f: F = if yi >= 0 { F::one() } else { -F::one() };
            let margin = y_f * self.decision(x);
            if margin >= F::one() {
                // loss has zero hinge; gradients are from regularizer only
                for i in 0..N {
                    grad_w[i] = grad_w[i] + F::zero();
                }
                grad_b = grad_b + F::zero();
            } else {
                // hinge active
                for i in 0..N {
                    grad_w[i] = grad_w[i] + -y_f * x[i];
                }
                grad_b = grad_b + -y_f;
            }
        }
        // apply average + regularization gradient (L2): grad = (1/m)*grad + (2/mC) * w
        for i in 0..N {
            let two = F::one() + F::one();
            let reg = (two / (self.c * m)) * self.weights[i];
            self.weights[i] = self.weights[i] - self.lr * ((grad_w[i] / m) + reg);
        }
        self.bias = self.bias - self.lr * (grad_b / m);
    }

    /// Fit using crate gradient descent over augmented parameters [w..., bias]
    pub fn fit_with_optimizer<const M: usize>(
        &mut self,
        x_data: &[Vector<F, N>],
        y: &[i8],
        mut opt: GradientDescentGeneric<F, M>,
    ) {
        assert_eq!(M, N + 1, "fit_with_optimizer: expected M == N+1");
        let m_usize = x_data.len();
        if m_usize == 0 {
            return;
        }
        let m: F = cast(m_usize).unwrap();
        let mut p_arr: [F; M] = [F::zero(); M];
        for i in 0..N {
            p_arr[i] = self.weights[i];
        }
        p_arr[N] = self.bias;
        let x0 = aether::math::Vector::new(p_arr);
        let c = self.c;

        let f = move |pv: &aether::math::Vector<F, M>| -> F {
            let mut loss = F::zero();
            for (x, &yi) in x_data.iter().zip(y.iter()) {
                let mut s = pv[N];
                for i in 0..N {
                    s = s + pv[i] * x[i];
                }
                let y_f: F = if yi >= 0 { F::one() } else { -F::one() };
                let margin = y_f * s;
                if margin < F::one() {
                    loss = loss + (F::one() - margin);
                }
            }
            // add L2 reg (1/(2C) * ||w||^2)
            let mut rw = F::zero();
            for i in 0..N {
                rw = rw + pv[i] * pv[i];
            }
            let two = F::one() + F::one();
            loss / m + (rw / (two * c))
        };

        let g = move |pv: &aether::math::Vector<F, M>| -> aether::math::Vector<F, M> {
            let mut grads: [F; M] = [F::zero(); M];
            for (x, &yi) in x_data.iter().zip(y.iter()) {
                let y_f: F = if yi >= 0 { F::one() } else { -F::one() };
                let mut s = pv[N];
                for i in 0..N {
                    s = s + pv[i] * x[i];
                }
                let margin = y_f * s;
                if margin < F::one() {
                    for i in 0..N {
                        grads[i] = grads[i] + -y_f * x[i];
                    }
                    grads[N] = grads[N] + -y_f;
                }
            }
            for i in 0..N {
                grads[i] = grads[i] / m + pv[i] / c;
            }
            grads[N] = grads[N] / m;
            aether::math::Vector::new(grads)
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
    use aether::math::Vector;

    #[test]
    fn svm_simple() {
        // Linearly separable
        let xs: [Vector<f64, 1>; 4] = [
            Vector::new([0.0]),
            Vector::new([1.0]),
            Vector::new([2.0]),
            Vector::new([3.0]),
        ];
        let ys: [i8; 4] = [-1, -1, 1, 1];
        let mut s = LinearSvm::<f64, 1> {
            weights: Vector::new([0.0]),
            bias: 0.0,
            lr: 0.1,
            c: 1.0,
        };
        for _ in 0..200 {
            s.train_epoch(&xs, &ys);
        }
        let preds: Vec<i8> = xs
            .iter()
            .map(|x| if s.decision(x) >= 0.0 { 1 } else { -1 })
            .collect();
        assert_eq!(preds, vec![-1, -1, 1, 1]);
    }
}
