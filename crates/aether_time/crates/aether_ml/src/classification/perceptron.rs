#![cfg(feature = "std")]
use aether::math::Vector;
use aether_opt::gradient_descent::GradientDescentGeneric;
use num_traits::{cast::cast, Float};

/// Classic perceptron for binary classification on fixed-size feature vectors.
#[derive(Clone, Copy, Debug)]
pub struct Perceptron<F: Float + Copy, const N: usize> {
    pub weights: Vector<F, N>,
    pub bias: F,
    pub lr: F,
}

impl<F: Float + Copy, const N: usize> Perceptron<F, N> {
    pub fn new() -> Self {
        Self {
            weights: Vector::new([F::zero(); N]),
            bias: F::zero(),
            lr: F::one(),
        }
    }

    pub fn predict_raw(&self, x: &Vector<F, N>) -> F {
        let mut s = self.bias;
        for i in 0..N {
            s = s + self.weights[i] * x[i];
        }
        s
    }
    pub fn predict(&self, x: &Vector<F, N>) -> i8 {
        if self.predict_raw(x) >= F::zero() {
            1
        } else {
            -1
        }
    }

    /// One epoch of perceptron updates using labels in {-1, +1}
    pub fn train_epoch(&mut self, x_data: &[Vector<F, N>], y: &[i8]) {
        for (x, &label) in x_data.iter().zip(y.iter()) {
            let pred = self.predict_raw(x);
            let lab_f: F = if label >= 0 { F::one() } else { -F::one() };
            if lab_f * pred <= F::zero() {
                // mistaken
                for i in 0..N {
                    self.weights[i] = self.weights[i] + self.lr * lab_f * x[i];
                }
                self.bias = self.bias + self.lr * lab_f;
            }
        }
    }

    /// Fit using GradientDescent over augmented param vector [w..., bias] with perceptron loss
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

        let f = move |pv: &aether::math::Vector<F, M>| -> F {
            // perceptron loss: sum(max(0, -y * (w·x + b)))
            let mut loss = F::zero();
            for (x, &yi) in x_data.iter().zip(y.iter()) {
                let mut s = pv[N];
                for i in 0..N {
                    s = s + pv[i] * x[i];
                }
                let yi_f: F = if yi >= 0 { F::one() } else { -F::one() };
                let val = -yi_f * s;
                if val > F::zero() {
                    loss = loss + val;
                }
            }
            loss / m
        };

        let g = move |pv: &aether::math::Vector<F, M>| -> aether::math::Vector<F, M> {
            let mut grads: [F; M] = [F::zero(); M];
            for (x, &yi) in x_data.iter().zip(y.iter()) {
                let mut s = pv[N];
                for i in 0..N {
                    s = s + pv[i] * x[i];
                }
                let yi_f: F = if yi >= 0 { F::one() } else { -F::one() };
                let val = -yi_f * s;
                if val > F::zero() {
                    for i in 0..N {
                        grads[i] = grads[i] + -yi_f * x[i];
                    }
                    grads[N] = grads[N] + -yi_f;
                }
            }
            for i in 0..M {
                grads[i] = grads[i] / m;
            }
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
    fn perceptron_learns_and() {
        let xs: [Vector<f64, 2>; 4] = [
            Vector::new([0.0, 0.0]),
            Vector::new([0.0, 1.0]),
            Vector::new([1.0, 0.0]),
            Vector::new([1.0, 1.0]),
        ];
        // AND mapping to {-1,+1}
        let ys: [i8; 4] = [-1, -1, -1, 1];
        let mut p = Perceptron::<f64, 2> {
            weights: Vector::new([0.0, 0.0]),
            bias: 0.0,
            lr: 0.5,
        };
        for _ in 0..100 {
            p.train_epoch(&xs, &ys);
        }
        let preds: Vec<i8> = xs.iter().map(|x| p.predict(x)).collect();
        assert_eq!(preds, vec![-1, -1, -1, 1]);
    }
}
