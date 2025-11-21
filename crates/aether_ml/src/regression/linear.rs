#![cfg(feature = "std")]
use aether_core::math::Vector;
use aether_opt::gradient_descent::GradientDescentGeneric;
use aether_core::real::Real;

/// Ordinary Least Squares via batch gradient descent for fixed-size features
#[derive(Clone, Copy, Debug)]
pub struct LinearRegression<F: Real + Copy, const N: usize> {
    pub weights: Vector<F, N>,
    pub bias: F,
    pub lr: F,
}

impl<F: Real + Copy, const N: usize> LinearRegression<F, N> {
    pub fn new() -> Self {
        // Neutral default learning rate: 1 (callers/tests should set a smaller lr if desired)
        Self {
            weights: Vector::new([F::ZERO; N]),
            bias: F::ZERO,
            lr: F::ONE,
        }
    }

    pub fn predict(&self, x: &Vector<F, N>) -> F {
        let mut s = self.bias;
        for i in 0..N {
            s = s + self.weights[i] * x[i];
        }
        s
    }

    /// One epoch of batch gradient descent
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
            self.weights[i] = self.weights[i] - self.lr * (two * grad_w[i] / m);
        }
        self.bias = self.bias - self.lr * (two * grad_b / m);
    }

    /// Fit using the crate's GradientDescent optimizer. Packs weights and bias into an
    /// augmented parameter vector of size N+1: [w_0..w_{N-1}, bias].
    pub fn fit_with_optimizer<const M: usize>(
        &mut self,
        x_data: &[Vector<F, N>],
        y: &[F],
        opt: GradientDescentGeneric<F, M>,
    ) {
        assert_eq!(M, N + 1, "fit_with_optimizer: expected M == N+1");
        let m_usize = x_data.len();
        if m_usize == 0 {
            return;
        }
        let m: F = F::from_usize(m_usize);

        // initial param vector
        let mut p_arr: [F; M] = [F::ZERO; M];
        for i in 0..N {
            p_arr[i] = self.weights[i];
        }
        p_arr[N] = self.bias;
        let x0 = Vector::new(p_arr);
        let two = F::ONE + F::ONE;

        let f = |pv: &Vector<F, M>| -> F {
            // unpack
            let mut loss = F::ZERO;
            for (x, &yi) in x_data.iter().zip(y.iter()) {
                let mut pred = pv[N];
                for i in 0..N {
                    pred = pred + pv[i] * x[i];
                }
                let e = pred - yi;
                loss = loss + e * e;
            }
            loss / m
        };

        let g = |pv: &Vector<F, M>| -> Vector<F, M> {
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

pub type LinearRegressionF64<const N: usize> = LinearRegression<f64, N>;

#[cfg(test)]
mod tests {
    use super::*;
    use aether_core::math::Vector;

    #[test]
    fn fit_simple_line() {
        // y = 2*x + 1
        let xs: [Vector<f64, 1>; 5] = [
            Vector::new([0.0]),
            Vector::new([1.0]),
            Vector::new([2.0]),
            Vector::new([3.0]),
            Vector::new([4.0]),
        ];
        let ys: [f64; 5] = [1.0, 3.0, 5.0, 7.0, 9.0];

        let mut model = LinearRegression::<f64, 1> {
            weights: Vector::new([0.0]),
            bias: 0.0,
            lr: 0.1,
        };
        for _ in 0..2000 {
            model.train_epoch(&xs, &ys);
        }

        let p = model.predict(&xs[2]);
        assert!((p - 5.0).abs() < 1e-1, "pred={}", p);
    }
}
