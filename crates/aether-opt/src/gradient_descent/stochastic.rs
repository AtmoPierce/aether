#![cfg_attr(not(test), no_std)]

use aether_core::math::Vector;

// Use the shared RNG implementation from `crate::randomizers`.
use aether_rand::randomizers::XorShift64Star;

// ------------------ Stochastic Gradient Descent ------------------

#[derive(Debug, Clone, Copy)]
pub struct GDStochastic<const N: usize> {
    // Learning-rate schedule: alpha_k = lr0 / (1 + decay * k)
    pub lr0: f64,
    pub decay: f64,

    // Momentum (Polyak): x_{k+1} = x_k - a_k * g_k + mu * (x_k - x_{k-1})
    pub momentum: f64, // 0 disables

    // Mini-batch size (calls to the gradient sampler per iteration)
    pub batch_size: u32, // 1 = pure SGD

    // Stopping (EMA of grad-norm) and step-norm
    pub ema_beta: f64,     // e.g., 0.9
    pub tol_grad_ema: f64, // stop if ema(||g||) <= tol_grad_ema
    pub tol_step: f64,     // stop if ||x_{k+1}-x_k|| <= tol_step
    pub max_iters: u32,
    pub tol_f: f64,
    pub min_iters: u32,
    pub use_bias_corrected_ema: bool,

    // Projection box
    has_bounds: bool,
    min_bounds: [f64; N],
    max_bounds: [f64; N],

    // RNG
    rng: XorShift64Star,
}

impl<const N: usize> GDStochastic<N> {
    /// Defaults:
    /// - lr0 = 1e-2, decay = 1e-3 (gentle)
    /// - momentum = 0.0
    /// - batch_size = 1
    /// - ema_beta = 0.9, tol_grad_ema = 1e-4, tol_step = 1e-8
    /// - max_iters = 100_000
    /// - no bounds
    pub const fn new(seed: u64) -> Self {
        Self {
            lr0: 1e-2,
            decay: 1e-3,
            momentum: 0.0,
            batch_size: 1,
            ema_beta: 0.9,
            tol_grad_ema: 1e-4,
            tol_step: 1e-8,
            max_iters: 100_000,
            tol_f: 0.0,
            min_iters: 0,
            use_bias_corrected_ema: true,
            has_bounds: false,
            min_bounds: [f64::NEG_INFINITY; N],
            max_bounds: [f64::INFINITY; N],
            rng: XorShift64Star::new(seed),
        }
    }

    // ---- mutating config setters (chainable) ----
    #[inline]
    pub fn learning_rate(&mut self, lr0: f64) -> &mut Self {
        self.lr0 = lr0.max(0.0);
        self
    }
    #[inline]
    pub fn decay(&mut self, decay: f64) -> &mut Self {
        self.decay = if decay < 0.0 { 0.0 } else { decay };
        self
    }
    #[inline]
    pub fn momentum(&mut self, mu: f64) -> &mut Self {
        self.momentum = clamp(mu, 0.0, 0.999999);
        self
    }
    #[inline]
    pub fn batch_size(&mut self, bs: u32) -> &mut Self {
        self.batch_size = if bs == 0 { 1 } else { bs };
        self
    }
    #[inline]
    pub fn ema(&mut self, beta: f64) -> &mut Self {
        self.ema_beta = clamp(beta, 0.0, 0.999999999);
        self
    }
    #[inline]
    pub fn tolerances(&mut self, tol_grad_ema: f64, tol_step: f64) -> &mut Self {
        self.tol_grad_ema = tol_grad_ema.max(0.0);
        self.tol_step = tol_step.max(0.0);
        self
    }
    #[inline]
    pub fn max_iters(&mut self, iters: u32) -> &mut Self {
        self.max_iters = iters;
        self
    }
    #[inline]
    pub fn set_bounds(&mut self, min_bounds: [f64; N], max_bounds: [f64; N]) -> &mut Self {
        self.has_bounds = true;
        self.min_bounds = min_bounds;
        self.max_bounds = max_bounds;
        self
    }
    #[inline]
    pub fn clear_bounds(&mut self) -> &mut Self {
        self.has_bounds = false;
        self
    }

    // ---- utilities ----
    #[inline]
    fn project(&self, mut x: Vector<f64, N>) -> Vector<f64, N> {
        if self.has_bounds {
            for i in 0..N {
                let xi = x[i];
                let lo = self.min_bounds[i];
                let hi = self.max_bounds[i];
                x[i] = if xi < lo {
                    lo
                } else if xi > hi {
                    hi
                } else {
                    xi
                };
            }
        }
        x
    }

    #[inline]
    fn norm2(v: &Vector<f64, N>) -> f64 {
        let mut s: f64 = 0.0;
        for i in 0..N {
            s += v[i] * v[i];
        }
        s.sqrt()
    }

    #[inline]
    fn lr_at(&self, k: u32) -> f64 {
        self.lr0 / (1.0 + self.decay * (k as f64))
    }

    // setters
    #[inline]
    pub fn tol_f(&mut self, v: f64) -> &mut Self {
        self.tol_f = if v < 0.0 { 0.0 } else { v };
        self
    }
    #[inline]
    pub fn min_iters(&mut self, it: u32) -> &mut Self {
        self.min_iters = it;
        self
    }
    #[inline]
    pub fn bias_corrected_ema(&mut self, on: bool) -> &mut Self {
        self.use_bias_corrected_ema = on;
        self
    }

    /// Stochastic minimize.
    ///
    /// - `x0` initial point
    /// - `f` objective value (used for returning the final value; can be cheap/approx)
    /// - `sample_grad` returns a **stochastic** gradient sample at `x` given `rng`
    ///
    /// If `batch_size > 1`, we call `sample_grad` that many times and average.
    ///
    /// Returns: (x_best, f_best, iters_used, converged)
    pub fn minimize<F, G>(
        &mut self,
        mut x: Vector<f64, N>,
        f: F,
        sample_grad: G,
    ) -> (Vector<f64, N>, f64, u32, bool)
    where
        F: Fn(&Vector<f64, N>) -> f64,
        G: Fn(&Vector<f64, N>, &mut XorShift64Star) -> Vector<f64, N>,
    {
        let mut v = Vector::new([0.0; N]);
        let mut ema_g = 0.0;
        let mut fx = f(&x);

        for k in 0..self.max_iters {
            // --- batch grad ---
            let mut g = Vector::new([0.0; N]);
            let bs = self.batch_size.max(1);
            for _ in 0..bs {
                let gi = sample_grad(&x, &mut self.rng);
                for i in 0..N {
                    g[i] += gi[i];
                }
            }
            let inv_bs = 1.0 / (bs as f64);
            for i in 0..N {
                g[i] *= inv_bs;
            }

            // --- stopping checks ---
            let gnorm = Self::norm2(&g);
            ema_g = self.ema_beta * ema_g + (1.0 - self.ema_beta) * gnorm;

            // bias-corrected EMA (like Adam)
            let ema_hat = if self.use_bias_corrected_ema {
                let bpow = self.ema_beta.powi((k as i32) + 1);
                ema_g / (1.0 - bpow)
            } else {
                ema_g
            };

            if k >= self.min_iters {
                if ema_hat <= self.tol_grad_ema || fx <= self.tol_f {
                    return (x, fx, k, true);
                }
            }

            // --- step ---
            let alpha = self.lr_at(k);
            let mut x_next = x;
            for i in 0..N {
                x_next[i] = x[i] + alpha * (-g[i] + self.momentum * v[i]);
            }
            x_next = self.project(x_next);
            for i in 0..N {
                v[i] = x_next[i] - x[i];
            }

            let step_norm = Self::norm2(&v);
            x = x_next;
            fx = f(&x);

            if k >= self.min_iters && step_norm <= self.tol_step {
                return (x, fx, k + 1, true);
            }
        }
        (x, fx, self.max_iters, false)
    }
}

// -------------- helpers --------------
#[inline]
const fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    if x < lo {
        lo
    } else if x > hi {
        hi
    } else {
        x
    }
}
