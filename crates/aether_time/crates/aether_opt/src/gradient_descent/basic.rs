#![cfg_attr(not(test), no_std)]

use aether::math::Vector;
use num_traits::{cast::cast, Float};

/// No-alloc, no_std gradient descent with optional momentum, projection,
/// and Armijo backtracking. Generic over floating-point type F.
///
/// Strategy:
///   x_{k+1} = P[ x_k - α_k * g_k + μ * v_k ]
///   v_{k+1} = x_{k+1} - x_k        (Polyak momentum)
#[derive(Debug, Clone, Copy)]
pub struct GradientDescentGeneric<F: Float + Copy, const N: usize> {
    // Step-size strategy
    pub lr: F, // base learning rate (used when backtracking is off, or as initial guess when on)
    pub use_backtracking: bool,
    pub armijo_c1: F,     // typical ~1e-4
    pub backtrack_tau: F, // multiplicative shrink factor (0<tau<1), e.g. 0.5

    // Momentum (Polyak): x_{k+1} = x_k - α g_k + μ v_k
    pub momentum: F, // 0 disables momentum

    // Stopping criteria
    pub tol_grad: F, // stop if ||grad|| <= tol_grad
    pub tol_step: F, // or if ||Δx|| <= tol_step
    pub max_iters: u32,

    // Optional box constraints (projection)
    has_bounds: bool,
    min_bounds: [F; N],
    max_bounds: [F; N],
}

impl<F: Float + Copy, const N: usize> GradientDescentGeneric<F, N> {
    /// Sensible defaults (non-const because F is generic):
    ///  - lr = 1e-2
    ///  - Armijo backtracking on (c1 = 1e-4, tau = 0.5)
    ///  - momentum = 0 (off)
    ///  - tolerances = 1e-6
    ///  - max_iters = 10_000
    pub fn new() -> Self {
        // initialize arrays by filling
        let mut min_bounds: [F; N] = [F::zero(); N];
        let mut max_bounds: [F; N] = [F::zero(); N];
        let neg_inf = -F::infinity();
        let pos_inf = F::infinity();
        for i in 0..N {
            min_bounds[i] = neg_inf;
            max_bounds[i] = pos_inf;
        }

        Self {
            lr: cast(1e-2_f64).unwrap(),
            use_backtracking: true,
            armijo_c1: cast(1e-4_f64).unwrap(),
            backtrack_tau: cast(0.5_f64).unwrap(),
            momentum: F::zero(),
            tol_grad: cast(1e-6_f64).unwrap(),
            tol_step: cast(1e-6_f64).unwrap(),
            max_iters: 10_000,
            has_bounds: false,
            min_bounds,
            max_bounds,
        }
    }

    /// Disable/enable backtracking.
    pub fn backtracking(&mut self, on: bool) -> &mut Self {
        self.use_backtracking = on;
        self
    }
    pub fn armijo(&mut self, c1: F, tau: F) -> &mut Self {
        self.armijo_c1 = c1;
        self.backtrack_tau = if tau <= F::zero() {
            cast(0.5_f64).unwrap()
        } else if tau >= F::one() {
            cast(0.9_f64).unwrap()
        } else {
            tau
        };
        self
    }
    pub fn learning_rate(&mut self, lr: F) -> &mut Self {
        self.lr = if lr <= F::zero() { F::zero() } else { lr };
        self
    }
    pub fn momentum(&mut self, mu: F) -> &mut Self {
        let cap = cast(0.999999_f64).unwrap();
        self.momentum = if mu < F::zero() {
            F::zero()
        } else if mu > cap {
            cap
        } else {
            mu
        };
        self
    }
    pub fn tolerances(&mut self, tol_grad: F, tol_step: F) -> &mut Self {
        self.tol_grad = if tol_grad <= F::zero() {
            F::zero()
        } else {
            tol_grad
        };
        self.tol_step = if tol_step <= F::zero() {
            F::zero()
        } else {
            tol_step
        };
        self
    }
    pub fn max_iters(&mut self, iters: u32) -> &mut Self {
        self.max_iters = iters;
        self
    }

    /// Set per-axis [min,max] bounds; projection is applied every step.
    pub fn set_bounds(&mut self, min_bounds: [F; N], max_bounds: [F; N]) -> &mut Self {
        self.has_bounds = true;
        self.min_bounds = min_bounds;
        self.max_bounds = max_bounds;
        self
    }
    /// Clear bounds (disable projection).
    pub fn clear_bounds(&mut self) -> &mut Self {
        self.has_bounds = false;
        self
    }

    /// Project each component of `x` into [min_i, max_i] if bounds are enabled.
    #[inline]
    fn project(&self, mut x: Vector<F, N>) -> Vector<F, N> {
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

    /// Euclidean norm (L2) using only core ops.
    #[inline]
    fn norm2(v: &Vector<F, N>) -> F {
        let mut s = F::zero();
        for i in 0..N {
            s = s + v[i] * v[i];
        }
        s.sqrt()
    }

    /// One optimization run.
    ///
    /// `f`  — objective value at x
    /// `g`  — gradient at x (must return same-dimension vector)
    ///
    /// Returns: (x_best, f_best, iters_used, converged)
    pub fn minimize<FUN, GRAD>(
        &self,
        mut x: Vector<F, N>,
        f: FUN,
        g: GRAD,
    ) -> (Vector<F, N>, F, u32, bool)
    where
        FUN: Fn(&Vector<F, N>) -> F,
        GRAD: Fn(&Vector<F, N>) -> Vector<F, N>,
    {
        // Momentum state (v = x_k - x_{k-1}); initial 0
        let mut v = Vector::new([F::zero(); N]);

        let mut fx = f(&x);
        let mut converged = false;

        for k in 0..self.max_iters {
            let grad = g(&x);
            let grad_norm = Self::norm2(&grad);
            if grad_norm <= self.tol_grad {
                converged = true;
                return (x, fx, k, converged);
            }

            // Search direction: steepest descent with momentum
            // d = -grad + μ * v
            let mut d = Vector::new([F::zero(); N]);
            for i in 0..N {
                d[i] = -grad[i] + self.momentum * v[i];
            }

            let mut step = self.lr;
            let x_old = x;
            let fx_old = fx;

            if self.use_backtracking {
                // Armijo condition: f(x + a d) <= f(x) + c1 * a * grad·d
                let mut gd = F::zero();
                for i in 0..N {
                    gd = gd + grad[i] * d[i];
                }
                // If d accidentally isn’t a descent direction (gd >= 0),
                // fall back to pure -grad.
                if gd >= F::zero() {
                    for i in 0..N {
                        d[i] = -grad[i];
                    }
                    gd = -grad_norm * grad_norm;
                }

                loop {
                    let mut trial = x_old;
                    for i in 0..N {
                        trial[i] = x_old[i] + step * d[i];
                    }
                    trial = self.project(trial);
                    let f_trial = f(&trial);
                    let rhs = fx_old + self.armijo_c1 * step * gd;

                    if f_trial <= rhs {
                        x = trial;
                        fx = f_trial;
                        break;
                    }

                    step = step * self.backtrack_tau;
                    if step <= cast(1e-20_f64).unwrap() {
                        // give up this iteration; behave like tiny step
                        x = trial;
                        fx = f_trial;
                        break;
                    }
                }
            } else {
                // Fixed step
                for i in 0..N {
                    x[i] = x[i] + step * d[i];
                }
                x = self.project(x);
                fx = f(&x);
            }

            // Update momentum v = x_{k+1} - x_k
            for i in 0..N {
                v[i] = x[i] - x_old[i];
            }

            // Step-norm stopping
            let step_norm = Self::norm2(&v);
            if step_norm <= self.tol_step {
                converged = true;
                return (x, fx, k + 1, converged);
            }
        }

        (x, fx, self.max_iters, converged)
    }
}

// No compatibility aliases: users should instantiate `GradientDescentGeneric<F, N>` directly.
