#![cfg(feature = "std")]

// Requires std (parallelism via std::thread)
use std::thread;

use aether::math::Vector;
use aether_rand::randomizers::XorShift64Star;
use num_traits::{cast::cast, Float};

#[inline]
fn clamp<F: Float>(x: F, lo: F, hi: F) -> F {
    if x < lo {
        lo
    } else if x > hi {
        hi
    } else {
        x
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GDStochasticParallel<F: Float + Copy + Send, const N: usize> {
    // Learning-rate schedule α_k = lr0 / (1 + decay k)
    pub lr0: F,
    pub decay: F,
    // Momentum (Polyak) on parameter steps
    pub momentum: F, // 0 disables
    // Mini-batch size per iteration (total, across threads)
    pub batch_size: u32,
    // Stopping
    pub ema_beta: F, // EMA(||g||), e.g., 0.9–0.99
    pub tol_grad_ema: F,
    pub tol_step: F,
    pub tol_f: F,       // optional stop on objective value
    pub min_iters: u32, // patience before we allow stopping
    pub max_iters: u32,
    // Projection box
    has_bounds: bool,
    min_bounds: [F; N],
    max_bounds: [F; N],
    // Threads
    threads: Option<usize>, // None => auto
    // RNG base seed (per-iter/thread seeds derived from this)
    seed: u64,
}

impl<F: Float + Copy + Send, const N: usize> GDStochasticParallel<F, N> {
    pub fn new(seed: u64) -> Self {
        // initialize bounds arrays
        let mut min_bounds: [F; N] = [F::zero(); N];
        let mut max_bounds: [F; N] = [F::zero(); N];
        let neg_inf = -F::infinity();
        let pos_inf = F::infinity();
        for i in 0..N {
            min_bounds[i] = neg_inf;
            max_bounds[i] = pos_inf;
        }

        Self {
            lr0: cast(1e-2_f64).unwrap(),
            decay: cast(1e-3_f64).unwrap(),
            momentum: F::zero(),
            batch_size: 32,
            ema_beta: cast(0.95_f64).unwrap(),
            tol_grad_ema: cast(1e-4_f64).unwrap(),
            tol_step: cast(1e-8_f64).unwrap(),
            tol_f: F::zero(),
            min_iters: 0,
            max_iters: 100_000,
            has_bounds: false,
            min_bounds,
            max_bounds,
            threads: None,
            seed,
        }
    }

    // --------- config (chainable, mutating) ----------
    #[inline]
    pub fn learning_rate(&mut self, lr0: F) -> &mut Self {
        self.lr0 = if lr0 <= F::zero() { F::zero() } else { lr0 };
        self
    }
    #[inline]
    pub fn decay(&mut self, d: F) -> &mut Self {
        self.decay = if d < F::zero() { F::zero() } else { d };
        self
    }
    #[inline]
    pub fn momentum(&mut self, mu: F) -> &mut Self {
        self.momentum = clamp(mu, F::zero(), F::one());
        self
    }
    #[inline]
    pub fn batch_size(&mut self, bs: u32) -> &mut Self {
        self.batch_size = if bs == 0 { 1 } else { bs };
        self
    }
    #[inline]
    pub fn ema(&mut self, beta: F) -> &mut Self {
        self.ema_beta = clamp(beta, F::zero(), F::one());
        self
    }
    #[inline]
    pub fn tolerances(&mut self, tol_grad_ema: F, tol_step: F) -> &mut Self {
        self.tol_grad_ema = if tol_grad_ema <= F::zero() {
            F::zero()
        } else {
            tol_grad_ema
        };
        self.tol_step = if tol_step <= F::zero() {
            F::zero()
        } else {
            tol_step
        };
        self
    }
    #[inline]
    pub fn tol_f(&mut self, tf: F) -> &mut Self {
        self.tol_f = if tf < F::zero() { F::zero() } else { tf };
        self
    }
    #[inline]
    pub fn min_iters(&mut self, it: u32) -> &mut Self {
        self.min_iters = it;
        self
    }
    #[inline]
    pub fn max_iters(&mut self, it: u32) -> &mut Self {
        self.max_iters = it;
        self
    }
    #[inline]
    pub fn set_bounds(&mut self, min_bounds: [F; N], max_bounds: [F; N]) -> &mut Self {
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
    #[inline]
    pub fn threads(&mut self, n: usize) -> &mut Self {
        self.threads = Some(n.max(1));
        self
    }

    // -------------- internals ---------------
    #[inline]
    fn project(&self, mut x: Vector<F, N>) -> Vector<F, N> {
        if self.has_bounds {
            for i in 0..N {
                let (lo, hi) = (self.min_bounds[i], self.max_bounds[i]);
                let xi = x[i];
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
    fn norm2(v: &Vector<F, N>) -> F {
        let mut s = F::zero();
        for i in 0..N {
            s = s + v[i] * v[i];
        }
        s.sqrt()
    }

    #[inline]
    fn lr_at(&self, k: u32) -> F {
        // lr0 / (1 + decay * k)
        let one: F = cast::<u32, F>(1u32).unwrap();
        let k_f: F = cast::<u32, F>(k as u32).unwrap();
        self.lr0 / (one + self.decay * k_f)
    }

    // Deterministic per-thread seed for iteration k and thread t
    #[inline]
    fn thread_seed(base: u64, iter: u32, t: usize) -> u64 {
        base ^ ((iter as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15))
            ^ ((t as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9))
    }

    /// Parallel stochastic minimize.
    ///
    /// - `f(x)` returns objective value (cheap; used for return/optionally tol_f)
    /// - `sample_grad(x, rng)` returns **one** stochastic gradient sample
    ///
    /// Returns: (x_best, f_best, iters_used, converged)
    pub fn minimize<FuncF, G>(
        &self,
        mut x: Vector<F, N>,
        f: FuncF,
        sample_grad: G,
    ) -> (Vector<F, N>, F, u32, bool)
    where
        FuncF: Fn(&Vector<F, N>) -> F + Sync,
        G: Fn(&Vector<F, N>, &mut XorShift64Star) -> Vector<F, N> + Sync,
    {
        let threads = self
            .threads
            .unwrap_or_else(|| {
                thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(1)
            })
            .max(1);
        std::println!("Using {threads} threads for batch size {}", self.batch_size);
        let mut v = Vector::new([F::zero(); N]); // momentum buffer (Δx)
        let mut ema_g = F::zero();
        let mut fx = f(&x);

        for k in 0..self.max_iters {
            let bs_total = self.batch_size.max(1) as usize;
            let t_used = threads.min(bs_total);
            let base = self.seed;

            // Take a shared reference to the gradient sampler so threads don't move G.
            let sg = &sample_grad;

            // Spawn threads that each return (partial_sum, count).
            // Scoped join handles cannot escape the scope, so collect joined
            // results (owned) inside the scope into `partials`, which we
            // then reduce after the scope.
            let mut partials: Vec<(Vector<F, N>, usize)> = Vec::with_capacity(t_used);

            thread::scope(|scope| {
                // create handles inside the scope
                let mut handles = Vec::with_capacity(t_used);

                for t in 0..t_used {
                    let start = (bs_total * t) / t_used;
                    let end = (bs_total * (t + 1)) / t_used;
                    let seed_t = Self::thread_seed(base, k, t);
                    let x_local = x; // Copy by value into the thread
                    let sg_ref = sg; // &'scope G; Send because G: Sync

                    let handle = scope.spawn(move || {
                        let mut rng = XorShift64Star::new(seed_t);
                        let mut acc = Vector::new([F::zero(); N]);
                        for _ in start..end {
                            let gi = sg_ref(&x_local, &mut rng);
                            for i in 0..N {
                                acc[i] = acc[i] + gi[i];
                            }
                        }
                        (acc, end - start)
                    });

                    handles.push(handle);
                }

                // Join inside the scope and move owned results into `partials`.
                for h in handles {
                    let res = h.join().expect("thread panicked");
                    partials.push(res);
                }
            });

            // Reduce partials (sum)
            let mut g = Vector::new([F::zero(); N]);
            let mut total = 0usize;
            for (acc, cnt) in partials {
                total += cnt;
                for i in 0..N {
                    g[i] = g[i] + acc[i];
                }
            }
            // average gradient over total samples
            let inv_bs: F = cast::<u32, F>(1u32).unwrap() / cast::<u32, F>(total as u32).unwrap();
            for i in 0..N {
                g[i] = g[i] * inv_bs;
            }

            // EMA of grad norm (bias-corrected)
            let gnorm = Self::norm2(&g);
            let one: F = cast(1u8).unwrap();
            ema_g = self.ema_beta * ema_g + (one - self.ema_beta) * gnorm;
            let ema_hat = {
                let bpow = self.ema_beta.powi((k as i32) + 1);
                ema_g / (one - bpow)
            };

            if k >= self.min_iters {
                if ema_hat <= self.tol_grad_ema || (self.tol_f > F::zero() && fx <= self.tol_f) {
                    return (x, fx, k, true);
                }
            }

            // Parameter update: x_{k+1} = x_k - α g + μ v
            let alpha = self.lr_at(k);
            let mut x_next = x;
            for i in 0..N {
                let d_i = -g[i] + self.momentum * v[i];
                x_next[i] = x[i] + alpha * d_i;
            }
            x_next = self.project(x_next);

            // Momentum buffer v = Δx
            for i in 0..N {
                v[i] = x_next[i] - x[i];
            }

            // Step stop
            let step_norm = Self::norm2(&v);
            x = x_next;
            fx = f(&x);

            if k >= self.min_iters && step_norm <= self.tol_step {
                return (x, fx, k + 1, true);
            }
        }

        (x, fx, self.max_iters, F::zero() == F::zero())
    }
}
