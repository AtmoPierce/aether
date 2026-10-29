use aether_core::math::{Matrix, Vector};
use aether_core::numerical_methods::solvers::lu::*;
use aether_core::numerical_methods::solvers::newton::*;
use num_traits::{Float, One, Zero};
use core::fmt::Debug;



/// N = number of coefficients (degree = N-1)
/// K = number of alternation points / abscissas (must be N+1)
pub struct RemezApprox<T: Float, const N: usize, const K: usize> {
    pub coefficients: [T; N],      // monomial basis c0..c_{N-1}
    pub lower_bound: T,
    pub upper_bound: T,
    pub max_error_approximation: T, // |E|
}

impl<T: Float + Debug, const N: usize, const K: usize> RemezApprox<T, N, K> {
    pub fn new(lower: T, upper: T) -> Self {
        debug_assert!(K == N + 1, "Remez requires K = N + 1 (n+2 points)");
        Self {
            coefficients: [T::zero(); N],
            lower_bound: lower,
            upper_bound: upper,
            max_error_approximation: T::zero(),
        }
    }

    /// Analytic Remez (needs f and f'); single-point exchange per iteration.
    pub fn solve(&mut self, f: fn(T) -> T, df: fn(T) -> T) {
        const MAX_ITERS: usize = 10000;              // safeguard
        let tol = T::from(1e-7).unwrap();
        let mut xs = self.init_cheb_extrema();     // K nodes in [a,b], ascending
        let mut newton = NewtonRaphson::new(tol, 50);
        // match the paper’s sign pattern: -E at x0, then + - + ...
        let start_neg = true;

        for _ in 0..MAX_ITERS {
            let E = self.solve_coeffs_and_E(&xs, f, start_neg);

            let mut worst_x = xs[0];
            let mut worst_e = T::zero();

            self.scan_alternating_gaps(&xs, start_neg, |left, right| {
                // g(x) = P'(x) - f'(x);
                let g  = |x: T| self.eval_poly_deriv(x) - df(x);
                let h0 = T::from(1e-5).unwrap();
                let gp = |x: T| {
                    let hh = h0 * (T::one() + x.abs());
                    (g(x + hh) - g(x - hh)) / (T::from(2).unwrap() * hh)
                };

                let mid = (left + right) / T::from(2).unwrap();
                if let Some(xstar) = newton.solve(mid, g, gp).ok() {
                    if xstar > left && xstar < right {
                        let e = self.eval_poly(xstar) - f(xstar);
                        if e.abs() > worst_e.abs() {
                            worst_e = e;
                            worst_x = xstar;
                        }
                    }
                }
            });

            let delta = T::from(1e-3).unwrap();
            if self.alternates_close_to_E(&xs, f, E, delta)
                && worst_e.abs() <= (T::one() + delta) * E.abs()
            {
                self.max_error_approximation = E.abs();
                break;
            }

            self.exchange_in_place(&mut xs, worst_x, f);
        }
    }

    #[inline]
    pub fn eval_poly(&self, x: T) -> T {
        // Convert real x into normalized t
        let a = self.lower_bound;
        let b = self.upper_bound;
        let t = ((T::one()+T::one()) * x - (a + b)) / (b - a);

        // Horner evaluation in t (not x)
        let coeffs = &self.coefficients;
        let mut acc = T::zero();
        for &c in coeffs.iter().rev() {
            acc = acc * t + c;
        }
        acc
    }

    #[inline]
    fn eval_poly_deriv(&self, x: T) -> T {
        let mut acc = T::zero();
        for k in (1..N).rev() {
            let ck = self.coefficients[k] * T::from(k).unwrap();
            acc = acc * x + ck;
        }
        acc
    }

    // ----------------- core linear solve -----------------
    /// Solve the (KxK) system:
    ///   sum_{j=0..N-1} c_j x_i^j + s_i * E = f(x_i),  s_i alternates with start_neg.
    #[allow(non_snake_case)]
    fn solve_coeffs_and_E(&mut self, xs: &[T; K], f: fn(T) -> T, start_neg: bool) -> T {
        let mut A = Matrix::<T, K, K>::zeros();
        let mut b = Vector::<T, K>::zeros();

        for i in 0..K {
            // monomial columns
            let mut xpow = T::one();
            for j in 0..N {
                A[(i, j)] = xpow;
                xpow = xpow * xs[i];
            }
            // E column
            let base = if (i & 1) == 0 { T::one() } else { -T::one() };
            A[(i, K - 1)] = if start_neg { -base } else { base };

            b[i] = f(xs[i]);
        }

        let sol: Vector<T, K> = A.solve(&b).expect("Remez system solve failed");
        for j in 0..N { self.coefficients[j] = sol[j]; }
        sol[K - 1] // E
    }

    fn init_cheb_extrema(&self) -> [T; K] {
        let a = self.lower_bound;
        let b = self.upper_bound;
        let m = (a + b) / T::from(2).unwrap();
        let r = (b - a) / T::from(2).unwrap();
        let pi = T::from(core::f64::consts::PI).unwrap();

        let mut xs = [T::zero(); K];
        // theta from PI down to 0 -> cos in [-1..+1] -> x in [a..b], ascending (like paper)
        for i in 0..K {
            let theta = pi * T::from(K - 1 - i).unwrap() / T::from(K - 1).unwrap();
            xs[i] = m + r * T::cos(theta);
        }
        xs
    }

    #[allow(non_snake_case)]
    fn alternates_close_to_E(
        &self,
        xs: &[T; K],
        f: fn(T) -> T,
        E: T,
        rel_tol: T,
    ) -> bool {
        let denom = E.abs().max(T::one());
        for i in 0..K {
            let e = self.eval_poly(xs[i]) - f(xs[i]);
            let want = if (i & 1) == 0 { E } else { -E };
            if (e - want).abs() > rel_tol * denom { return false; }
        }
        true
    }

    fn scan_alternating_gaps<F: FnMut(T, T)>(&self, xs: &[T; K], start_neg: bool, mut cb: F) {
        let want_even = start_neg;
        let mut done = 0usize;
        for i in 0..(K - 1) {
            if ((i & 1) == 0) != want_even { continue; }
            if done >= K - 2 { break; }
            cb(xs[i], xs[i + 1]);
            done += 1;
        }
    }

    fn exchange_in_place(&self, xs: &mut [T; K], x_new: T, f: fn(T) -> T) {
        // locate insertion index
        let mut idx = 0;
        while idx < K && xs[idx] < x_new { idx += 1; }
        let left  = if idx == 0 { 0 } else { idx - 1 };
        let right = if idx >= K { K - 1 } else { idx };

        // never drop endpoints
        let drop_candidate = if left == 0 { right }
                             else if right == K - 1 { left }
                             else {
                                 let e_left  = (self.eval_poly(xs[left])  - f(xs[left])).abs();
                                 let e_right = (self.eval_poly(xs[right]) - f(xs[right])).abs();
                                 if e_left <= e_right { left } else { right }
                             };

        let mut put = drop_candidate;
        xs[put] = x_new;

        while put > 0 && xs[put] < xs[put - 1] { xs.swap(put, put - 1); put -= 1; }
        while put + 1 < K && xs[put] > xs[put + 1] { xs.swap(put, put + 1); put += 1; }
    }
}

#[derive(Clone, Copy)]
pub struct CycleReport<T: Float, const N: usize, const K: usize> {
    pub xs: [T; K],         // x_i (input nodes)
    pub coeffs: [T; N],     // b_i
    pub m: [T; K],          // next extrema (endpoints pinned)
    pub E: T,               // common error magnitude/sign
}

impl<T: Float + core::fmt::Display, const N: usize, const K: usize> CycleReport<T, N, K> {
    pub fn print(&self, cycle: usize) {
        println!("Cycle {}", cycle);
        println!("E = {:+.10}", self.E);
        println!("{:<3} {:>14} {:>16} {:>16} {:>16}", "i", "x_i", "b_i", "m_i", "E_i");
        for i in 0..K {
            let bi = if i < N { format!("{:+.10}", self.coeffs[i]) } else { String::new() };
            let mi = format!("{:+.10}", self.m[i]);
            let Ei = if (i & 1) == 0 { self.E } else { -self.E };
            println!("{:<3} {:>14.10} {:>16} {:>16} {:>16.10}", i, self.xs[i], bi, mi, Ei);
        }
        println!();
    }
}

impl<T: Float + Debug, const N: usize, const K: usize> RemezApprox<T, N, K> {
    // recreate paper's output table
    pub fn do_cycle(&mut self, xs: &[T; K], f: fn(T) -> T, df: fn(T) -> T, start_neg: bool)
        -> CycleReport<T, N, K>
    {
        let E = self.solve_coeffs_and_E(xs, f, start_neg);

        let a = self.lower_bound;
        let b = self.upper_bound;
        let mut m = [T::zero(); K];
        m[0] = a;
        m[K - 1] = b;

        let tol = T::from(1e-7).unwrap();
        let mut newton = NewtonRaphson::new(tol, 50);

        self.scan_alternating_gaps(xs, start_neg, |left, right| {
            let g  = |x: T| self.eval_poly_deriv(x) - df(x);
            let h0 = T::from(1e-5).unwrap();
            let gp = |x: T| {
                let hh = h0 * (T::one() + x.abs());
                (g(x + hh) - g(x - hh)) / (T::from(2).unwrap() * hh)
            };

            let mid = (left + right) / T::from(2).unwrap();
            let xstar = newton.solve(mid, g, gp).ok()
                              .filter(|&xr| xr > left && xr < right)
                              .unwrap_or(mid);

            let mut j = 1usize;
            while j < K-1 && m[j] != T::zero() { j += 1; }
            if j < K-1 { m[j] = xstar; }
        });
        CycleReport { xs: *xs, coeffs: self.coefficients, m, E }
    }
}

pub fn max_error_on_grid<const N: usize, const K: usize>(
    rz: &RemezApprox<f64, N, K>,
    f: fn(f64) -> f64,
    a: f64,
    b: f64,
    samples: usize,
) -> f64 {
    let mut max_e = 0.0;
    for i in 0..=samples {
        let t = i as f64 / samples as f64;
        let x = a + (b - a) * t;
        let e = (rz.eval_poly(x) - f(x)).abs();
        if e > max_e { max_e = e; }
    }
    max_e
}
