use aether_core::math::{Matrix, Vector};
use aether_core::numerical_methods::solvers::lu::*;
use aether_core::numerical_methods::solvers::newton::*;
use aether_core::real::Real;
use core::fmt::Debug;

/// N = number of coefficients (degree = N-1)
/// K = number of alternation points / abscissas (must be N+1)
pub struct RemezApprox<T: Real, const N: usize, const K: usize> {
    pub coefficients: [T; N], // monomial basis c0..c_{N-1}
    pub lower_bound: T,
    pub upper_bound: T,
    pub max_error_approximation: T, // |E|
}

impl<T: Real + Debug, const N: usize, const K: usize> RemezApprox<T, N, K> {
    pub fn new(lower: T, upper: T) -> Self {
        debug_assert!(K == N + 1, "Remez requires K = N + 1 (n+2 points)");
        Self {
            coefficients: [T::ZERO; N],
            lower_bound: lower,
            upper_bound: upper,
            max_error_approximation: T::ZERO,
        }
    }

    /// Analytic Remez (needs f and f'); single-point exchange per iteration.
    pub fn solve(&mut self, f: fn(T) -> T, df: fn(T) -> T) {
        const MAX_AXB_ITERS: usize = 1000;
        const MAX_NEWRAPHSON_ITERS: usize = 1000;

        let tol = T::from_f64(1e-7);
        let delta = T::from_f64(1e-3);
        let h0 = T::from_f64(1e-5);
        let two = T::from_f64(2.0);

        let mut xs = self.init_cheb_extrema();
        let newton = NewtonRaphson::new(tol, MAX_NEWRAPHSON_ITERS);

        // Match the paper’s sign pattern: -E at x0, then +E, -E, ...
        let start_neg = true;
        for _iter in 0..MAX_AXB_ITERS {
            let e_signed = self.solve_coeffs_and_E(&xs, f, start_neg);

            let mut worst_e = self.eval_poly(self.lower_bound) - f(self.lower_bound);

            // Check upper endpoint.
            {
                let x = self.upper_bound;
                let e = self.eval_poly(x) - f(x);

                if e.abs() > worst_e.abs() {
                    worst_e = e;
                }
            }

            // Check current support points.
            for &x in xs.iter() {
                let e = self.eval_poly(x) - f(x);

                if e.abs() > worst_e.abs() {
                    worst_e = e;
                }
            }

            // Check interior stationary points of the error in each support gap.
            self.scan_gaps(&xs, |left, right| {
                let g = |x: T| self.eval_poly_deriv(x) - df(x);

                let gp = |x: T| {
                    let hh = (T::ONE + x.abs()) * h0;
                    (g(x + hh) - g(x - hh)) / (two * hh)
                };

                let mid = (left + right) / two;

                if let Some(xstar) = newton
                    .solve(mid, g, gp)
                    .ok()
                    .filter(|&xr| xr > left && xr < right)
                {
                    let e = self.eval_poly(xstar) - f(xstar);

                    if e.abs() > worst_e.abs() {
                        worst_e = e;
                    }
                }
            });

            let e_abs = e_signed.abs();

            // Do not use max(|E|, 1). For typical approximations |E| << 1,
            // so using 1.0 makes the convergence test effectively meaningless.
            let floor = T::from_f64(1e-14);
            let denom = if e_abs > floor { e_abs } else { floor };

            if self.alternates_close_to_E(&xs, f, e_signed, delta, start_neg)
                && worst_e.abs() <= (T::ONE + delta) * denom
            {
                self.max_error_approximation = e_abs;
                break;
            }

            let prev_xs = xs;
            xs = self.next_extrema(&xs, f, df, tol, h0, MAX_NEWRAPHSON_ITERS);

            let mut moved = false;
            for i in 0..K {
                if xs[i] != prev_xs[i] {
                    moved = true;
                    break;
                }
            }

            if !moved {
                self.max_error_approximation = e_abs;
                break;
            }
        }
    }

    #[inline]
    pub fn eval_poly(&self, x: T) -> T {
        // Horner form
        let mut acc = T::ZERO;
        for &c in self.coefficients.iter().rev() {
            acc = acc * x + c;
        }
        acc
    }

    #[inline]
    fn eval_poly_deriv(&self, x: T) -> T {
        let mut acc = T::ZERO;
        for k in (1..N).rev() {
            let k_t = T::from_usize(k);
            let ck = self.coefficients[k] * k_t;
            acc = acc * x + ck;
        }
        acc
    }

    // ----------------- core linear solve -----------------
    /// Solve the (KxK) system:
    ///   sum_{j=0..N-1} c_j x_i^j + s_i * E = f(x_i),  s_i alternates with start_neg.
    #[allow(non_snake_case)]
    fn solve_coeffs_and_E(&mut self, xs: &[T; K], f: fn(T) -> T, start_neg: bool) -> T {
        if !strictly_increasing(xs) {
            panic!("Remez support points are not strictly increasing: {:?}", xs);
        }

        // Build A in Chebyshev basis on t = map(x)->[-1,1]
        let mut A = Matrix::<T, K, K>::new([[T::ZERO; K]; K]);
        let mut b = Vector::<T, K>::new([T::ZERO; K]);

        let two = T::from_f64(2.0);

        for i in 0..K {
            let x = xs[i];
            let t = self.to_t(x);

            // Fill columns 0..N-1 with T_j(t)
            // T0=1, T1=t, T_{j+1}=2 t T_j - T_{j-1}
            let mut Tjm1 = T::ONE; // T0
            A[(i, 0)] = Tjm1;

            if N > 1 {
                let mut Tj = t; // T1
                A[(i, 1)] = Tj;

                for j in 2..N {
                    let Tjp1 = two * t * Tj - Tjm1;
                    A[(i, j)] = Tjp1;
                    Tjm1 = Tj;
                    Tj = Tjp1;
                }
            }

            // E column (same as before)
            let base = if (i & 1) == 0 { T::ONE } else { -T::ONE };
            A[(i, K - 1)] = if start_neg { -base } else { base };

            b[i] = f(x);
        }

        let sol: Vector<T, K> = A
            .solve(&b)
            .unwrap_or_else(|| panic!("Remez system solve failed with xs={:?}", xs));

        // Chebyshev coeffs a_k are sol[0..N)
        let mut a_cheb = [T::ZERO; N];
        for j in 0..N {
            a_cheb[j] = sol[j];
        }

        // Convert Chebyshev series in t -> power series in t: p_t(t)=sum of( p_t[j] t^j)
        let p_t = cheb_to_power::<T, N>(&a_cheb);

        // Compose t = (x - m)/r to get power series in x (monomials),
        // so your existing eval/export remains unchanged.
        let (m, r) = self.m_r();
        let p_x = power_affine_compose::<T, N>(&p_t, m, r);

        self.coefficients = p_x;

        sol[K - 1] // E
    }

    fn init_cheb_extrema(&self) -> [T; K] {
        let a = self.lower_bound;
        let b = self.upper_bound;

        let two = T::from_f64(2.0);

        let m = (a + b) / two;
        let r = (b - a) / two;
        let pi = T::PI;

        let mut xs = [T::ZERO; K];

        // theta from PI down to 0 -> cos in [-1..+1] -> x in [a..b], ascending
        for i in 0..K {
            let num = T::from_usize(K - 1 - i);
            let den = T::from_usize(K - 1);
            let theta = pi * num / den;
            xs[i] = m + r * theta.cos();
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
        start_neg: bool,
    ) -> bool {
        let e_abs = E.abs();
        let floor = T::from_f64(1e-14);
        let denom = if e_abs > floor { e_abs } else { floor };
        for i in 0..K {
            let e = self.eval_poly(xs[i]) - f(xs[i]);

            // sign convention matches solve_coeffs_and_E:
            //   want = (-base * E) if start_neg else (base * E)
            let base = if (i & 1) == 0 { T::ONE } else { -T::ONE };
            let want = if start_neg { -base * E } else { base * E };

            if (e - want).abs() > rel_tol * denom {
                return false;
            }
        }
        true
    }

    fn scan_gaps<F: FnMut(T, T)>(&self, xs: &[T; K], mut cb: F) {
        for i in 0..(K - 1) {
            cb(xs[i], xs[i + 1]);
        }
    }

    fn next_extrema(
        &self,
        xs: &[T; K],
        f: fn(T) -> T,
        df: fn(T) -> T,
        tol: T,
        h0: T,
        max_newraphson_iters: usize,
    ) -> [T; K] {
        let two = T::from_f64(2.0);
        let mut drop_last = [T::ZERO; K];
        let mut drop_first = [T::ZERO; K];
        drop_last[0] = self.lower_bound;
        drop_first[K - 1] = self.upper_bound;

        let newton = NewtonRaphson::new(tol, max_newraphson_iters);
        let mut j = 0usize;

        self.scan_gaps(xs, |left, right| {
            if j >= K - 1 {
                return;
            }

            let g = |x: T| self.eval_poly_deriv(x) - df(x);
            let gp = |x: T| {
                let hh = (T::ONE + x.abs()) * h0;
                (g(x + hh) - g(x - hh)) / (two * hh)
            };

            let mid = (left + right) / two;
            let xstar = newton
                .solve(mid, g, gp)
                .ok()
                .filter(|&xr| xr > left && xr < right)
                .unwrap_or(mid);

            drop_last[j + 1] = xstar;
            drop_first[j] = xstar;
            j += 1;
        });

        if self.alternation_score(&drop_last, f) >= self.alternation_score(&drop_first, f) {
            drop_last
        } else {
            drop_first
        }
    }

    fn exchange_in_place(&self, xs: &mut [T; K], x_new: T, f: fn(T) -> T) {
        let two = T::from_f64(2.0);
        let mut x_new = x_new;

        // Keep inserted point strictly inside the support interval.
        if x_new <= xs[0] {
            x_new = (xs[0] + xs[1]) / two;
        } else if x_new >= xs[K - 1] {
            x_new = (xs[K - 2] + xs[K - 1]) / two;
        } else {
            // Avoid replacing a support point with itself.
            for i in 0..K {
                if x_new == xs[i] {
                    x_new = if i == 0 {
                        (xs[0] + xs[1]) / two
                    } else if i == K - 1 {
                        (xs[K - 2] + xs[K - 1]) / two
                    } else {
                        let left_gap = xs[i] - xs[i - 1];
                        let right_gap = xs[i + 1] - xs[i];

                        if left_gap >= right_gap {
                            (xs[i - 1] + xs[i]) / two
                        } else {
                            (xs[i] + xs[i + 1]) / two
                        }
                    };
                    break;
                }
            }
        }

        let e_new = self.eval_poly(x_new) - f(x_new);
        let s_new = signum_nonzero(e_new);

        // Locate insertion position in sorted xs.
        let mut idx = 0usize;
        while idx < K && xs[idx] < x_new {
            idx += 1;
        }

        // Candidate 1: drop left neighbor.
        let drop_left = if idx > 0 { Some(idx - 1) } else { None };

        // Candidate 2: drop right neighbor.
        let drop_right = if idx < K { Some(idx) } else { None };

        let mut best_drop: Option<usize> = None;
        let mut best_score = T::from_f64(-1.0);

        for drop in [drop_left, drop_right].into_iter().flatten() {
            // Never drop endpoints.
            if drop == 0 || drop == K - 1 {
                continue;
            }

            let mut trial = *xs;
            trial[drop] = x_new;
            sort_fixed(&mut trial);

            if !strictly_increasing(&trial) {
                continue;
            }

            let score = self.alternation_score(&trial, f);

            if score > best_score {
                best_score = score;
                best_drop = Some(drop);
            }
        }

        // If the sign-aware candidates are unavailable, fall back to endpoint-safe
        // local replacement.
        let drop = if let Some(drop) = best_drop {
            drop
        } else if idx == 0 {
            1
        } else if idx >= K {
            K - 2
        } else {
            let left = idx - 1;
            let right = idx;

            if left == 0 {
                right
            } else if right == K - 1 {
                left
            } else {
                let e_left = self.eval_poly(xs[left]) - f(xs[left]);
                let e_right = self.eval_poly(xs[right]) - f(xs[right]);

                let s_left = signum_nonzero(e_left);
                let s_right = signum_nonzero(e_right);

                // Prefer dropping a neighbor with the same sign as the new point.
                if s_left == s_new && s_right != s_new {
                    left
                } else if s_right == s_new && s_left != s_new {
                    right
                } else {
                    // Tie-break by dropping the smaller absolute error.
                    if e_left.abs() <= e_right.abs() {
                        left
                    } else {
                        right
                    }
                }
            }
        };

        xs[drop] = x_new;
        sort_fixed(xs);
    }

    fn alternation_score(&self, xs: &[T; K], f: fn(T) -> T) -> T {
        let mut sign_score = T::ZERO;
        let mut mag_score = T::ZERO;

        let mut prev_sign = signum_nonzero(self.eval_poly(xs[0]) - f(xs[0]));

        for i in 1..K {
            let e = self.eval_poly(xs[i]) - f(xs[i]);
            let sign = signum_nonzero(e);

            if sign != prev_sign {
                sign_score = sign_score + T::ONE;
            }

            mag_score = mag_score + e.abs();
            prev_sign = sign;
        }

        // Sign alternation dominates; magnitude is only a tie-breaker.
        sign_score * T::from_f64(1.0e6) + mag_score
    }
}

#[derive(Clone, Copy)]
pub struct CycleReport<T: Real, const N: usize, const K: usize> {
    pub xs: [T; K],     // x_i (input nodes)
    pub coeffs: [T; N], // b_i
    pub m: [T; K],      // next extrema (endpoints pinned)
    pub E: T,           // common error magnitude/sign
    pub start_neg: bool,
}

impl<T: Real + core::fmt::Display, const N: usize, const K: usize> CycleReport<T, N, K> {
    pub fn print(&self, cycle: usize) {
        println!("Cycle {}", cycle);
        println!("E = {:+.10}", self.E);
        println!("start_neg = {}", self.start_neg);
        println!(
            "{:<3} {:>14} {:>16} {:>16} {:>16}",
            "i", "x_i", "b_i", "m_i", "E_i"
        );

        for i in 0..K {
            let bi = if i < N {
                format!("{:+.10}", self.coeffs[i])
            } else {
                String::new()
            };
            let mi = format!("{:+.10}", self.m[i]);

            let base = if (i & 1) == 0 { T::ONE } else { -T::ONE };
            let Ei = if self.start_neg {
                -base * self.E
            } else {
                base * self.E
            };

            println!(
                "{:<3} {:>14.10} {:>16} {:>16} {:>16.10}",
                i, self.xs[i], bi, mi, Ei
            );
        }
        println!();
    }
}

impl<T: Real + Debug, const N: usize, const K: usize> RemezApprox<T, N, K> {
    // recreate paper's output table
    pub fn do_cycle(
        &mut self,
        xs: &[T; K],
        f: fn(T) -> T,
        df: fn(T) -> T,
        start_neg: bool,
    ) -> CycleReport<T, N, K> {
        let E = self.solve_coeffs_and_E(xs, f, start_neg);

        let tol = T::from_f64(1e-7);
        let h0 = T::from_f64(1e-5);
        let m = self.next_extrema(xs, f, df, tol, h0, 50);

        CycleReport {
            xs: *xs,
            coeffs: self.coefficients,
            m,
            E,
            start_neg,
        }
    }
    #[inline]
    fn to_t(&self, x: T) -> T {
        // map x in [a,b] -> t in [-1,1]
        let a = self.lower_bound;
        let b = self.upper_bound;
        let two = T::from_f64(2.0);
        (two * x - (a + b)) / (b - a)
    }

    #[inline]
    fn m_r(&self) -> (T, T) {
        // x = m + r*t
        let two = T::from_f64(2.0);
        let m = (self.lower_bound + self.upper_bound) / two;
        let r = (self.upper_bound - self.lower_bound) / two;
        (m, r)
    }
}

pub fn max_error_on_grid<T, const N: usize, const K: usize>(
    rz: &RemezApprox<T, N, K>,
    f: fn(T) -> T,
    a: T,
    b: T,
    samples: usize,
) -> T
where
    T: Real,
{
    let mut max_e = T::ZERO;
    for i in 0..=samples {
        let t = T::from_usize(i) / T::from_usize(samples);
        let x = a + (b - a) * t;
        let e = (rz.eval_poly(x) - f(x)).abs();
        if e > max_e {
            max_e = e;
        }
    }
    max_e
}
#[inline]
fn poly_shift_mul_t<T: Real, const N: usize>(p: &[T; N]) -> [T; N] {
    // returns t * p(t)
    let mut out = [T::ZERO; N];
    for j in (0..(N - 1)).rev() {
        out[j + 1] = p[j];
    }
    out
}

fn cheb_to_power<T: Real, const N: usize>(a: &[T; N]) -> [T; N] {
    // Build Σ a_k T_k(t) in power basis.
    // T0 = 1
    // T1 = t
    // Tk = 2 t T_{k-1} - T_{k-2}
    let two = T::from_f64(2.0);

    let mut sum = [T::ZERO; N];

    // T0
    let mut Tkm2 = [T::ZERO; N];
    Tkm2[0] = T::ONE;

    // add a0*T0
    sum[0] = sum[0] + a[0] * T::ONE;

    if N == 1 {
        return sum;
    }

    // T1
    let mut Tkm1 = [T::ZERO; N];
    Tkm1[1] = T::ONE;

    // add a1*T1
    for j in 0..N {
        sum[j] = sum[j] + a[1] * Tkm1[j];
    }

    // Tk for k>=2
    for k in 2..N {
        let tTkm1 = poly_shift_mul_t::<T, N>(&Tkm1);
        let mut Tk = [T::ZERO; N];
        for j in 0..N {
            Tk[j] = two * tTkm1[j] - Tkm2[j];
        }

        // sum += a[k]*Tk
        for j in 0..N {
            sum[j] = sum[j] + a[k] * Tk[j];
        }

        Tkm2 = Tkm1;
        Tkm1 = Tk;
    }

    sum
}

fn power_affine_compose<T: Real, const N: usize>(p_t: &[T; N], m: T, r: T) -> [T; N] {
    // Given p_t(t)=sum of p_t[j] t^j, compute p_x(x)=p_t((x-m)/r) in monomial x basis.
    //
    // Horner:
    //   y = p_t[n-1]
    //   for k=n-2..0: y = y * ((x-m)/r) + p_t[k]
    //
    // where y is a polynomial in x represented by coefficients [c0..c_{n-1}]
    let inv_r = T::ONE / r;

    let mut y = [T::ZERO; N];
    y[0] = p_t[N - 1];

    for kk in (0..(N - 1)).rev() {
        // y = y * (x - m)
        // (1) y*x : shift
        let mut yx = [T::ZERO; N];
        for j in (0..(N - 1)).rev() {
            yx[j + 1] = y[j];
        }

        // (2) yx - m*y
        for j in 0..N {
            yx[j] = yx[j] - m * y[j];
        }

        // (3) scale by 1/r
        for j in 0..N {
            yx[j] = yx[j] * inv_r;
        }

        // (4) add p_t[kk] to constant term
        yx[0] = yx[0] + p_t[kk];

        y = yx;
    }

    y
}

#[inline]
fn signum_nonzero<T: Real>(x: T) -> i8 {
    if x < T::ZERO {
        -1
    } else {
        1
    }
}

fn sort_fixed<T: Real, const K: usize>(xs: &mut [T; K]) {
    for i in 1..K {
        let mut j = i;
        while j > 0 && xs[j] < xs[j - 1] {
            xs.swap(j, j - 1);
            j -= 1;
        }
    }
}

fn strictly_increasing<T: Real, const K: usize>(xs: &[T; K]) -> bool {
    for i in 1..K {
        if !(xs[i - 1] < xs[i]) {
            return false;
        }
    }

    true
}
