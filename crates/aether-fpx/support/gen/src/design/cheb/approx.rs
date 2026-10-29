pub struct Chebyshev {}

impl Chebyshev {
    pub fn identity(k: usize, x: f64) -> Result<f64, &'static str> {
        if x.is_nan() {
            return Err("x is NaN");
        }
        let kf = k as f64;
        if x.abs() <= 1.0 {
            Ok((kf * x.acos()).cos())
        } else if x > 1.0 {
            Ok((kf * x.acosh()).cosh())
        } else {
            // x < -1
            Ok((kf * (-x).acosh()).cosh() * (if k % 2 == 0 { 1.0 } else { -1.0 }))
            // since T_k(-x) = (-1)^k T_k(x)
        }
    }
}

pub fn chebyshev_t(n: usize, x: f64) -> f64 {
    // Trivial cases
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return x;
    }

    // Exact endpoints
    if x == 1.0 {
        return 1.0;
    }
    if x == -1.0 {
        return if n % 2 == 0 { 1.0 } else { -1.0 };
    }

    // Exact center avoids divide-by-zero in the series stepping
    if x == 0.0 {
        return if n % 2 == 0 {
            // n = 2k => T_n(0) = (-1)^k
            if (n / 2) % 2 == 0 {
                1.0
            } else {
                -1.0
            }
        } else {
            0.0
        };
    }

    // Series with coefficient ratio
    // A_0 = 2^(n-1)
    let mut a_m = pow2(n as i32 - 1);
    // x_pow = x^n
    let mut x_pow = powi(x, n as i32);
    let mut sum = a_m * x_pow;

    let m_max = n / 2;
    for m in 0..m_max {
        // A_{m+1} = -A_m * ((n-2m)(n-2m-1)) / (4 (m+1)(n-m-1))
        let nm2m = (n as f64) - 2.0 * (m as f64);
        let num = nm2m * (nm2m - 1.0);
        let den = 4.0 * ((m as f64) + 1.0) * ((n as f64) - (m as f64) - 1.0);
        a_m = -a_m * (num / den);
        // Step exponent down by 2 safely (we know x != 0 here)
        x_pow /= x * x;
        sum += a_m * x_pow;
    }
    sum
}

/// Cheap 2^k using exponent bits; handles typical degrees.
fn pow2(k: i32) -> f64 {
    if k >= -1022 && k <= 1023 {
        let bits = ((k + 1023) as u64) << 52;
        f64::from_bits(bits)
    } else if k > 1023 {
        let mut y = f64::from_bits((2046u64) << 52); // 2^1023
        for _ in 0..(k - 1023) {
            y *= 2.0;
        }
        y
    } else {
        // subnormal range: 2^k with k in [-1074, -1023)
        let shift = (k + 1074) as u32; // 0..=51
        let mant = 1u64 << shift;
        f64::from_bits(mant) // exponent=0
    }
}

/// x^k via fast powi
fn powi(mut x: f64, mut k: i32) -> f64 {
    if k == 0 {
        return 1.0;
    }
    let inv = k < 0;
    if inv {
        k = -k;
    }
    let mut a = 1.0;
    while k > 0 {
        if (k & 1) == 1 {
            a *= x;
        }
        x *= x;
        k >>= 1;
    }
    if inv {
        1.0 / a
    } else {
        a
    }
}

pub fn to_unit_interval(x: f64, a: f64, b: f64) -> f64 {
    (2.0 * x - (b + a)) / (b - a)
}
pub fn from_unit_interval(t: f64, a: f64, b: f64) -> f64 {
    0.5 * ((b - a) * t + (b + a))
}

// Stable Chebyshev Clenshaw: f(t) = sum of a[k] T_k(t), t in [-1,1]
#[inline]
pub fn cheb_eval<const N: usize>(a: &[f64; N], t: f64) -> f64 {
    if N == 0 {
        return 0.0;
    }
    if N == 1 {
        return a[0];
    }
    let mut b1 = 0.0;
    let mut b2 = 0.0;
    let mut k = N - 1;
    while k >= 1 {
        let bk = 2.0 * t * b1 - b2 + a[k];
        b2 = b1;
        b1 = bk;
        if k == 1 {
            break;
        }
        k -= 1;
    }
    t * b1 - b2 + 0.5 * a[0]
}

/// A Chebyshev approximation with compile-time length N on [a,b].
pub struct ChebApprox<const N: usize> {
    pub coeffs: [f64; N], // a[0..N-1] for T_0..T_{N-1}
    pub a: f64,
    pub b: f64,
}
impl<const N: usize> ChebApprox<N> {
    pub const fn new(coeffs: [f64; N], a: f64, b: f64) -> Self {
        Self { coeffs, a, b }
    }
    #[inline]
    pub fn eval_on_ab(&self, x: f64) -> f64 {
        let t = to_unit_interval(x, self.a, self.b);
        cheb_eval::<N>(&self.coeffs, t)
    }
}

/// Compute Chebyshev coefficients a_k for f on [a,b] using the DCT formula.
/// a_k = (2/N) S_{j=0}^{N-1} f(x_j) cos(k theta_j), where theta_j = PI(j+0.5)/N,
/// and x_j = from_unit_interval(cos theta_j).
pub fn cheb_coeffs<const N: usize>(mut f: impl FnMut(f64) -> f64, a: f64, b: f64) -> [f64; N] {
    let mut a_k = [0.0f64; N];
    // Precompute f at Chebyshev nodes mapped into [a,b]
    let mut fj = [0.0f64; N];
    for j in 0..N {
        let theta = core::f64::consts::PI * (j as f64 + 0.5) / (N as f64);
        let t = theta.cos(); // node in [-1,1]
        let x = from_unit_interval(t, a, b);
        fj[j] = f(x);
    }
    // DCT-II style sum
    let scale = 2.0 / (N as f64);
    for k in 0..N {
        let mut s = 0.0;
        for j in 0..N {
            let theta = core::f64::consts::PI * (j as f64 + 0.5) / (N as f64);
            s += fj[j] * (k as f64 * theta).cos();
        }
        a_k[k] = scale * s;
    }
    a_k
}
