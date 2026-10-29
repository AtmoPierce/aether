use aether_core::real::Real;

// ---- tiny RNG reused from no_std version ----
#[derive(Debug, Clone, Copy)]
pub struct XorShift64Star {
    state: u64,
}
impl XorShift64Star {
    pub const fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 0x9E3779B97F4A7C15 } else { seed },
        }
    }
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }

    #[inline]
    pub fn next_real<T: Real>(&mut self) -> T {
        const DEN: f64 = (1u64 << 53) as f64;
        let unit = ((self.next_u64() >> 11) as f64) / DEN;
        T::from_f64(unit)
    }

    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        self.next_real::<f64>()
    }

    /// Approximate standard normal via Box–Muller.
    pub fn normal01_real<T: Real>(&mut self) -> T {
        let mut u1 = self.next_real::<T>();
        if u1 <= T::ZERO {
            u1 = T::EPSILON;
        }
        let u2 = self.next_real::<T>();
        let two = T::ONE + T::ONE;
        let r = (-(two * u1.ln())).sqrt();
        let theta = two * T::PI * u2;
        r * theta.cos()
    }

    /// Approximate standard normal via Box–Muller (for tests)
    pub fn normal01(&mut self) -> f64 {
        self.normal01_real::<f64>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_sequence() {
        let mut a = XorShift64Star::new(42);
        let mut b = XorShift64Star::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn next_f64_range() {
        let mut r = XorShift64Star::new(12345);
        for _ in 0..1000 {
            let v = r.next_f64();
            assert!(v >= 0.0 && v < 1.0);
        }
    }

    #[test]
    fn next_real_range_f32() {
        let mut r = XorShift64Star::new(12345);
        for _ in 0..1000 {
            let v: f32 = r.next_real();
            assert!(v >= 0.0 && v < 1.0);
        }
    }
}
