// ---- tiny RNG reused from your no_std version ----
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
    pub fn next_f64(&mut self) -> f64 {
        const DEN: f64 = (1u64 << 53) as f64;
        ((self.next_u64() >> 11) as f64) / DEN
    }

    /// Approximate standard normal via Box–Muller (for tests)
    pub fn normal01(&mut self) -> f64 {
        let mut u1 = self.next_f64();
        if u1 <= 0.0 {
            u1 = std::f64::EPSILON;
        }
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        r * theta.cos()
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
}
