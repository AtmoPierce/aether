use aether_core::real::Real;

use super::WhiteNoise;

#[derive(Debug, Clone, Copy)]
pub struct RandomWalk<F: Real> {
    state: F,
    driving_noise: WhiteNoise<F>,
}

impl<F: Real> RandomWalk<F> {
    pub fn new(initial_state: F, sigma_per_sqrt_second: F, seed: u64) -> Self {
        Self {
            state: initial_state,
            driving_noise: WhiteNoise::new(F::ZERO, sigma_per_sqrt_second, seed),
        }
    }

    pub fn state(&self) -> F {
        self.state
    }

    pub fn reset(&mut self, state: F) {
        self.state = state;
    }

    pub fn advance(&mut self, dt: F) -> F {
        if dt <= F::ZERO {
            return self.state;
        }

        self.state = self.state + self.driving_noise.sample() * dt.sqrt();
        self.state
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_for_same_seed() {
        let mut a = RandomWalk::new(0.0_f64, 1.0_f64, 42);
        let mut b = RandomWalk::new(0.0_f64, 1.0_f64, 42);

        for _ in 0..16 {
            assert_eq!(a.advance(0.1), b.advance(0.1));
        }
    }

    #[test]
    fn zero_dt_keeps_state() {
        let mut walk = RandomWalk::new(1.5_f64, 0.3_f64, 7);
        assert_eq!(walk.advance(0.0), 1.5);
    }
}