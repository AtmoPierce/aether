use aether_core::real::Real;

use super::WhiteNoise;

#[derive(Debug, Clone, Copy)]
pub struct GaussMarkov1<F: Real> {
    state: F,
    tau: F,
    steady_state_sigma: F,
    unit_noise: WhiteNoise<F>,
}

impl<F: Real> GaussMarkov1<F> {
    pub fn new(initial_state: F, tau: F, steady_state_sigma: F, seed: u64) -> Self {
        Self {
            state: initial_state,
            tau,
            steady_state_sigma,
            unit_noise: WhiteNoise::new(F::ZERO, F::ONE, seed),
        }
    }

    pub fn state(&self) -> F {
        self.state
    }

    pub fn tau(&self) -> F {
        self.tau
    }

    pub fn steady_state_sigma(&self) -> F {
        self.steady_state_sigma
    }

    pub fn reset(&mut self, state: F) {
        self.state = state;
    }

    pub fn advance(&mut self, dt: F) -> F {
        if dt <= F::ZERO {
            return self.state;
        }

        if self.tau <= F::EPSILON {
            self.state = self.steady_state_sigma * self.unit_noise.sample();
            return self.state;
        }

        let phi = (-(dt / self.tau)).exp();
        let sigma_drive = (F::ONE - phi * phi).max(F::ZERO).sqrt() * self.steady_state_sigma;
        self.state = phi * self.state + sigma_drive * self.unit_noise.sample();
        self.state
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_for_same_seed() {
        let mut a = GaussMarkov1::new(0.0_f64, 10.0_f64, 0.5_f64, 42);
        let mut b = GaussMarkov1::new(0.0_f64, 10.0_f64, 0.5_f64, 42);

        for _ in 0..16 {
            assert_eq!(a.advance(0.1), b.advance(0.1));
        }
    }

    #[test]
    fn zero_dt_keeps_state() {
        let mut process = GaussMarkov1::new(1.5_f64, 10.0_f64, 0.3_f64, 7);
        assert_eq!(process.advance(0.0), 1.5);
    }
}