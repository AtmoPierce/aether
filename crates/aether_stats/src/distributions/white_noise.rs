use aether_core::real::Real;
use aether_rand::randomizers::XorShift64Star;

use super::Gaussian;

#[derive(Debug, Clone, Copy)]
pub struct WhiteNoise<F: Real> {
    distribution: Gaussian<F>,
    rng: XorShift64Star,
}

impl<F: Real> WhiteNoise<F> {
    pub fn new(mean: F, sigma: F, seed: u64) -> Self {
        Self {
            distribution: Gaussian::new(mean, sigma),
            rng: XorShift64Star::new(seed),
        }
    }

    pub fn from_gaussian(distribution: Gaussian<F>, seed: u64) -> Self {
        Self {
            distribution,
            rng: XorShift64Star::new(seed),
        }
    }

    pub fn sample(&mut self) -> F {
        self.distribution.sample(&mut self.rng)
    }

    pub fn mean(&self) -> F {
        self.distribution.mu
    }

    pub fn sigma(&self) -> F {
        self.distribution.sigma
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_for_same_seed() {
        let mut a = WhiteNoise::new(0.0_f64, 1.0_f64, 42);
        let mut b = WhiteNoise::new(0.0_f64, 1.0_f64, 42);

        for _ in 0..32 {
            assert_eq!(a.sample(), b.sample());
        }
    }

    #[test]
    fn exposes_distribution_parameters() {
        let noise = WhiteNoise::new(0.5_f64, 2.0_f64, 7);
        assert_eq!(noise.mean(), 0.5);
        assert_eq!(noise.sigma(), 2.0);
    }
}