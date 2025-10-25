#![cfg(feature = "std")]
use aether::math::Vector;
use crate::randomizers::XorShift64Star;
use num_traits::{Float, FloatConst, cast};

/// Univariate Gaussian distribution helpers parameterized over a floating type F.
#[derive(Clone, Copy, Debug)]
pub struct Gaussian<F: Float> {
    pub mu: F,
    pub sigma: F,
}

impl<F> Gaussian<F>
where
    F: Float + FloatConst,
{
    /// Create a new Gaussian with mean `mu` and standard deviation `sigma`.
    pub fn new(mu: F, sigma: F) -> Self { Self { mu, sigma } }

    /// Draw one sample from this Gaussian using Box–Muller and the provided RNG.
    /// The RNG yields f64 uniforms which are converted to `F`.
    pub fn sample(&self, rng: &mut XorShift64Star) -> F {
        // Box–Muller: two uniforms in (0,1]
        let mut u1 = rng.next_f64();
        if u1 <= 0.0 { u1 = core::f64::EPSILON; }
        let u2 = rng.next_f64();
        // do math in f64 (RNG gives f64) then cast final z to F
        let r = (-2.0_f64 * u1.ln()).sqrt();
        let theta: f64 = 2.0_f64 * core::f64::consts::PI * u2;
        let z0 = r * theta.cos();
        let z: F = num_traits::cast(z0).expect("conversion from f64 to float failed");
        z * self.sigma + self.mu
    }

    /// Probability density function evaluated at x
    pub fn pdf(&self, x: F) -> F {
        let two = F::one() + F::one();
        let half = F::one() / two;
        let tiny = F::zero();
        let pi = num_traits::FloatConst::PI();
        let var = (self.sigma * self.sigma).max(tiny);
        let coeff = F::one() / ((two * pi * var).sqrt());
        let diff = x - self.mu;
        coeff * (-(half * diff * diff / var)).exp()
    }

    /// Fit MLE for a set of scalar samples (samples are of type F)
    pub fn fit_mle(samples: &[F]) -> Self {
        let n_f: F = cast(samples.len() as u32).unwrap();
        if samples.is_empty() { return Self::new(F::zero(), F::one()); }
        let mut sum = F::zero();
        for &s in samples { sum = sum + s; }
        let mu = sum / n_f;
        let mut s2 = F::zero();
        for &s in samples { let d = s - mu; s2 = s2 + d * d; }
        let sigma_squared = s2 / n_f; // MLE uses 1/n
        let sigma = sigma_squared.sqrt().max(F::epsilon());
        Self { mu, sigma }
    }
}

// Keep a convenient alias for f64 users
pub type GaussianF64 = Gaussian<f64>;

// Backward-compatible re-export: allow `Gaussian::new` usage with f64 by using type inference
impl GaussianF64 {
    /// Convenience constructor specialized for f64
    pub fn new_f64(mu: f64, sigma: f64) -> Self { Self::new(mu, sigma) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gaussian_pdf_center() {
        let g = Gaussian::new(0.0, 1.0);
        let p = g.pdf(0.0);
        assert!( (p - (1.0 / (2.0 * std::f64::consts::PI).sqrt())).abs() < 1e-12 );
    }

    #[test]
    fn fit_simple() {
        let samples = [1.0, 2.0, 3.0, 2.0, 1.0];
        let g = Gaussian::fit_mle(&samples);
        assert!((g.mu - 1.8).abs() < 1e-12);
    }
}

#[cfg(test)]
mod wn_tests {
    use super::*;

    #[test]
    fn white_noise_deterministic() {
        let mu = 0.5;
        let sigma = 2.0;
        let n = 100;
        let seed = 42u64;
        let v1 = (0..n).map(|_| {
            let mut rng = XorShift64Star::new(seed);
            let g = Gaussian::new(mu, sigma);
            // produce only first sample repeatedly to check determinism per-seed
            g.sample(&mut rng)
        }).collect::<Vec<_>>();
        let v2 = (0..n).map(|_| {
            let mut rng = XorShift64Star::new(seed);
            let g = Gaussian::new(mu, sigma);
            g.sample(&mut rng)
        }).collect::<Vec<_>>();
        assert_eq!(v1, v2);
    }
}
