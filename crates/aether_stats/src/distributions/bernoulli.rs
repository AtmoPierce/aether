use aether_core::real::Real;

/// Bernoulli distribution: P(X=1)=p, P(X=0)=1-p
#[derive(Clone, Debug)]
pub struct Bernoulli<F: Real> {
    pub p: F,
}

impl<F> Bernoulli<F>
where
    F: Real,
{
    pub fn new(p: F) -> Self {
        Self { p }
    }

    pub fn pmf(&self, x: u8) -> F {
        match x {
            0 => F::ONE - self.p,
            1 => self.p,
            _ => F::ZERO,
        }
    }

    /// MLE fit for bernoulli is just sample mean (assumes x in {0,1})
    pub fn fit_mle(samples: &[u8]) -> Self {
        let mut sum = 0usize;
        for &s in samples {
            sum += s as usize;
        }
        let total = samples.len().max(1);
        let p = F::from_usize(sum) / F::from_usize(total);
        Self::new(p)
    }
}

pub type BernoulliF64 = Bernoulli<f64>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bernoulli_fit() {
        let samples = [1u8, 0, 1, 1, 0];
        let b: Bernoulli<f64> = Bernoulli::fit_mle(&samples);
        assert!((b.p - 0.6).abs() < 1e-12);
    }
}
