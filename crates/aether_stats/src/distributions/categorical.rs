#![cfg(feature = "std")]

use aether_core::real::Real;

/// Categorical distribution over K categories.
#[derive(Clone, Debug)]
pub struct Categorical<F: Real + > {
    pub probs: Vec<F>,
}

impl<F> Categorical<F>
where
    F: Real + ,
{
    pub fn new(probs: Vec<F>) -> Self { Self { probs } }

    pub fn pmf(&self, k: usize) -> F {
        if k < self.probs.len() { self.probs[k] } else { F::ZERO }
    }

    /// MLE fit: counts normalized
    pub fn fit_mle(counts: &[usize]) -> Self {
        let total: usize = counts.iter().sum();
        let total_f = F::from_usize(total.max(1));
        let probs: Vec<F> = counts
            .iter()
            .map(|&c| F::from_usize(c) / total_f)
            .collect();
        Self::new(probs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn categorical_fit() {
        let counts = [2usize, 3, 5];
        let cat: Categorical<f64> = Categorical::fit_mle(&counts);
        assert!((cat.probs[0] - 0.2).abs() < 1e-12);
        assert!((cat.probs[1] - 0.3).abs() < 1e-12);
        assert!((cat.probs[2] - 0.5).abs() < 1e-12);
    }
}
