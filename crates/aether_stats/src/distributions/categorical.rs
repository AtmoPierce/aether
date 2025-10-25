#![cfg(feature = "std")]

use num_traits::{Float, cast::cast };

/// Categorical distribution over K categories.
#[derive(Clone, Debug)]
pub struct Categorical<F: Float + > {
    pub probs: Vec<F>,
}

impl<F> Categorical<F>
where
    F: Float + ,
{
    pub fn new(probs: Vec<F>) -> Self { Self { probs } }

    pub fn pmf(&self, k: usize) -> F {
        if k < self.probs.len() { self.probs[k] } else { F::zero() }
    }

    /// MLE fit: counts normalized
    pub fn fit_mle(counts: &[usize]) -> Self {
    let total: usize = counts.iter().sum();
    let total = total.max(1);
    let probs = counts.iter().map(|&c| cast::<u32, F>(c as u32).unwrap() / cast::<u32, F>(total as u32).unwrap()).collect();
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
