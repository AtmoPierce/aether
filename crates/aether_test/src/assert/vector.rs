use aether_core::math::Vector;
use approx::{abs_diff_eq, AbsDiffEq};
use thiserror::Error;

/// Errors that can occur when comparing vectors
#[derive(Debug, Error)]
pub enum VectorAssertError {
    #[error("Vector length mismatch: left={left}, right={right}")]
    LengthMismatch { left: usize, right: usize },

    #[error("Value mismatch at index {index}: left={left}, right={right}, |Δ|={diff}")]
    ValueMismatch {
        index: usize,
        left: f64,
        right: f64,
        diff: f64,
    },
}

/// Assert that two vectors are approximately equal within a tolerance.
///
/// Panics on error unless compiled with `--no-default-features` (where it returns a Result).
///
/// # Example
/// ```
/// use aether_core::math::Vector;
/// use aether_test::assert::vector::assert_vector_approx_eq;
///
/// let v1 = Vector::from([1.0, 2.0, 3.0]);
/// let v2 = Vector::from([1.0, 2.0001, 2.9999]);
/// assert_vector_approx_eq(&v1, &v2, 1e-3);
/// ```
pub fn assert_vector_approx_eq<const N: usize>(
    a: &Vector<f64, N>,
    b: &Vector<f64, N>,
    tol: f64,
) {
    if let Err(err) = vector_approx_eq(a, b, tol) {
        panic!("{:?}", err);
    }
}

/// Non-panicking comparison helper for programmatic use.
///
/// Returns `Ok(())` if all elements are within `tol`, or `Err(VectorAssertError)` otherwise.
pub fn vector_approx_eq<const N: usize>(
    a: &Vector<f64, N>,
    b: &Vector<f64, N>,
    tol: f64,
) -> Result<(), VectorAssertError> {
    if a.data.len() != b.data.len() {
        return Err(VectorAssertError::LengthMismatch {
            left: a.data.len(),
            right: b.data.len(),
        });
    }

    for (i, (&x, &y)) in a.data.iter().zip(b.data.iter()).enumerate() {
        if !abs_diff_eq!(x, y, epsilon = tol) {
            return Err(VectorAssertError::ValueMismatch {
                index: i,
                left: x,
                right: y,
                diff: (x - y).abs(),
            });
        }
    }
    Ok(())
}

/// Compute the maximum absolute difference between two vectors.
pub fn vector_max_abs_diff<const N: usize>(a: &Vector<f64, N>, b: &Vector<f64, N>) -> f64 {
    a.data
        .iter()
        .zip(b.data.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

/// Assert that all elements of a vector are finite (no NaN or Inf).
pub fn assert_vector_finite<const N: usize>(v: &Vector<f64, N>) {
    for (i, &x) in v.data.iter().enumerate() {
        assert!(
            x.is_finite(),
            "Vector element at index {} is not finite: {}",
            i,
            x
        );
    }
}
