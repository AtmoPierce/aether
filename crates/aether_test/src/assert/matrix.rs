use aether_core::math::Matrix;
use approx::abs_diff_eq;
use thiserror::Error;

/* ----------------------------- Matrix Asserts ----------------------------- */

#[derive(Debug, Error)]
pub enum MatrixAssertError {
    #[error("Matrix shape mismatch: left={left_rows}x{left_cols}, right={right_rows}x{right_cols}")]
    ShapeMismatch {
        left_rows: usize,
        left_cols: usize,
        right_rows: usize,
        right_cols: usize,
    },

    #[error(
        "Value mismatch at (row={row}, col={col}): left={left}, right={right}, |Δ|={diff}"
    )]
    ValueMismatch {
        row: usize,
        col: usize,
        left: f64,
        right: f64,
        diff: f64,
    },
}

/// Panic-on-failure approximate equality for matrices (elementwise).
///
/// # Example
/// ```
/// use aether_core::math::Matrix;
/// use aether_test::assert::matrix::assert_matrix_approx_eq;
///
/// let a = Matrix::new([[1.0, 2.0], [3.0, 4.0]]);
/// let b = Matrix::new([[1.0, 2.0001], [2.9999, 4.0]]);
/// assert_matrix_approx_eq(&a, &b, 1e-3);
/// ```
pub fn assert_matrix_approx_eq<const M: usize, const N: usize>(
    a: &Matrix<f64, M, N>,
    b: &Matrix<f64, M, N>,
    tol: f64,
) {
    if let Err(err) = matrix_approx_eq(a, b, tol) {
        panic!("{:?}", err);
    }
}

/// Non-panicking version that returns a `Result`.
pub fn matrix_approx_eq<const M: usize, const N: usize>(
    a: &Matrix<f64, M, N>,
    b: &Matrix<f64, M, N>,
    tol: f64,
) -> Result<(), MatrixAssertError> {
    // Shapes are compile-time equal for M,N, but keep a defensive check for
    // situations where users may adapt this to dynamic shapes later.
    let (lr, lc) = (M, N);
    let (rr, rc) = (M, N);
    if (lr, lc) != (rr, rc) {
        return Err(MatrixAssertError::ShapeMismatch {
            left_rows: lr,
            left_cols: lc,
            right_rows: rr,
            right_cols: rc,
        });
    }

    for r in 0..M {
        for c in 0..N {
            let x = a.data[r][c];
            let y = b.data[r][c];
            if !abs_diff_eq!(x, y, epsilon = tol) {
                return Err(MatrixAssertError::ValueMismatch {
                    row: r,
                    col: c,
                    left: x,
                    right: y,
                    diff: (x - y).abs(),
                });
            }
        }
    }
    Ok(())
}

/// Maximum absolute elementwise difference (Linf elementwise norm).
pub fn matrix_max_abs_diff<const M: usize, const N: usize>(
    a: &Matrix<f64, M, N>,
    b: &Matrix<f64, M, N>,
) -> f64 {
    let mut maxd = 0.0_f64;
    for r in 0..M {
        for c in 0..N {
            let d = (a.data[r][c] - b.data[r][c]).abs();
            if d > maxd {
                maxd = d;
            }
        }
    }
    maxd
}

/// Frobenius norm of the matrix (sqrt of sum of squares).
pub fn matrix_frobenius_norm<const M: usize, const N: usize>(a: &Matrix<f64, M, N>) -> f64 {
    let mut sum = 0.0_f64;
    for r in 0..M {
        for c in 0..N {
            let v = a.data[r][c];
            sum += v * v;
        }
    }
    sum.sqrt()
}

/// Frobenius norm of the difference (||A - B||F).
pub fn matrix_diff_frobenius<const M: usize, const N: usize>(
    a: &Matrix<f64, M, N>,
    b: &Matrix<f64, M, N>,
) -> f64 {
    let mut sum = 0.0_f64;
    for r in 0..M {
        for c in 0..N {
            let d = a.data[r][c] - b.data[r][c];
            sum += d * d;
        }
    }
    sum.sqrt()
}

/// Assert that all elements are finite (no NaN/Inf).
pub fn assert_matrix_finite<const M: usize, const N: usize>(a: &Matrix<f64, M, N>) {
    for r in 0..M {
        for c in 0..N {
            let v = a.data[r][c];
            assert!(
                v.is_finite(),
                "Matrix element at (row={}, col={}) is not finite: {}",
                r,
                c,
                v
            );
        }
    }
}

/// Assert symmetry for square matrices within a tolerance:
/// checks |A_ij - A_ji| <= tol for all i,j.
///
/// # Example
/// ```
/// use aether_core::math::Matrix;
/// use aether_test::assert::matrix::assert_matrix_symmetric;
/// let a = Matrix::new([[1.0, 2.0], [2.0, 3.0]]);
/// assert_matrix_symmetric(&a, 1e-12);
/// ```
pub fn assert_matrix_symmetric<const N: usize>(a: &Matrix<f64, N, N>, tol: f64) {
    for i in 0..N {
        for j in (i + 1)..N {
            let d = (a.data[i][j] - a.data[j][i]).abs();
            assert!(
                d <= tol,
                "Matrix not symmetric at ({},{}) vs ({},{}): |{} - {}| = {} > {}",
                i,
                j,
                j,
                i,
                a.data[i][j],
                a.data[j][i],
                d,
                tol
            );
        }
    }
}


