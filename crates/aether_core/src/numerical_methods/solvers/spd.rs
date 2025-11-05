use num_traits::Float;
use crate::math::{Matrix, Vector};

#[derive(Debug, Clone, Copy)]
pub struct CholLower<T: Float + Copy, const N: usize> {
    pub l: Matrix<T, N, N>, // lower-triangular with positive diagonal
}

impl<T: Float + Copy, const N: usize> Matrix<T, N, N> {
    /// Compute lower Cholesky: A = L * L^t. Returns None if not SPD.
    pub fn cholesky_lower(&self) -> Option<CholLower<T, N>> {
        let mut l = Matrix::<T, N, N>::zeros();

        for i in 0..N {
            for j in 0..=i {
                let mut s = self[(i, j)];
                for k in 0..j {
                    s = s - l[(i, k)] * l[(j, k)];
                }
                if i == j {
                    if s <= T::zero() { return None; }
                    l[(i, j)] = s.sqrt();
                } else {
                    l[(i, j)] = s / l[(j, j)];
                }
            }
            // zero the strict upper part
            for j in (i + 1)..N {
                l[(i, j)] = T::zero();
            }
        }

        Some(CholLower { l })
    }
}

impl<T: Float + Copy, const N: usize> CholLower<T, N> {
    /// Solve A x = b using Cholesky factors (A = L L^t).
    pub fn solve(&self, b: &Vector<T, N>) -> Vector<T, N> {
        let l = &self.l;

        // Forward solve: L y = b
        let mut y = *b;
        for i in 0..N {
            let mut s = y[i];
            for k in 0..i {
                s = s - l[(i, k)] * y[k];
            }
            y[i] = s / l[(i, i)];
        }

        // Backward solve: L^t x = y
        let mut x = y;
        for i_rev in 0..N {
            let i = N - 1 - i_rev;
            let mut s = x[i];
            for k in (i + 1)..N {
                s = s - l[(k, i)] * x[k];
            }
            x[i] = s / l[(i, i)];
        }

        x
    }
}

/// Convenience: A.cholesky_solve(b)
impl<T: Float + Copy, const N: usize> Matrix<T, N, N> {
    pub fn cholesky_solve(&self, b: &Vector<T, N>) -> Option<Vector<T, N>> {
        let chol = self.cholesky_lower()?;
        Some(chol.solve(b))
    }
}
