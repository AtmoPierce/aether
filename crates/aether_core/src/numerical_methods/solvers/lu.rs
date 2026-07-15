use crate::math::{Matrix, Vector};
use crate::real::Real;

#[derive(Debug, Clone, Copy)]
pub struct LuDecomp<T: Real + Copy, const N: usize> {
    lu: Matrix<T, N, N>, // Combined L (unit diag) and U
    piv: [usize; N],     // Permutation such that (P*b)[i] = b[piv[i]]
    sign: i8,            // +1 or -1 depending on row-swap parity
}

impl<T: Real + Copy, const N: usize> LuDecomp<T, N> {
    #[inline]
    pub fn lu(&self) -> &Matrix<T, N, N> {
        &self.lu
    }

    #[inline]
    pub fn pivots(&self) -> &[usize; N] {
        &self.piv
    }

    #[inline]
    pub fn sign(&self) -> i8 {
        self.sign
    }

    /// Determinant of the original matrix A.
    ///
    /// Since P*A = L*U and det(L)=1:
    ///
    /// det(A) = det(P)^(-1) * det(U)
    ///
    /// Because det(P) is +/- 1, this is equivalent to sign * prod(diag(U)).
    pub fn determinant(&self) -> T {
        let mut det = if self.sign >= 0 {
            T::ONE
        } else {
            T::ZERO - T::ONE
        };

        for i in 0..N {
            det = det * self.lu[(i, i)];
        }

        det
    }

    /// Solve A*x = b using the LU factors.
    ///
    /// The stored factorization satisfies:
    ///
    ///     P*A = L*U
    ///
    /// so this solves:
    ///
    ///     L*U*x = P*b
    pub fn solve(&self, b: &Vector<T, N>) -> Vector<T, N> {
        let mut x = Vector { data: [T::ZERO; N] };

        // Apply permutation: x = P*b.
        for i in 0..N {
            x[i] = b[self.piv[i]];
        }

        // Forward solve: L*y = P*b.
        // L has unit diagonal and is stored below the diagonal of lu.
        for i in 0..N {
            let mut s = x[i];

            for j in 0..i {
                s = s - self.lu[(i, j)] * x[j];
            }

            x[i] = s;
        }

        // Backward solve: U*x = y.
        // U is stored on and above the diagonal of lu.
        for i_rev in 0..N {
            let i = N - 1 - i_rev;
            let mut s = x[i];

            for j in (i + 1)..N {
                s = s - self.lu[(i, j)] * x[j];
            }

            x[i] = s / self.lu[(i, i)];
        }

        x
    }

    /// Solve A*X = B where B is an N x K matrix of K right-hand sides.
    ///
    /// The stored factorization satisfies:
    ///
    ///     P*A = L*U
    ///
    /// so this solves:
    ///
    ///     L*U*X = P*B
    pub fn solve_multi<const K: usize>(&self, b: &Matrix<T, N, K>) -> Matrix<T, N, K> {
        let mut x = *b;

        // Apply permutation P to each RHS column.
        for col in 0..K {
            let mut tmp = [T::ZERO; N];

            for i in 0..N {
                tmp[i] = x[(self.piv[i], col)];
            }

            for i in 0..N {
                x[(i, col)] = tmp[i];
            }
        }

        // Forward solve: L*Y = P*B.
        // L has unit diagonal.
        for i in 0..N {
            for col in 0..K {
                let mut s = x[(i, col)];

                for j in 0..i {
                    s = s - self.lu[(i, j)] * x[(j, col)];
                }

                x[(i, col)] = s;
            }
        }

        // Backward solve: U*X = Y.
        for i_rev in 0..N {
            let i = N - 1 - i_rev;

            for col in 0..K {
                let mut s = x[(i, col)];

                for j in (i + 1)..N {
                    s = s - self.lu[(i, j)] * x[(j, col)];
                }

                x[(i, col)] = s / self.lu[(i, i)];
            }
        }

        x
    }
}

impl<T: Real + Copy, const N: usize> Matrix<T, N, N> {
    /// Compute LU decomposition with partial pivoting.
    ///
    /// Produces:
    ///
    ///     P*A = L*U
    ///
    /// where:
    ///
    /// - P is the row-permutation matrix,
    /// - L is lower triangular with unit diagonal,
    /// - U is upper triangular.
    ///
    /// The combined `lu` matrix stores:
    ///
    /// - L below the diagonal,
    /// - U on and above the diagonal.
    ///
    /// The permutation vector is stored such that:
    ///
    ///     (P*b)[i] = b[piv[i]]
    ///
    /// Returns `None` if the matrix is numerically singular.
    pub fn lu_decompose(&self) -> Option<LuDecomp<T, N>> {
        let mut lu = *self;

        let mut piv = [0usize; N];
        for i in 0..N {
            piv[i] = i;
        }

        let mut sign: i8 = 1;

        // Compute matrix scale for relative singularity detection.
        //
        // An absolute check like `pval <= T::EPSILON` incorrectly rejects
        // perfectly valid small-scale matrices.
        let mut scale = T::ZERO;
        for r in 0..N {
            for c in 0..N {
                let v = self[(r, c)].abs();
                if v > scale {
                    scale = v;
                }
            }
        }

        if scale <= T::ZERO {
            return None;
        }

        let singular_tol = scale * T::EPSILON * T::from_f64(N as f64);

        for k in 0..N {
            // Find pivot row for column k.
            let mut p = k;
            let mut pval = lu[(k, k)].abs();

            for r in (k + 1)..N {
                let v = lu[(r, k)].abs();

                if v > pval {
                    p = r;
                    pval = v;
                }
            }

            if pval <= singular_tol {
                return None;
            }

            // Swap rows if needed.
            if p != k {
                for c in 0..N {
                    let tmp = lu[(k, c)];
                    lu[(k, c)] = lu[(p, c)];
                    lu[(p, c)] = tmp;
                }

                piv.swap(k, p);
                sign = -sign;
            }

            // Eliminate entries below the pivot.
            let pivot = lu[(k, k)];

            for i in (k + 1)..N {
                lu[(i, k)] = lu[(i, k)] / pivot;
                let lik = lu[(i, k)];

                for j in (k + 1)..N {
                    lu[(i, j)] = lu[(i, j)] - lik * lu[(k, j)];
                }
            }
        }

        Some(LuDecomp { lu, piv, sign })
    }

    /// Convenience solve for A*x = b.
    pub fn solve(&self, b: &Vector<T, N>) -> Option<Vector<T, N>> {
        let lu = self.lu_decompose()?;
        Some(lu.solve(b))
    }

    /// Convenience solve for A*X = B where B has K right-hand sides.
    pub fn solve_multi<const K: usize>(&self, b: &Matrix<T, N, K>) -> Option<Matrix<T, N, K>> {
        let lu = self.lu_decompose()?;
        Some(lu.solve_multi(b))
    }
}