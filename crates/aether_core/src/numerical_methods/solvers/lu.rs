use crate::real::Real;
use crate::math::{Matrix, Vector};

#[derive(Debug, Clone, Copy)]
pub struct LuDecomp<T: Real + Copy, const N: usize> {
    pub lu: Matrix<T, N, N>,       // Combined L (unit diag) and U
    pub piv: [usize; N],           // Row permutations applied to A
    pub sign: i8,                  // +1 or -1 depending on parity
}

impl<T: Real + Copy, const N: usize> Matrix<T, N, N> {
    /// Compute LU with partial pivoting. A = P*L*U
    /// Returns None if numerically singular.
    pub fn lu_decompose(&self) -> Option<LuDecomp<T, N>> {
        let mut lu = *self;
        let mut piv = [0usize; N];
        for i in 0..N { piv[i] = i; }
        let mut sign: i8 = 1;

        for k in 0..N {
            // Find pivot row
            let mut p = k;
            let mut pval = lu[(k, k)].abs();
            for r in (k + 1)..N {
                let v = lu[(r, k)].abs();
                if v > pval { p = r; pval = v; }
            }
            if pval <= T::EPSILON { return None; } // singular

            // Swap rows if needed
            if p != k {
                for c in 0..N {
                    let tmp = lu[(k, c)];
                    lu[(k, c)] = lu[(p, c)];
                    lu[(p, c)] = tmp;
                }
                piv.swap(k, p);
                sign = -sign;
            }

            // Factorization step
            let pivot = lu[(k, k)];
            for i in (k + 1)..N {
                lu[(i, k)] = lu[(i, k)] / pivot;
                let lik = lu[(i, k)];
                // rank-1 update of the trailing block
                for j in (k + 1)..N {
                    lu[(i, j)] = lu[(i, j)] - lik * lu[(k, j)];
                }
            }
        }

        Some(LuDecomp { lu, piv, sign })
    }
}

impl<T: Real + Copy, const N: usize> LuDecomp<T, N> {
    /// Solve Ax=b using the LU factors (single RHS).
    pub fn solve(&self, b: &Vector<T, N>) -> Vector<T, N> {
        // Apply permutation to b -> Pb
        let mut x = Vector { data: [T::ZERO; N] };
        for i in 0..N {
            x[i] = b[self.piv[i]];
        }

        // Forward solve: Ly = Pb  (L is unit diagonal, stored below diag in lu)
        for i in 0..N {
            let mut s = x[i];
            for j in 0..i {
                s = s - self.lu[(i, j)] * x[j];
            }
            x[i] = s; // since diag(L)=1
        }

        // Backward solve: Ux=y  (U is upper triangular, stored on/above diag in lu)
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

    /// Solve AX=B where B is (N×K) matrix of K RHS.
    pub fn solve_multi<const K: usize>(&self, b: &Matrix<T, N, K>) -> Matrix<T, N, K> {
        let mut x = *b; // will hold the permuted RHS and eventually the solution

        // Apply permutation P to each RHS column
        for col in 0..K {
            let mut tmp = [T::ZERO; N];
            for i in 0..N {
                tmp[i] = x[self.piv[i]][col];
            }
            for i in 0..N {
                x[i][col] = tmp[i];
            }
        }

        // Forward solve: L Y = P B
        for i in 0..N {
            for col in 0..K {
                let mut s = x[i][col];
                for j in 0..i {
                    s = s - self.lu[(i, j)] * x[j][col];
                }
                x[i][col] = s; // diag(L)=1
            }
        }

        // Backward solve: U X = Y
        for i_rev in 0..N {
            let i = N - 1 - i_rev;
            for col in 0..K {
                let mut s = x[i][col];
                for j in (i + 1)..N {
                    s = s - self.lu[(i, j)] * x[j][col];
                }
                x[i][col] = s / self.lu[(i, i)];
            }
        }

        x
    }
}

/// Convenience: A.solve(b) -> Option<x>
impl<T: Real + Copy, const N: usize> Matrix<T, N, N> {
    pub fn solve(&self, b: &Vector<T, N>) -> Option<Vector<T, N>> {
        let lu = self.lu_decompose()?;
        Some(lu.solve(b))
    }

    pub fn solve_multi<const K: usize>(&self, b: &Matrix<T, N, K>) -> Option<Matrix<T, N, K>> {
        let lu = self.lu_decompose()?;
        Some(lu.solve_multi(b))
    }
}
