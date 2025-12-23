use super::vector::Vector;
use crate::real::{Real};
use core::ops::{Add, Sub, Mul, Div, Neg, AddAssign};
use core::ops::{Index, IndexMut};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix<T, const M: usize, const N: usize> {
    pub data: [[T; N]; M],
}

impl<T: Copy, const M: usize, const N: usize> Matrix<T, M, N> {
    #[inline]
    pub fn new(data: [[T; N]; M]) -> Self {
        Self { data }
    }
}

/* -------------------- Default -------------------- */

impl<T: Real, const M: usize, const N: usize> Default for Matrix<T, M, N> {
    fn default() -> Self {
        Self {
            data: [[T::ZERO; N]; M],
        }
    }
}

/* -------------------- Basic ops -------------------- */

// Matrix addition
impl<T: Real, const M: usize, const N: usize> Add for Matrix<T, M, N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut result = self;
        for r in 0..M {
            for c in 0..N {
                result.data[r][c] = self.data[r][c] + rhs.data[r][c];
            }
        }
        result
    }
}

// Matrix subtraction
impl<T: Real, const M: usize, const N: usize> Sub for Matrix<T, M, N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let mut result = self;
        for r in 0..M {
            for c in 0..N {
                result.data[r][c] = self.data[r][c] - rhs.data[r][c];
            }
        }
        result
    }
}

// Scalar division
impl<T: Real, const M: usize, const N: usize> Div<T> for Matrix<T, M, N> {
    type Output = Self;

    fn div(self, rhs: T) -> Self {
        let mut result = self;
        for r in 0..M {
            for c in 0..N {
                result.data[r][c] = self.data[r][c] / rhs;
            }
        }
        result
    }
}

// Negation
impl<T: Real, const M: usize, const N: usize> Neg for Matrix<T, M, N> {
    type Output = Self;

    fn neg(self) -> Self {
        let mut result = self;
        for r in 0..M {
            for c in 0..N {
                result.data[r][c] = -self.data[r][c];
            }
        }
        result
    }
}

// Scalar multiplication: Matrix<T> * T
impl<T: Real, const M: usize, const N: usize> Mul<T> for Matrix<T, M, N> {
    type Output = Matrix<T, M, N>;

    fn mul(self, rhs: T) -> Self::Output {
        let mut result = self;
        for r in 0..M {
            for c in 0..N {
                result.data[r][c] = self.data[r][c] * rhs;
            }
        }
        result
    }
}

/* -------------------- Matrix × Matrix / Vector -------------------- */

// Matrix multiplication: (M×N) * (N×P) -> (M×P)
impl<T: Real, const M: usize, const N: usize, const P: usize> Mul<Matrix<T, N, P>>
    for Matrix<T, M, N>
{
    type Output = Matrix<T, M, P>;

    fn mul(self, rhs: Matrix<T, N, P>) -> Matrix<T, M, P> {
        let mut result = Matrix {
            data: [[T::ZERO; P]; M],
        };

        for i in 0..M {
            for j in 0..P {
                let mut sum = T::ZERO;
                for k in 0..N {
                    sum = sum + self.data[i][k] * rhs.data[k][j];
                }
                result.data[i][j] = sum;
            }
        }

        result
    }
}

// Matrix × Vector: (M×N) * (N) -> (M)
impl<T: Real, const M: usize, const N: usize> Mul<Vector<T, N>> for Matrix<T, M, N> {
    type Output = Vector<T, M>;

    fn mul(self, rhs: Vector<T, N>) -> Vector<T, M> {
        let mut result = Vector {
            data: [T::ZERO; M],
        };

        for i in 0..M {
            let mut sum = T::ZERO;
            for j in 0..N {
                sum = sum + self.data[i][j] * rhs.data[j];
            }
            result.data[i] = sum;
        }

        result
    }
}

/* -------------------- Constructors -------------------- */

impl<T: Real, const M: usize, const N: usize> Matrix<T, M, N> {
    #[inline]
    pub fn zeros() -> Self {
        Self {
            data: [[T::ZERO; N]; M],
        }
    }

    #[inline]
    pub fn ones() -> Self {
        Self {
            data: [[T::ONE; N]; M],
        }
    }
}

impl<T: Real, const N: usize> Matrix<T, N, N> {
    /// Identity matrix (NxN)
    #[inline]
    pub fn identity() -> Self {
        let mut data = [[T::ZERO; N]; N];
        for i in 0..N {
            data[i][i] = T::ONE;
        }
        Self { data }
    }

    /// Diagonal matrix from a slice/array of length N
    #[inline]
    pub fn diag(diag: &[T; N]) -> Self {
        let mut data = [[T::ZERO; N]; N];
        for i in 0..N {
            data[i][i] = diag[i];
        }
        Self { data }
    }

    #[inline]
    pub fn diag_from_vector(diag: &Vector<T, N>) -> Self {
        let mut data = [[T::ZERO; N]; N];
        for i in 0..N {
            data[i][i] = diag[i];
        }
        Self { data }
    }
}

/* -------------------- Outer Product -------------------- */

impl<T: Real, const M: usize, const N: usize> Matrix<T, M, N> {
    /// Outer (dyadic) product: a (Mx1) * b^T (1xN) -> (MxN)
    #[inline]
    pub fn outer(a: &Vector<T, M>, b: &Vector<T, N>) -> Self {
        let mut m = Self::zeros();
        for i in 0..M {
            for j in 0..N {
                m.data[i][j] = a[i] * b[j];
            }
        }
        m
    }
}

/* -------------------- Trace -------------------- */

impl<T: Real + AddAssign, const N: usize> Matrix<T, N, N> {
    #[inline]
    pub fn trace(&self) -> T {
        let mut sum = T::ZERO;
        for i in 0..N {
            sum += self.data[i][i];
        }
        sum
    }
}

/* -------------------- Determinant (generic N×N) -------------------- */

impl<T: Real, const N: usize> Matrix<T, N, N> {
    /// Compute the determinant using an LU-style elimination with partial pivoting.
    /// Returns zero if the matrix is (near-)singular up to a small tolerance.
    pub fn determinant(&self) -> T {
        let mut a = *self;
        let mut det = T::ONE;

        // crude tolerance;
        let eps = T::from_f64(1e-12);

        for k in 0..N {
            // pivot selection
            let mut piv = k;
            let mut piv_val = a[(k, k)].abs();
            for r in (k + 1)..N {
                let v = a[(r, k)].abs();
                if v > piv_val {
                    piv = r;
                    piv_val = v;
                }
            }

            if piv_val <= eps {
                return T::ZERO;
            }

            // swap rows
            if piv != k {
                for c in 0..N {
                    let ta = a[(k, c)];
                    a[(k, c)] = a[(piv, c)];
                    a[(piv, c)] = ta;
                }
                det = -det;
            }

            let pivot = a[(k, k)];
            det = det * pivot;

            // eliminate below
            for r in (k + 1)..N {
                let f = a[(r, k)] / pivot;
                for c in (k + 1)..N {
                    a[(r, c)] = a[(r, c)] - f * a[(k, c)];
                }
            }
        }

        det
    }
}

/* -------------------- 3×3 inverse -------------------- */

impl<T: Real> Matrix<T, 3, 3> {
    pub fn inverse(&self) -> Option<Self> {
        let m = &self.data;

        let det = self.determinant();
        let eps = T::from_f64(1e-12);
        if det.abs() <= eps {
            return None;
        }

        let inv_det = T::ONE / det;

        let mut inv = [[T::ZERO; 3]; 3];

        inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
        inv[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * inv_det;
        inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;

        inv[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * inv_det;
        inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
        inv[1][2] = -(m[0][0] * m[1][2] - m[1][0] * m[0][2]) * inv_det;

        inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
        inv[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * inv_det;
        inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;

        Some(Matrix::new(inv))
    }
}

/* -------------------- 4×4 helpers -------------------- */

impl<T: Real> Matrix<T, 4, 4> {
    pub fn as_flat_array(&self) -> [T; 16] {
        let mut flat = [self.data[0][0]; 16];
        for i in 0..4 {
            for j in 0..4 {
                flat[i * 4 + j] = self.data[i][j];
            }
        }
        flat
    }
}

/* -------------------- Transpose -------------------- */

impl<T: Real, const M: usize, const N: usize> Matrix<T, M, N> {
    pub fn transpose(&self) -> Matrix<T, N, M> {
        let mut out = Matrix::<T, N, M>::zeros();
        for r in 0..M {
            for c in 0..N {
                out[(c, r)] = self[(r, c)];
            }
        }
        out
    }

    pub fn max_abs(&self) -> T {
        let mut m = T::ZERO;
        for r in 0..M {
            for c in 0..N {
                let v = self[(r, c)].abs();
                if v > m {
                    m = v;
                }
            }
        }
        m
    }
}

/* -------------------- Gauss–Jordan inverse (generic N×N) -------------------- */

impl<T: Real, const N: usize> Matrix<T, N, N> {
    /// Invert via Gauss–Jordan with partial pivoting.
    /// Returns None if (near-)singular.
    pub fn inverse_gauss_jordan(&self) -> Option<Self> {
        let mut a = *self;
        let mut i = Matrix::<T, N, N>::identity();

        let eps = T::from_f64(1e-12);

        for k in 0..N {
            // pivot row
            let mut piv = k;
            let mut piv_val = a[(k, k)].abs();
            for r in (k + 1)..N {
                let v = a[(r, k)].abs();
                if v > piv_val {
                    piv = r;
                    piv_val = v;
                }
            }
            if piv_val <= eps {
                return None;
            }

            // swap rows in A and I
            if piv != k {
                for c in 0..N {
                    let ta = a[(k, c)];
                    a[(k, c)] = a[(piv, c)];
                    a[(piv, c)] = ta;

                    let ti = i[(k, c)];
                    i[(k, c)] = i[(piv, c)];
                    i[(piv, c)] = ti;
                }
            }

            // normalize pivot row
            let invd = T::ONE / a[(k, k)];
            for c in 0..N {
                a[(k, c)] = a[(k, c)] * invd;
                i[(k, c)] = i[(k, c)] * invd;
            }

            // eliminate others
            for r in 0..N {
                if r == k {
                    continue;
                }
                let f = a[(r, k)];
                if f != T::ZERO {
                    for c in 0..N {
                        a[(r, c)] = a[(r, c)] - f * a[(k, c)];
                        i[(r, c)] = i[(r, c)] - f * i[(k, c)];
                    }
                }
            }
        }
        Some(i)
    }
}

/* -------------------- Indexing -------------------- */

// [m][n]
impl<T, const M: usize, const N: usize> Index<usize> for Matrix<T, M, N> {
    type Output = [T; N];

    #[inline]
    fn index(&self, row: usize) -> &Self::Output {
        &self.data[row]
    }
}

impl<T, const M: usize, const N: usize> IndexMut<usize> for Matrix<T, M, N> {
    #[inline]
    fn index_mut(&mut self, row: usize) -> &mut Self::Output {
        &mut self.data[row]
    }
}

// [(m, n)]
impl<T, const M: usize, const N: usize> Index<(usize, usize)> for Matrix<T, M, N> {
    type Output = T;

    #[inline]
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (r, c) = index;
        &self.data[r][c]
    }
}

impl<T, const M: usize, const N: usize> IndexMut<(usize, usize)> for Matrix<T, M, N> {
    #[inline]
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (r, c) = index;
        &mut self.data[r][c]
    }
}

/* -------------------- Casting ----------------------------- */
impl<const M: usize, const N: usize> From<Matrix<f32, M, N>>
    for Matrix<f64, M, N>
{
    #[inline]
    fn from(src: Matrix<f32, M, N>) -> Self {
        let mut out = Matrix::<f64, M, N>::zeros();
        for r in 0..M {
            for c in 0..N {
                out[(r, c)] = src[(r, c)] as f64;
            }
        }
        out
    }
}

impl<const M: usize, const N: usize> From<Matrix<f64, M, N>>
    for Matrix<f32, M, N>
{
    #[inline]
    fn from(src: Matrix<f64, M, N>) -> Self {
        let mut out = Matrix::<f32, M, N>::zeros();
        for r in 0..M {
            for c in 0..N {
                out[(r, c)] = src[(r, c)] as f32;
            }
        }
        out
    }
}

impl<T: Real, const M: usize, const N: usize> Matrix<T, M, N> {
    #[inline]
    pub fn cast<U: Real>(self) -> Matrix<U, M, N> {
        let mut out = Matrix::<U, M, N>::zeros();
        for r in 0..M {
            for c in 0..N {
                out[(r, c)] = U::from_f64(self[(r, c)].to_f64());
            }
        }
        out
    }
}

/* -------------------- std-only Display -------------------- */

#[cfg(feature = "std")]
impl<T, const M: usize, const N: usize> std::fmt::Display for Matrix<T, M, N>
where
    T: Real + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[")?;
        for (i, row) in self.data.iter().enumerate() {
            for (j, val) in row.iter().enumerate() {
                write!(f, "{}", val)?;
                if j < N - 1 {
                    write!(f, " ")?;
                }
            }
            if i < M - 1 {
                write!(f, "; ")?;
            }
        }
        write!(f, "]")
    }
}
