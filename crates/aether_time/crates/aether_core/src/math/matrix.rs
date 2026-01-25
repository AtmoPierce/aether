use super::vector::Vector;
use core::ops::{Add, Div, Mul, Neg, Sub};
use num_traits::Float;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix<T, const M: usize, const N: usize> {
    pub data: [[T; N]; M],
}

impl<T: Copy, const M: usize, const N: usize> Matrix<T, M, N> {
    pub fn new(data: [[T; N]; M]) -> Self {
        Self { data }
    }
}
impl<T: Default + Copy + num_traits::Zero, const M: usize, const N: usize> Default
    for Matrix<T, M, N>
{
    fn default() -> Self {
        Self {
            data: [[T::default(); N]; M],
        }
    }
}

// Matrix addition
impl<T, const M: usize, const N: usize> Add for Matrix<T, M, N>
where
    T: Float + Add<Output = T> + Copy,
{
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
impl<T, const M: usize, const N: usize> Sub for Matrix<T, M, N>
where
    T: Float + Sub<Output = T> + Copy,
{
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
impl<T, const M: usize, const N: usize> Div<T> for Matrix<T, M, N>
where
    T: Float + Div<Output = T> + Copy,
{
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
impl<T, const M: usize, const N: usize> Neg for Matrix<T, M, N>
where
    T: Float + Neg<Output = T> + Copy,
{
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

// Scalar multiplication
impl<U, S, const M: usize, const N: usize> Mul<S> for Matrix<U, M, N>
where
    U: Float + Mul<S, Output = U> + Copy,
    S: Copy,
    S: num_traits::Num + core::marker::Sized,
    // Prevent overlap: S must not be a Matrix
    // This uses a negative trait bound, which is not yet stable in Rust,
    // so we use a trait bound that will not be satisfied for Matrix types.
    // For practical purposes, we can add a bound that S: num_traits::Num, which Matrix does not implement.
{
    type Output = Matrix<U, M, N>;

    fn mul(self, rhs: S) -> Self::Output {
        let mut result = self;
        for r in 0..M {
            for c in 0..N {
                result.data[r][c] = self.data[r][c] * rhs;
            }
        }
        result
    }
}

// Matrix multiplication
impl<T, const M: usize, const N: usize, const P: usize> Mul<Matrix<T, N, P>> for Matrix<T, M, N>
where
    T: Float + Mul<Output = T> + Add<Output = T> + Copy + num_traits::Zero,
{
    type Output = Matrix<T, M, P>;

    fn mul(self, rhs: Matrix<T, N, P>) -> Matrix<T, M, P> {
        let mut result = Matrix {
            data: [[T::zero(); P]; M],
        };

        for i in 0..M {
            for j in 0..P {
                let mut sum = T::zero();
                for k in 0..N {
                    sum = sum + self.data[i][k] * rhs.data[k][j];
                }
                result.data[i][j] = sum;
            }
        }

        result
    }
}

// Matrix × Vector multiplication
impl<T, const M: usize, const N: usize> Mul<Vector<T, N>> for Matrix<T, M, N>
where
    T: Float + Mul<Output = T> + Add<Output = T> + Copy + num_traits::Zero,
{
    type Output = Vector<T, M>;

    fn mul(self, rhs: Vector<T, N>) -> Vector<T, M> {
        let mut result = Vector {
            data: [T::zero(); M],
        };

        for i in 0..M {
            let mut sum = T::zero();
            for j in 0..N {
                sum = sum + self.data[i][j] * rhs.data[j];
            }
            result.data[i] = sum;
        }

        result
    }
}

// Generic Matrix Implementations
impl<T, const M: usize, const N: usize> Matrix<T, M, N>
where
    T: Float + Copy,
{
    pub fn zeros() -> Self {
        Self {
            data: [[T::zero(); N]; M],
        }
    }
    pub fn ones() -> Self {
        Self {
            data: [[T::one(); N]; M],
        }
    }
}

impl<T: Float, const N: usize> Matrix<T, N, N> {
    /// Identity matrix (NxN)
    #[inline]
    pub fn identity() -> Self {
        let mut data = [[T::zero(); N]; N];
        for i in 0..N {
            data[i][i] = T::one();
        }
        Self { data }
    }

    /// Diagonal matrix from a slice/array of length N
    #[inline]
    pub fn diag(diag: &[T; N]) -> Self {
        let mut data = [[T::zero(); N]; N];
        for i in 0..N {
            data[i][i] = diag[i];
        }
        Self { data }
    }
    pub fn diag_from_vector(diag: &Vector<T, N>) -> Self {
        let mut data = [[T::zero(); N]; N];
        for i in 0..N {
            data[i][i] = diag[i];
        }
        Self { data }
    }
}

// ---------- outer product ----------
impl<F: Float, const M: usize, const N: usize> Matrix<F, M, N> {
    /// Outer (dyadic) product: a (Mx1) * b^T (1xN) -> (MxN)
    #[inline]
    pub fn outer(a: &Vector<F, M>, b: &Vector<F, N>) -> Self {
        let mut m = Self::zeros();
        for i in 0..M {
            for j in 0..N {
                m.data[i][j] = a[i] * b[j];
            }
        }
        m
    }
}

impl<T: Copy + core::ops::AddAssign + Default + num_traits::Zero, const N: usize> Matrix<T, N, N> {
    pub fn trace(&self) -> T {
        let mut sum = T::zero();
        for i in 0..N {
            sum += self.data[i][i];
        }
        sum
    }
}

// Determinant: generic implementation for any square matrix size.
impl<T: Float + Copy, const N: usize> Matrix<T, N, N> {
    /// Compute the determinant using an LU-style elimination with partial pivoting.
    /// Returns zero if the matrix is (near-)singular up to machine epsilon.
    pub fn determinant(&self) -> T {
        // make a working copy
        let mut a = *self;
        let mut det = T::one();

        for k in 0..N {
            // pivot selection (partial pivoting)
            let mut piv = k;
            let mut piv_val = a[(k, k)].abs();
            for r in (k + 1)..N {
                let v = a[(r, k)].abs();
                if v > piv_val {
                    piv = r;
                    piv_val = v;
                }
            }

            if piv_val <= T::epsilon() {
                // singular (or numerically zero)
                return T::zero();
            }

            // swap rows if needed
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
impl<T: Float + Copy + Default> Matrix<T, 3, 3> {
    pub fn inverse(&self) -> Option<Self> {
        let m = &self.data;

        let det = self.determinant();
        if det.abs() <= T::epsilon() {
            return None; // Singular matrix
        }

        let inv_det = T::one() / det;

        let mut inv = [[T::zero(); 3]; 3];

        inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
        inv[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * inv_det;
        inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;

        inv[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * inv_det;
        inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
        inv[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * inv_det;

        inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
        inv[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * inv_det;
        inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;

        Some(Matrix::new(inv))
    }
}

impl<T: Float + Copy> Matrix<T, 4, 4> {
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
impl<T: Float + Copy, const M: usize, const N: usize> Matrix<T, M, N> {
    pub fn transpose(&self) -> Matrix<T, N, M> {
        let mut out = Matrix::<T, N, M>::zeros();
        for r in 0..M {
            for c in 0..N {
                out[(c, r)] = self[(r, c)];
            }
        }
        out
    }
}

pub fn mat_block_set<T: Copy, const M: usize, const N: usize, const RM: usize, const CN: usize>(
    src: &Matrix<T, RM, CN>,
    r0: usize,
    c0: usize,
    dst: &mut Matrix<T, M, N>,
) {
    for r in 0..RM {
        for c in 0..CN {
            dst[(r0 + r, c0 + c)] = src[(r, c)];
        }
    }
}

impl<T: Float + Copy, const M: usize, const N: usize> Matrix<T, M, N> {
    pub fn max_abs(&self) -> T {
        let mut m = T::zero();
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

// Numerical Methods
impl<T: Float + Copy, const N: usize> Matrix<T, N, N> {
    /// Invert via Gauss–Jordan with partial pivoting.
    /// Returns None if (near-)singular.
    pub fn inverse_gauss_jordan(&self) -> Option<Self> {
        let mut a = *self; // working copy of A
        let mut i = Matrix::<T, N, N>::identity(); // identity → inverse

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
            if piv_val <= T::epsilon() {
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
            let invd = T::one() / a[(k, k)];
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
                if f != T::zero() {
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

// Behavior
use core::ops::{Index, IndexMut};
// [m][n]
impl<T, const M: usize, const N: usize> Index<usize> for Matrix<T, M, N> {
    type Output = [T; N];

    fn index(&self, row: usize) -> &Self::Output {
        &self.data[row]
    }
}

impl<T, const M: usize, const N: usize> IndexMut<usize> for Matrix<T, M, N> {
    fn index_mut(&mut self, row: usize) -> &mut Self::Output {
        &mut self.data[row]
    }
}
// [(m, n)]
impl<T, const M: usize, const N: usize> Index<(usize, usize)> for Matrix<T, M, N> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (r, c) = index;
        &self.data[r][c]
    }
}

impl<T, const M: usize, const N: usize> IndexMut<(usize, usize)> for Matrix<T, M, N> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (r, c) = index;
        &mut self.data[r][c]
    }
}
use num_traits::{NumCast, ToPrimitive};
impl<T: ToPrimitive + Copy, const M: usize, const N: usize> Matrix<T, M, N> {
    pub fn try_cast<U: NumCast + Copy>(self) -> Option<Matrix<U, M, N>> {
        let mut out = [[U::from(0.0)?; N]; M]; // seed (any U value works)
        for r in 0..M {
            for c in 0..N {
                out[r][c] = U::from(self.data[r][c])?;
            }
        }
        Some(Matrix { data: out })
    }

    pub fn cast<U: NumCast + Copy>(self) -> Matrix<U, M, N> {
        self.try_cast()
            .expect("Matrix::cast: not able to cast matrix....")
    }
}

// std
#[cfg(feature = "std")]
impl<T, const M: usize, const N: usize> std::fmt::Display for Matrix<T, M, N>
where
    T: Float + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[")?;
        for (i, row) in self.data.iter().enumerate() {
            for (j, val) in row.iter().enumerate() {
                write!(f, "{}", val)?;
                if j < N - 1 {
                    write!(f, " ")?; // space between columns
                }
            }
            if i < M - 1 {
                write!(f, "; ")?; // semicolon + space between rows
            }
        }
        write!(f, "]")
    }
}
