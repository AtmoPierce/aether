use super::vector::Vector;
use num_traits::{Float, Zero, zero};
use core::ops::{Add, Sub, Mul, Div, Neg};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix<T, const M: usize, const N: usize> {
    pub data: [[T; N]; M],
}

impl<T: Copy, const M: usize, const N: usize> Matrix<T, M, N> {
    pub fn new(data: [[T; N]; M]) -> Self {
        Self { data }
    }
}
impl<T: Default + Copy, const M: usize, const N: usize> Default for Matrix<T, M, N> {
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
        let mut result = Vector { data: [T::zero(); M] };

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
    T: Float + Default + Copy,
{
    pub fn zeros() -> Self {
        Self { data: [[T::zero(); N]; M] }
    }
    pub fn ones() -> Self {
        Self { data: [[T::one(); N]; M] }
    }
}

impl<T, const N: usize> Matrix<T, N, N>
where
    T: Float + Default + Copy,
{
    // Returns an identity matrix
    pub fn identity() -> Self {
        let mut data = [[T::zero(); N]; N];
        for i in 0..N {
            data[i][i] = T::one();
        }
        Self { data }
    }
    // Returns a diagonal matrix with the elements of `diag` on the diagonal.
    pub fn diag(diag: &[T; N]) -> Self {
        let mut data = [[T::zero(); N]; N];
        for i in 0..N {
            data[i][i] = diag[i];
        }
        Self { data }
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

// Determinant
impl<T: Float + Copy> Matrix<T, 2, 2> {
    pub fn determinant(&self) -> T {
        self.data[0][0] * self.data[1][1] - self.data[0][1] * self.data[1][0]
    }
}
impl<T: Float + Copy> Matrix<T, 3, 3> {
    pub fn determinant(&self) -> T {
        let m = &self.data;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
      - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
      + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }
}
impl<T: Float + Copy + Default> Matrix<T, 3, 3> {
    pub fn transpose(&self) -> Self {
        let mut result = Matrix::zeros();
        for i in 0..3 {
            for j in 0..3 {
                result[(i, j)] = self[(j, i)];
            }
        }
        result
    }

    pub fn inverse(&self) -> Option<Self> {
        let m = &self.data;

        let det = self.determinant();
        if det.abs() <= T::epsilon() {
            return None; // Singular matrix
        }

        let inv_det = T::one() / det;

        let mut inv = [[T::zero(); 3]; 3];

        inv[0][0] =  (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
        inv[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * inv_det;
        inv[0][2] =  (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;

        inv[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * inv_det;
        inv[1][1] =  (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
        inv[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * inv_det;

        inv[2][0] =  (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
        inv[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * inv_det;
        inv[2][2] =  (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;

        Some(Matrix::new(inv))
    }
}

// Behavior
use core::ops::{Index, IndexMut};
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

// std
#[cfg(feature = "std")]
use std::fmt;
#[cfg(feature = "std")]
impl<T, const M: usize, const N: usize> fmt::Display for Matrix<T, M, N>
where
    T: Float + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
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

