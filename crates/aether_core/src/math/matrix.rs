use super::algorithms::MatrixAlgorithms;
#[cfg(all(feature = "simd", any(target_arch = "x86", target_arch = "x86_64")))]
use super::arch::x86::matrix_simd;
#[cfg(all(feature = "simd", target_arch = "aarch64"))]
use super::arch::arm::neon as arm_neon;
#[cfg(all(feature = "simd", target_arch = "arm"))]
use super::arch::arm::m33_dsp;
use super::vector::Vector;
use crate::real::{Real};
use core::ops::{Add, Sub, Mul, Div, Neg, AddAssign};
use core::ops::{Index, IndexMut};
use core::mem::MaybeUninit;
use core::ptr;

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
        self.mul_matrix_strided(&rhs)
    }
}

impl<T: Real, const M: usize, const N: usize, const P: usize> Mul<&Matrix<T, N, P>>
    for &Matrix<T, M, N>
{
    type Output = Matrix<T, M, P>;

    fn mul(self, rhs: &Matrix<T, N, P>) -> Matrix<T, M, P> {
        self.mul_matrix_strided(rhs)
    }
}

impl<T: Real, const M: usize, const N: usize, const P: usize> Mul<&Matrix<T, N, P>>
    for Matrix<T, M, N>
{
    type Output = Matrix<T, M, P>;

    fn mul(self, rhs: &Matrix<T, N, P>) -> Matrix<T, M, P> {
        (&self) * rhs
    }
}

impl<T: Real, const M: usize, const N: usize, const P: usize> Mul<Matrix<T, N, P>>
    for &Matrix<T, M, N>
{
    type Output = Matrix<T, M, P>;

    fn mul(self, rhs: Matrix<T, N, P>) -> Matrix<T, M, P> {
        self * (&rhs)
    }
}

// Matrix × Vector: (M×N) * (N) -> (M)
impl<T: Real, const M: usize, const N: usize> Mul<Vector<T, N>> for Matrix<T, M, N> {
    type Output = Vector<T, M>;

    fn mul(self, rhs: Vector<T, N>) -> Vector<T, M> {
        self.mul_vector_unrolled(&rhs)
    }
}

impl<T: Real, const M: usize, const N: usize> Mul<&Vector<T, N>> for &Matrix<T, M, N> {
    type Output = Vector<T, M>;

    fn mul(self, rhs: &Vector<T, N>) -> Vector<T, M> {
        self.mul_vector_unrolled(rhs)
    }
}

impl<T: Real, const M: usize, const N: usize> Mul<&Vector<T, N>> for Matrix<T, M, N> {
    type Output = Vector<T, M>;

    fn mul(self, rhs: &Vector<T, N>) -> Vector<T, M> {
        (&self) * rhs
    }
}

impl<T: Real, const M: usize, const N: usize> Mul<Vector<T, N>> for &Matrix<T, M, N> {
    type Output = Vector<T, M>;

    fn mul(self, rhs: Vector<T, N>) -> Vector<T, M> {
        self * (&rhs)
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

impl<const M: usize> Matrix<f64, M, 6> {
    #[inline(always)]
    pub fn mul_vec6_simd(&self, rhs: &Vector<f64, 6>) -> Vector<f64, M> {
        #[cfg(all(
            feature = "simd",
            feature = "fma",
            feature = "std",
            any(target_arch = "x86", target_arch = "x86_64")
        ))]
        {
            if std::is_x86_feature_detected!("fma") && std::is_x86_feature_detected!("sse2") {
                unsafe {
                    return matrix_simd::mul_vec6_fma128_f64(self, rhs);
                }
            }
        }

        #[cfg(all(
            feature = "simd",
            feature = "fma",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma",
            target_feature = "sse2"
        ))]
        unsafe {
            return matrix_simd::mul_vec6_fma128_f64(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            feature = "std",
            any(target_arch = "x86", target_arch = "x86_64")
        ))]
        {
            if std::is_x86_feature_detected!("avx") && std::is_x86_feature_detected!("fma") {
                unsafe {
                    return matrix_simd::mul_vec6_avx_fma_f64(self, rhs);
                }
            }

            if std::is_x86_feature_detected!("avx") {
                unsafe {
                    return matrix_simd::mul_vec6_avx_f64(self, rhs);
                }
            }

            if std::is_x86_feature_detected!("sse2") {
                unsafe {
                    return matrix_simd::mul_vec6_sse2_f64(self, rhs);
                }
            }
        }

        #[cfg(all(
            feature = "simd",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "avx",
            target_feature = "fma"
        ))]
        unsafe {
            return matrix_simd::mul_vec6_avx_fma_f64(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "avx"
        ))]
        unsafe {
            return matrix_simd::mul_vec6_avx_f64(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "sse2"
        ))]
        unsafe {
            return matrix_simd::mul_vec6_sse2_f64(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_vec6_scalar(rhs)
        }
    }

    #[inline(always)]
    fn mul_vec6_scalar(&self, rhs: &Vector<f64, 6>) -> Vector<f64, M> {
        let mut result = Vector { data: [0.0; M] };
        for i in 0..M {
            let row = &self.data[i];
            result.data[i] = row[0] * rhs.data[0]
                + row[1] * rhs.data[1]
                + row[2] * rhs.data[2]
                + row[3] * rhs.data[3]
                + row[4] * rhs.data[4]
                + row[5] * rhs.data[5];
        }
        result
    }
}

impl<const M: usize> Matrix<f64, M, 4> {
    #[inline(always)]
    pub fn mul_vec4_simd(&self, rhs: &Vector<f64, 4>) -> Vector<f64, M> {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::mul_vec4_neon_f64(self, rhs);
        }

        #[cfg(all(feature = "simd", target_arch = "arm", target_feature = "dsp"))]
        {
            return m33_dsp::mul_vec4_m33_f64(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            feature = "fma",
            feature = "std",
            any(target_arch = "x86", target_arch = "x86_64")
        ))]
        {
            if std::is_x86_feature_detected!("avx") && std::is_x86_feature_detected!("fma") {
                unsafe {
                    return matrix_simd::mul_vec4_avx_fma_f64(self, rhs);
                }
            }
        }

        #[cfg(all(
            feature = "simd",
            feature = "std",
            any(target_arch = "x86", target_arch = "x86_64")
        ))]
        {
            if std::is_x86_feature_detected!("avx") {
                unsafe {
                    return matrix_simd::mul_vec4_avx_f64(self, rhs);
                }
            }
        }

        #[cfg(all(
            feature = "simd",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "avx",
            target_feature = "fma",
            feature = "fma"
        ))]
        unsafe {
            return matrix_simd::mul_vec4_avx_fma_f64(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "avx"
        ))]
        unsafe {
            return matrix_simd::mul_vec4_avx_f64(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_vec4_scalar(rhs)
        }
    }

    #[inline(always)]
    fn mul_vec4_scalar(&self, rhs: &Vector<f64, 4>) -> Vector<f64, M> {
        let mut result = Vector { data: [0.0; M] };
        for i in 0..M {
            let row = &self.data[i];
            result.data[i] = row[0] * rhs.data[0]
                + row[1] * rhs.data[1]
                + row[2] * rhs.data[2]
                + row[3] * rhs.data[3];
        }
        result
    }
}

impl<const M: usize> Matrix<f32, M, 6> {
    #[inline(always)]
    pub fn mul_vec6_simd(&self, rhs: &Vector<f32, 6>) -> Vector<f32, M> {
        #[cfg(all(
            feature = "simd",
            feature = "fma",
            feature = "std",
            any(target_arch = "x86", target_arch = "x86_64")
        ))]
        {
            if std::is_x86_feature_detected!("fma") && std::is_x86_feature_detected!("sse") {
                unsafe {
                    return matrix_simd::mul_vec6_fma128_f32(self, rhs);
                }
            }
        }

        #[cfg(all(
            feature = "simd",
            feature = "fma",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma",
            target_feature = "sse"
        ))]
        unsafe {
            return matrix_simd::mul_vec6_fma128_f32(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            feature = "std",
            any(target_arch = "x86", target_arch = "x86_64")
        ))]
        {
            if std::is_x86_feature_detected!("avx") && std::is_x86_feature_detected!("fma") {
                unsafe {
                    return matrix_simd::mul_vec6_avx_fma_f32(self, rhs);
                }
            }

            if std::is_x86_feature_detected!("avx") {
                unsafe {
                    return matrix_simd::mul_vec6_avx_f32(self, rhs);
                }
            }

            if std::is_x86_feature_detected!("sse") {
                unsafe {
                    return matrix_simd::mul_vec6_sse_f32(self, rhs);
                }
            }
        }

        #[cfg(all(
            feature = "simd",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "avx",
            target_feature = "fma"
        ))]
        unsafe {
            return matrix_simd::mul_vec6_avx_fma_f32(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "avx"
        ))]
        unsafe {
            return matrix_simd::mul_vec6_avx_f32(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "sse"
        ))]
        unsafe {
            return matrix_simd::mul_vec6_sse_f32(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_vec6_scalar(rhs)
        }
    }

    #[inline(always)]
    fn mul_vec6_scalar(&self, rhs: &Vector<f32, 6>) -> Vector<f32, M> {
        let mut result = Vector { data: [0.0; M] };
        for i in 0..M {
            let row = &self.data[i];
            result.data[i] = row[0] * rhs.data[0]
                + row[1] * rhs.data[1]
                + row[2] * rhs.data[2]
                + row[3] * rhs.data[3]
                + row[4] * rhs.data[4]
                + row[5] * rhs.data[5];
        }
        result
    }
}

impl<const M: usize> Matrix<f32, M, 4> {
    #[inline(always)]
    pub fn mul_vec4_simd(&self, rhs: &Vector<f32, 4>) -> Vector<f32, M> {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::mul_vec4_neon_f32(self, rhs);
        }

        #[cfg(all(feature = "simd", target_arch = "arm", target_feature = "dsp"))]
        {
            return m33_dsp::mul_vec4_m33_f32(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            feature = "fma",
            feature = "std",
            any(target_arch = "x86", target_arch = "x86_64")
        ))]
        {
            if std::is_x86_feature_detected!("fma") && std::is_x86_feature_detected!("sse") {
                unsafe {
                    return matrix_simd::mul_vec4_sse_fma_f32(self, rhs);
                }
            }
        }

        #[cfg(all(
            feature = "simd",
            feature = "std",
            any(target_arch = "x86", target_arch = "x86_64")
        ))]
        {
            if std::is_x86_feature_detected!("sse") {
                unsafe {
                    return matrix_simd::mul_vec4_sse_f32(self, rhs);
                }
            }
        }

        #[cfg(all(
            feature = "simd",
            feature = "fma",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "sse",
            target_feature = "fma"
        ))]
        unsafe {
            return matrix_simd::mul_vec4_sse_fma_f32(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "sse"
        ))]
        unsafe {
            return matrix_simd::mul_vec4_sse_f32(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_vec4_scalar(rhs)
        }
    }

    #[inline(always)]
    fn mul_vec4_scalar(&self, rhs: &Vector<f32, 4>) -> Vector<f32, M> {
        let mut result = Vector { data: [0.0; M] };
        for i in 0..M {
            let row = &self.data[i];
            result.data[i] = row[0] * rhs.data[0]
                + row[1] * rhs.data[1]
                + row[2] * rhs.data[2]
                + row[3] * rhs.data[3];
        }
        result
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

impl<const M: usize, const N: usize> Matrix<f32, M, N> {
    #[inline(always)]
    pub fn mul_matrix_simd<const P: usize>(&self, rhs: &Matrix<f32, N, P>) -> Matrix<f32, M, P> {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::mul_matrix_neon_f32(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_matrix_strided(rhs)
        }
    }
}

impl<const M: usize, const N: usize> Matrix<f64, M, N> {
    #[inline(always)]
    pub fn mul_matrix_simd<const P: usize>(&self, rhs: &Matrix<f64, N, P>) -> Matrix<f64, M, P> {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::mul_matrix_neon_f64(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_matrix_strided(rhs)
        }
    }
}

impl Matrix<f32, 3, 3> {
    #[inline(always)]
    pub fn mul_mat3_simd(&self, rhs: &Matrix<f32, 3, 3>) -> Matrix<f32, 3, 3> {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::mul_mat3_neon_f32(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_matrix_strided(rhs)
        }
    }
}

impl Matrix<f64, 3, 3> {
    #[inline(always)]
    pub fn mul_mat3_simd(&self, rhs: &Matrix<f64, 3, 3>) -> Matrix<f64, 3, 3> {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::mul_mat3_neon_f64(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_matrix_strided(rhs)
        }
    }
}

impl Matrix<f32, 4, 4> {
    #[inline(always)]
    pub fn mul_mat4_simd(&self, rhs: &Matrix<f32, 4, 4>) -> Matrix<f32, 4, 4> {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::mul_mat4_neon_f32(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_matrix_strided(rhs)
        }
    }
}

impl Matrix<f64, 4, 4> {
    #[inline(always)]
    pub fn mul_mat4_simd(&self, rhs: &Matrix<f64, 4, 4>) -> Matrix<f64, 4, 4> {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::mul_mat4_neon_f64(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_matrix_strided(rhs)
        }
    }
}

impl Matrix<f32, 6, 6> {
    #[inline(always)]
    pub fn mul_mat6_simd(&self, rhs: &Matrix<f32, 6, 6>) -> Matrix<f32, 6, 6> {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::mul_mat6_neon_f32(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_matrix_strided(rhs)
        }
    }
}

impl Matrix<f64, 6, 6> {
    #[inline(always)]
    pub fn mul_mat6_simd(&self, rhs: &Matrix<f64, 6, 6>) -> Matrix<f64, 6, 6> {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::mul_mat6_neon_f64(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.mul_matrix_strided(rhs)
        }
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

/* Serialization and Deserialization */
#[cfg(feature = "bincode")]
impl<T, const M: usize, const N: usize> bincode::Encode for Matrix<T, M, N>
where
    T: bincode::Encode,
{
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        for r in 0..M {
            for c in 0..N {
                self.data[r][c].encode(encoder)?;
            }
        }
        Ok(())
    }
}

#[cfg(feature = "bincode")]
impl<T, const M: usize, const N: usize, Ctx> bincode::Decode<Ctx> for Matrix<T, M, N>
where
    T: bincode::Decode<Ctx>,
{
    fn decode<D: bincode::de::Decoder<Context = Ctx>>(
        decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        let mut out = MaybeUninit::<[[T; N]; M]>::uninit();
        let out_ptr = out.as_mut_ptr() as *mut [T; N];
        let mut inited_rows = 0usize;

        for r in 0..M {
            let mut row = MaybeUninit::<[T; N]>::uninit();
            let row_ptr = row.as_mut_ptr() as *mut T;
            let mut inited_cols = 0usize;

            for c in 0..N {
                match T::decode(decoder) {
                    Ok(v) => {
                        unsafe { row_ptr.add(c).write(v) };
                        inited_cols += 1;
                    }
                    Err(e) => {
                        unsafe {
                            for i in 0..inited_cols {
                                ptr::drop_in_place(row_ptr.add(i));
                            }
                            for i in 0..inited_rows {
                                ptr::drop_in_place(out_ptr.add(i));
                            }
                        }
                        return Err(e);
                    }
                }
            }

            unsafe { out_ptr.add(r).write(row.assume_init()) };
            inited_rows += 1;
        }

        Ok(Matrix { data: unsafe { out.assume_init() } })
    }
}

#[cfg(feature = "bincode")]
impl<'de, T, const M: usize, const N: usize, Ctx> bincode::BorrowDecode<'de, Ctx>
    for Matrix<T, M, N>
where
    T: bincode::BorrowDecode<'de, Ctx>,
{
    fn borrow_decode<D: bincode::de::BorrowDecoder<'de, Context = Ctx>>(
        decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        let mut out = MaybeUninit::<[[T; N]; M]>::uninit();
        let out_ptr = out.as_mut_ptr() as *mut [T; N];
        let mut inited_rows = 0usize;

        for r in 0..M {
            let mut row = MaybeUninit::<[T; N]>::uninit();
            let row_ptr = row.as_mut_ptr() as *mut T;
            let mut inited_cols = 0usize;

            for c in 0..N {
                match T::borrow_decode(decoder) {
                    Ok(v) => {
                        unsafe { row_ptr.add(c).write(v) };
                        inited_cols += 1;
                    }
                    Err(e) => {
                        unsafe {
                            for i in 0..inited_cols {
                                ptr::drop_in_place(row_ptr.add(i));
                            }
                            for i in 0..inited_rows {
                                ptr::drop_in_place(out_ptr.add(i));
                            }
                        }
                        return Err(e);
                    }
                }
            }

            unsafe { out_ptr.add(r).write(row.assume_init()) };
            inited_rows += 1;
        }

        Ok(Matrix { data: unsafe { out.assume_init() } })
    }
}

#[cfg(feature = "serde")]
use core::marker::PhantomData;
#[cfg(feature = "serde")]
use serde::ser::SerializeSeq;

#[cfg(feature = "serde")]
impl<T, const M: usize, const N: usize> serde::Serialize for Matrix<T, M, N>
where
    T: serde::Serialize,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        struct RowRef<'a, T, const N: usize>(&'a [T; N]);

        impl<'a, T, const N: usize> serde::Serialize for RowRef<'a, T, N>
        where
            T: serde::Serialize,
        {
            fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
            where
                S: serde::Serializer,
            {
                let mut seq = serializer.serialize_seq(Some(N))?;
                for v in self.0.iter() {
                    seq.serialize_element(v)?;
                }
                seq.end()
            }
        }

        let mut seq = serializer.serialize_seq(Some(M))?;
        for row in self.data.iter() {
            seq.serialize_element(&RowRef::<T, N>(row))?;
        }
        seq.end()
    }
}

#[cfg(feature = "serde")]
impl<'de, T, const M: usize, const N: usize> serde::Deserialize<'de> for Matrix<T, M, N>
where
    T: serde::Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        struct RowSeed<T, const N: usize>(PhantomData<T>);

        impl<'de, T, const N: usize> serde::de::DeserializeSeed<'de> for RowSeed<T, N>
        where
            T: serde::Deserialize<'de>,
        {
            type Value = [T; N];

            fn deserialize<D>(self, deserializer: D) -> Result<Self::Value, D::Error>
            where
                D: serde::Deserializer<'de>,
            {
                struct RowVisitor<T, const N: usize>(PhantomData<T>);

                impl<'de, T, const N: usize> serde::de::Visitor<'de> for RowVisitor<T, N>
                where
                    T: serde::Deserialize<'de>,
                {
                    type Value = [T; N];

                    fn expecting(&self, formatter: &mut core::fmt::Formatter) -> core::fmt::Result {
                        write!(formatter, "a row with {} elements", N)
                    }

                    fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
                    where
                        A: serde::de::SeqAccess<'de>,
                    {
                        let mut row = MaybeUninit::<[T; N]>::uninit();
                        let row_ptr = row.as_mut_ptr() as *mut T;
                        let mut inited = 0usize;

                        for i in 0..N {
                            match seq.next_element()? {
                                Some(v) => {
                                    unsafe { row_ptr.add(i).write(v) };
                                    inited += 1;
                                }
                                None => {
                                    unsafe {
                                        for j in 0..inited {
                                            ptr::drop_in_place(row_ptr.add(j));
                                        }
                                    }
                                    return Err(serde::de::Error::invalid_length(i, &self));
                                }
                            }
                        }

                        if let Some(_) = seq.next_element::<serde::de::IgnoredAny>()? {
                            unsafe {
                                for j in 0..inited {
                                    ptr::drop_in_place(row_ptr.add(j));
                                }
                            }
                            return Err(serde::de::Error::invalid_length(N + 1, &self));
                        }

                        Ok(unsafe { row.assume_init() })
                    }
                }

                deserializer.deserialize_seq(RowVisitor::<T, N>(PhantomData))
            }
        }

        struct MatrixVisitor<T, const M: usize, const N: usize>(PhantomData<T>);

        impl<'de, T, const M: usize, const N: usize> serde::de::Visitor<'de>
            for MatrixVisitor<T, M, N>
        where
            T: serde::Deserialize<'de>,
        {
            type Value = Matrix<T, M, N>;

            fn expecting(&self, formatter: &mut core::fmt::Formatter) -> core::fmt::Result {
                write!(formatter, "a {M}x{N} matrix")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut out = MaybeUninit::<[[T; N]; M]>::uninit();
                let out_ptr = out.as_mut_ptr() as *mut [T; N];
                let mut inited = 0usize;

                for i in 0..M {
                    match seq.next_element_seed(RowSeed::<T, N>(PhantomData))? {
                        Some(row) => {
                            unsafe { out_ptr.add(i).write(row) };
                            inited += 1;
                        }
                        None => {
                            unsafe {
                                for j in 0..inited {
                                    ptr::drop_in_place(out_ptr.add(j));
                                }
                            }
                            return Err(serde::de::Error::invalid_length(i, &self));
                        }
                    }
                }

                if let Some(_) = seq.next_element::<serde::de::IgnoredAny>()? {
                    unsafe {
                        for j in 0..inited {
                            ptr::drop_in_place(out_ptr.add(j));
                        }
                    }
                    return Err(serde::de::Error::invalid_length(M + 1, &self));
                }

                Ok(Matrix { data: unsafe { out.assume_init() } })
            }
        }

        deserializer.deserialize_seq(MatrixVisitor::<T, M, N>(PhantomData))
    }
}
