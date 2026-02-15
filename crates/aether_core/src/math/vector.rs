use crate::utils::{ToDegrees, ToRadians};
use crate::math::Matrix;
use crate::real::Real;
use super::algorithms::VectorAlgorithms;
#[cfg(all(feature = "simd", any(target_arch = "x86", target_arch = "x86_64")))]
use super::arch::x86::vector_simd;
#[cfg(all(feature = "simd", target_arch = "aarch64"))]
use super::arch::arm::neon as arm_neon;
#[cfg(all(feature = "simd", target_arch = "arm"))]
use super::arch::arm::m33_dsp;

use core::array::IntoIter as ArrayIntoIter;
use core::iter::FromIterator;
use core::ops::{
    Add, AddAssign, Div, Mul, Neg, Sub, SubAssign, Index, IndexMut
};
use core::slice::{Iter, IterMut};

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub struct Vector<T, const N: usize> {
    pub data: [T; N],
}

/* -------------------- Constructors / Default -------------------- */

impl<T: Copy, const N: usize> Vector<T, N> {
    #[inline]
    pub fn new(data: [T; N]) -> Self {
        Self { data }
    }
}

impl<T: Default + Copy, const N: usize> Default for Vector<T, N> {
    fn default() -> Self {
        Self {
            data: [T::default(); N],
        }
    }
}

/* -------------------- Add -------------------- */

// Vector + Vector
impl<T: Add<Output = T> + Copy, const N: usize> Add for Vector<T, N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let mut result = self;
        for i in 0..N {
            result.data[i] = self.data[i] + rhs.data[i];
        }
        result
    }
}

// &Vector + &Vector
impl<'a, T: Add<Output = T> + Copy, const N: usize> Add for &'a Vector<T, N> {
    type Output = Vector<T, N>;
    fn add(self, rhs: Self) -> Self::Output {
        let mut result = *self;
        for i in 0..N {
            result.data[i] = self.data[i] + rhs.data[i];
        }
        result
    }
}

// Vector + &Vector
impl<'a, T: Add<Output = T> + Copy, const N: usize> Add<&'a Vector<T, N>> for Vector<T, N> {
    type Output = Vector<T, N>;
    fn add(self, rhs: &'a Vector<T, N>) -> Self::Output {
        self + *rhs
    }
}

// &Vector + Vector
impl<'a, T: Add<Output = T> + Copy, const N: usize> Add<Vector<T, N>> for &'a Vector<T, N> {
    type Output = Vector<T, N>;
    fn add(self, rhs: Vector<T, N>) -> Self::Output {
        *self + rhs
    }
}

// AddAssign: Vector += Vector
impl<T: AddAssign + Copy, const N: usize> AddAssign for Vector<T, N> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.data[i] += rhs.data[i];
        }
    }
}

// AddAssign: Vector += &Vector
impl<T: AddAssign + Copy, const N: usize> AddAssign<&Vector<T, N>> for Vector<T, N> {
    fn add_assign(&mut self, rhs: &Vector<T, N>) {
        for i in 0..N {
            self.data[i] += rhs.data[i];
        }
    }
}

/* -------------------- Sub -------------------- */

// Vector - Vector
impl<T: Sub<Output = T> + Copy, const N: usize> Sub for Vector<T, N> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let mut result = self;
        for i in 0..N {
            result.data[i] = self.data[i] - rhs.data[i];
        }
        result
    }
}

// &Vector - &Vector
impl<'a, T: Sub<Output = T> + Copy, const N: usize> Sub for &'a Vector<T, N> {
    type Output = Vector<T, N>;
    fn sub(self, rhs: Self) -> Self::Output {
        let mut result = *self;
        for i in 0..N {
            result.data[i] = self.data[i] - rhs.data[i];
        }
        result
    }
}

// Vector - &Vector
impl<'a, T: Sub<Output = T> + Copy, const N: usize> Sub<&'a Vector<T, N>> for Vector<T, N> {
    type Output = Vector<T, N>;
    fn sub(self, rhs: &'a Vector<T, N>) -> Self::Output {
        self - *rhs
    }
}

// &Vector - Vector
impl<'a, T: Sub<Output = T> + Copy, const N: usize> Sub<Vector<T, N>> for &'a Vector<T, N> {
    type Output = Vector<T, N>;
    fn sub(self, rhs: Vector<T, N>) -> Self::Output {
        *self - rhs
    }
}

// SubAssign: Vector -= Vector
impl<T: SubAssign + Copy, const N: usize> SubAssign for Vector<T, N> {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.data[i] -= rhs.data[i];
        }
    }
}

// SubAssign: Vector -= &Vector
impl<T: SubAssign + Copy, const N: usize> SubAssign<&Vector<T, N>> for Vector<T, N> {
    fn sub_assign(&mut self, rhs: &Vector<T, N>) {
        for i in 0..N {
            self.data[i] -= rhs.data[i];
        }
    }
}

/* -------------------- Neg -------------------- */

// -Vector
impl<T: Neg<Output = T> + Copy, const N: usize> Neg for Vector<T, N> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        let mut result = self;
        for i in 0..N {
            result.data[i] = -self.data[i];
        }
        result
    }
}

// -&Vector
impl<'a, T: Neg<Output = T> + Copy, const N: usize> Neg for &'a Vector<T, N> {
    type Output = Vector<T, N>;
    fn neg(self) -> Self::Output {
        let mut result = *self;
        for i in 0..N {
            result.data[i] = -self.data[i];
        }
        result
    }
}

/* -------------------- Scalar mul/div -------------------- */

// Vector * scalar
impl<T, const N: usize> Mul<T> for Vector<T, N>
where
    T: Mul<Output = T> + Copy,
{
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        let mut result = self;
        for i in 0..N {
            result.data[i] = self.data[i] * rhs;
        }
        result
    }
}

// Vector / scalar
impl<T, const N: usize> Div<T> for Vector<T, N>
where
    T: Div<Output = T> + Copy,
{
    type Output = Self;
    fn div(self, rhs: T) -> Self {
        let mut result = self;
        for i in 0..N {
            result.data[i] = self.data[i] / rhs;
        }
        result
    }
}

// By reference: &Vector * scalar
impl<'a, T: Real, const N: usize> Mul<T> for &'a Vector<T, N> {
    type Output = Vector<T, N>;
    fn mul(self, rhs: T) -> Self::Output {
        (*self).clone() * rhs
    }
}

// By reference: &Vector / scalar
impl<'a, T: Real, const N: usize> Div<T> for &'a Vector<T, N> {
    type Output = Vector<T, N>;
    fn div(self, rhs: T) -> Self::Output {
        (*self).clone() / rhs
    }
}

/* -------------------- Dot as operator -------------------- */

impl<T: Real, const N: usize> Mul<&Vector<T, N>> for &Vector<T, N> {
    type Output = T;

    fn mul(self, rhs: &Vector<T, N>) -> Self::Output {
        self.dot(rhs)
    }
}

impl<T: Real, const N: usize> Mul<&Vector<T, N>> for Vector<T, N> {
    type Output = T;

    fn mul(self, rhs: &Vector<T, N>) -> Self::Output {
        (&self) * rhs
    }
}

impl<T: Real, const N: usize> Mul<Vector<T, N>> for &Vector<T, N> {
    type Output = T;

    fn mul(self, rhs: Vector<T, N>) -> Self::Output {
        self * (&rhs)
    }
}

impl<T: Real, const N: usize> Mul<Vector<T, N>> for Vector<T, N> {
    type Output = T;

    fn mul(self, rhs: Vector<T, N>) -> Self::Output {
        (&self) * (&rhs)
    }
}

/* -------------------- Norms, dot, angles -------------------- */

impl<T, const N: usize> Vector<T, N>
where
    T: Real,
{
    pub fn dot(&self, rhs: &Self) -> T {
        <Self as VectorAlgorithms<T, N>>::dot_generic(self, rhs)
    }

    pub fn norm(&self) -> T {
        <Self as VectorAlgorithms<T, N>>::norm_generic(self)
    }

    pub fn normalize(&self) -> Self {
        self / self.norm()
    }

    /// Returns (normalized, norm). If norm is too small, returns (self, norm).
    pub fn try_normalize(&self) -> (Self, T) {
        let n = self.norm();
        let eps = T::from_f64(1e-12);
        if n > eps {
            (*self / n, n)
        } else {
            (*self, n)
        }
    }

    pub fn magnitude(&self) -> T {
        self.norm()
    }

    pub fn angle(&self, rhs: &Self) -> T {
        let dot = self.dot(rhs);
        let denom = self.norm() * rhs.norm();
        if denom == T::ZERO {
            T::ZERO
        } else {
            let mut cos_theta = dot / denom;
            let minus_one = T::from_f64(-1.0);
            let one = T::ONE;

            // clamp to [-1, 1] without max/min
            if cos_theta > one {
                cos_theta = one;
            }
            if cos_theta < minus_one {
                cos_theta = minus_one;
            }
            cos_theta.acos()
        }
    }
}

impl Vector<f64, 4> {
    #[inline(always)]
    pub fn dot4_simd(&self, rhs: &Self) -> f64 {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::dot4_neon_f64(self, rhs);
        }

        #[cfg(all(feature = "simd", target_arch = "arm", target_feature = "dsp"))]
        {
            return m33_dsp::dot4_m33_f64(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            feature = "std",
            any(target_arch = "x86", target_arch = "x86_64")
        ))]
        {
            type Kernel = fn(&Vector<f64, 4>, &Vector<f64, 4>) -> f64;
            static KERNEL: std::sync::OnceLock<Kernel> = std::sync::OnceLock::new();

            let kernel = KERNEL.get_or_init(|| {
                #[cfg(feature = "fma")]
                {
                    if std::is_x86_feature_detected!("avx") && std::is_x86_feature_detected!("fma") {
                        return Self::dot4_avx_fma_kernel;
                    }
                }

                if std::is_x86_feature_detected!("avx") {
                    return Self::dot4_avx_kernel;
                }

                Self::dot4_scalar_kernel
            });

            return kernel(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            feature = "fma",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "avx",
            target_feature = "fma"
        ))]
        unsafe {
            return vector_simd::dot4_avx_fma_f64(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "avx"
        ))]
        unsafe {
            return vector_simd::dot4_avx_f64(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.dot(rhs)
        }
    }

    #[inline(always)]
    fn dot4_scalar_kernel(a: &Vector<f64, 4>, b: &Vector<f64, 4>) -> f64 {
        a.dot(b)
    }

    #[cfg(all(feature = "simd", any(target_arch = "x86", target_arch = "x86_64")))]
    #[inline(always)]
    fn dot4_avx_kernel(a: &Vector<f64, 4>, b: &Vector<f64, 4>) -> f64 {
        unsafe { vector_simd::dot4_avx_f64(a, b) }
    }

    #[cfg(all(feature = "simd", any(target_arch = "x86", target_arch = "x86_64"), feature = "fma"))]
    #[inline(always)]
    fn dot4_avx_fma_kernel(a: &Vector<f64, 4>, b: &Vector<f64, 4>) -> f64 {
        unsafe { vector_simd::dot4_avx_fma_f64(a, b) }
    }
}

impl Vector<f32, 4> {
    #[inline(always)]
    pub fn dot4_simd(&self, rhs: &Self) -> f32 {
        #[cfg(all(feature = "simd", target_arch = "aarch64"))]
        unsafe {
            return arm_neon::dot4_neon_f32(self, rhs);
        }

        #[cfg(all(feature = "simd", target_arch = "arm", target_feature = "dsp"))]
        {
            return m33_dsp::dot4_m33_f32(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            feature = "std",
            any(target_arch = "x86", target_arch = "x86_64")
        ))]
        {
            type Kernel = fn(&Vector<f32, 4>, &Vector<f32, 4>) -> f32;
            static KERNEL: std::sync::OnceLock<Kernel> = std::sync::OnceLock::new();

            let kernel = KERNEL.get_or_init(|| {
                #[cfg(feature = "fma")]
                {
                    if std::is_x86_feature_detected!("fma") && std::is_x86_feature_detected!("sse") {
                        return Self::dot4_sse_fma_kernel;
                    }
                }

                if std::is_x86_feature_detected!("sse") {
                    return Self::dot4_sse_kernel;
                }

                Self::dot4_scalar_kernel
            });

            return kernel(self, rhs);
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
            return vector_simd::dot4_sse_fma_f32(self, rhs);
        }

        #[cfg(all(
            feature = "simd",
            not(feature = "std"),
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "sse"
        ))]
        unsafe {
            return vector_simd::dot4_sse_f32(self, rhs);
        }

        #[allow(unreachable_code)]
        {
            self.dot(rhs)
        }
    }

    #[inline(always)]
    fn dot4_scalar_kernel(a: &Vector<f32, 4>, b: &Vector<f32, 4>) -> f32 {
        a.dot(b)
    }

    #[cfg(all(feature = "simd", any(target_arch = "x86", target_arch = "x86_64")))]
    #[inline(always)]
    fn dot4_sse_kernel(a: &Vector<f32, 4>, b: &Vector<f32, 4>) -> f32 {
        unsafe { vector_simd::dot4_sse_f32(a, b) }
    }

    #[cfg(all(feature = "simd", any(target_arch = "x86", target_arch = "x86_64"), feature = "fma"))]
    #[inline(always)]
    fn dot4_sse_fma_kernel(a: &Vector<f32, 4>, b: &Vector<f32, 4>) -> f32 {
        unsafe { vector_simd::dot4_sse_fma_f32(a, b) }
    }
}

/* -------------------- Cross product (3D) -------------------- */

impl<T> Vector<T, 3>
where
    T: Copy + Sub<Output = T> + Mul<Output = T>,
{
    pub fn cross(&self, rhs: &Self) -> Self {
        let [a1, a2, a3] = self.data;
        let [b1, b2, b3] = rhs.data;
        Self {
            data: [
                a2 * b3 - a3 * b2,
                a3 * b1 - a1 * b3,
                a1 * b2 - a2 * b1,
            ],
        }
    }
}

/* -------------------- Units helpers -------------------- */

impl<T, const N: usize> ToRadians for Vector<T, N>
where
    T: ToRadians + Copy,
{
    type Output = Vector<T::Output, N>;
    fn to_radians(self) -> Self::Output {
        let data = self.data.map(|x| x.to_radians());
        Vector { data }
    }
}

impl<T, const N: usize> ToDegrees for Vector<T, N>
where
    T: ToDegrees + Copy,
{
    type Output = Vector<T::Output, N>;
    fn to_degrees(self) -> Self::Output {
        let data = self.data.map(|x| x.to_degrees());
        Vector { data }
    }
}

/* -------------------- Indexing -------------------- */

impl<T, const N: usize> Index<usize> for Vector<T, N> {
    type Output = T;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.data[idx]
    }
}

impl<T, const N: usize> IndexMut<usize> for Vector<T, N> {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.data[idx]
    }
}

/* -------------------- Matrix conversion -------------------- */

impl<T: Copy, const N: usize> Vector<T, N> {
    #[inline]
    pub fn row(self) -> Matrix<T, 1, N> {
        Matrix { data: [self.data] }
    }

    #[inline]
    pub fn col(self) -> Matrix<T, N, 1> {
        Matrix { data: self.data.map(|x| [x]) }
    }
}

/* -------------------- AsRef / AsMut, From, Iteration -------------------- */

impl<T, const N: usize> AsRef<[T; N]> for Vector<T, N> {
    #[inline]
    fn as_ref(&self) -> &[T; N] {
        &self.data
    }
}

impl<T, const N: usize> AsMut<[T; N]> for Vector<T, N> {
    #[inline]
    fn as_mut(&mut self) -> &mut [T; N] {
        &mut self.data
    }
}

impl<T, const N: usize> From<[T; N]> for Vector<T, N> {
    #[inline]
    fn from(data: [T; N]) -> Self {
        Self { data }
    }
}

impl<T, const N: usize> Vector<T, N> {
    #[inline]
    pub fn iter(&self) -> Iter<'_, T> {
        self.data.iter()
    }

    #[inline]
    pub fn iter_mut(&mut self) -> IterMut<'_, T> {
        self.data.iter_mut()
    }

    #[inline]
    pub fn iter_copied(&self) -> impl Iterator<Item = T> + '_
    where
        T: Copy,
    {
        self.data.iter().copied()
    }

    pub fn map<U, F: FnMut(T) -> U>(self, mut f: F) -> Vector<U, N> {
        let data = self.data.map(|x| f(x));
        Vector { data }
    }
}

// IntoIterator (by value)
impl<T, const N: usize> IntoIterator for Vector<T, N> {
    type Item = T;
    type IntoIter = ArrayIntoIter<T, N>;
    #[inline]
    fn into_iter(self) -> ArrayIntoIter<T, N> {
        self.data.into_iter()
    }
}

// IntoIterator for &Vector -> &T
impl<'a, T, const N: usize> IntoIterator for &'a Vector<T, N> {
    type Item = &'a T;
    type IntoIter = Iter<'a, T>;
    #[inline]
    fn into_iter(self) -> Iter<'a, T> {
        self.data.iter()
    }
}

// IntoIterator for &mut Vector -> &mut T
impl<'a, T, const N: usize> IntoIterator for &'a mut Vector<T, N> {
    type Item = &'a mut T;
    type IntoIter = IterMut<'a, T>;
    #[inline]
    fn into_iter(self) -> IterMut<'a, T> {
        self.data.iter_mut()
    }
}

/* -------------------- FromIterator -------------------- */

// Build a Vector from any iterator; panics if length != N.
impl<T, const N: usize> FromIterator<T> for Vector<T, N> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut it = iter.into_iter();
        let mut data: core::mem::MaybeUninit<[T; N]> = core::mem::MaybeUninit::uninit();
        let ptr = unsafe { &mut *data.as_mut_ptr() } as *mut [T; N] as *mut T;

        for i in 0..N {
            unsafe {
                ptr.add(i)
                    .write(it.next().expect("FromIterator: not enough elements"));
            }
        }
        if it.next().is_some() {
            panic!("FromIterator: too many elements");
        }
        Vector {
            data: unsafe { data.assume_init() },
        }
    }
}

/* -------------------- Serde -------------------- */

#[cfg(feature = "serde")]
use core::{fmt, marker::PhantomData};
#[cfg(feature = "serde")]
use serde::{
    de::{self, SeqAccess, Visitor},
    ser::{SerializeSeq,SerializeTuple},
    Deserialize, Deserializer, Serialize, Serializer,
};

#[cfg(feature = "serde")]
impl<T: Serialize, const N: usize> Serialize for Vector<T, N> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut tup = serializer.serialize_tuple(N)?;
        for elem in &self.data {
            tup.serialize_element(elem)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serde")]
impl<'de, T: Real + Deserialize<'de>, const N: usize> Deserialize<'de> for Vector<T, N> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct VectorVisitor<T, const N: usize> {
            marker: PhantomData<T>,
        }

        impl<'de, T: Real + Deserialize<'de>, const N: usize> Visitor<'de> for VectorVisitor<T, N> {
            type Value = Vector<T, N>;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                write!(formatter, "an array of length {}", N)
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let mut data = [T::ZERO; N];
                for i in 0..N {
                    data[i] = seq
                        .next_element()?
                        .ok_or_else(|| de::Error::invalid_length(i, &self))?;
                }
                Ok(Vector { data })
            }
        }

        deserializer.deserialize_tuple(
            N,
            VectorVisitor::<T, N> {
                marker: PhantomData,
            },
        )
    }
}

/* -------------------- Display (std) -------------------- */

#[cfg(feature = "std")]
impl<T: Real + std::fmt::Display, const N: usize> std::fmt::Display for Vector<T, N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[")?;
        for (i, val) in self.data.iter().enumerate() {
            write!(f, "{}", val)?;
            if i < N - 1 {
                write!(f, " ")?;
            }
        }
        write!(f, "]")
    }
}
