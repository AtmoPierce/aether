use core::ops::{Add, Sub, Mul, Div, Neg, AddAssign, SubAssign};
use num_traits::Float;
use crate::utils::{ToRadians, ToDegrees};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector<T, const N: usize> {
    pub data: [T; N],
}

impl<T: Default + Copy, const N: usize> Default for Vector<T, N> {
    fn default() -> Self {
        Self {
            data: [T::default(); N],
        }
    }
}

impl<T: Copy, const N: usize> Vector<T, N> {
    pub fn new(data: [T; N]) -> Self {
        Self { data }
    }
}

// ===== Add =====

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

// ===== Sub =====

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

// ===== Neg =====

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

// Scalar multiplication: Vector * scalar
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

// Scalar division: Vector / scalar
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

// By reference:
impl<'a, T: Float + Copy, const N: usize> Mul<T> for &'a Vector<T, N> {
    type Output = Vector<T, N>;
    fn mul(self, rhs: T) -> Self::Output {
        (*self).clone() * rhs
    }
}

impl<'a, T: Float + Copy, const N: usize> Div<T> for &'a Vector<T, N> {
    type Output = Vector<T, N>;
    fn div(self, rhs: T) -> Self::Output {
        (*self).clone() / rhs
    }
}

impl<T, const N: usize> Vector<T, N>
where
    T: Mul<Output = T> + Copy + Float,
{
    pub fn dot(&self, rhs: &Self) -> T {
        let mut acc = T::zero();
        for i in 0..N {
            acc = acc + self.data[i] * rhs.data[i];
        }
        acc
    }
    pub fn norm(&self) -> T {
        let mut acc = T::zero();
        for &x in &self.data {
            acc = acc + x * x;
        }
        acc.sqrt()
    }
    pub fn angle(&self, rhs: &Self) -> T {
        let dot = self.dot(rhs);
        let denom = self.norm() * rhs.norm();
        if denom == T::zero() {
            T::zero()
        } else {
            let cos_theta = (dot / denom).max(T::from(-1.0).unwrap()).min(T::one());
            cos_theta.acos()
        }
    }
}

// Cross-product: only for N=3
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

// units
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

// Behavior
impl<T: Copy, const N: usize> Vector<T, N> {
    pub fn iter(&self) -> impl Iterator<Item = T> + '_ {
        self.data.iter().copied()
    }
}
use core::ops::{Index, IndexMut};
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

// std 
#[cfg(feature = "std")]
use std::fmt;

#[cfg(feature = "std")]
impl<T: Float + fmt::Display, const N: usize> fmt::Display for Vector<T, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[")?;
        for (i, val) in self.data.iter().enumerate() {
            write!(f, "{}", val)?;
            if i < N - 1 {
                write!(f, " ")?; // space between elements
            }
        }
        write!(f, "]")
    }
}
