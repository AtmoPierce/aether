use core::ops::{Add, Sub, Mul, Div, Neg, AddAssign, SubAssign};
use num_traits::Float;
use crate::utils::{ToRadians, ToDegrees};
use core::slice::{Iter, IterMut};
use core::array::IntoIter as ArrayIntoIter;
use core::iter::FromIterator;

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
    pub fn normalize(&self) -> Self {
        self / self.norm()
    }
    pub fn try_normalize(&self) -> (Self, T) {
        let n = self.norm();
        if n > T::epsilon() {
            (*self / n, n)
        } else {
            (*self, n) // unchanged; caller must handle
        }
    }
    pub fn magnitude(&self) -> T {
        self.norm()
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

// Iterator
// --- AsRef / AsMut ---
impl<T, const N: usize> AsRef<[T; N]> for Vector<T, N> {
    #[inline] fn as_ref(&self) -> &[T; N] { &self.data }
}
impl<T, const N: usize> AsMut<[T; N]> for Vector<T, N> {
    #[inline] fn as_mut(&mut self) -> &mut [T; N] { &mut self.data }
}

// --- Convenience constructors ---
impl<T, const N: usize> From<[T; N]> for Vector<T, N> {
    #[inline] fn from(data: [T; N]) -> Self { Self { data } }
}

// --- Borrowed iteration ---
impl<T, const N: usize> Vector<T, N> {
    /// Iterate by reference (`&T`)
    #[inline] pub fn iter(&self) -> Iter<'_, T> { self.data.iter() }
    /// Iterate by mutable reference (`&mut T`)
    #[inline] pub fn iter_mut(&mut self) -> IterMut<'_, T> { self.data.iter_mut() }
    /// If you specifically want values (`T`), keep this helper (requires `T: Copy`)
    #[inline] pub fn iter_copied(&self) -> impl Iterator<Item = T> + '_
    where T: Copy { self.data.iter().copied() }
}

// --- IntoIterator (by value): moves out each T ---
impl<T, const N: usize> IntoIterator for Vector<T, N> {
    type Item = T;
    type IntoIter = ArrayIntoIter<T, N>;
    #[inline] fn into_iter(self) -> Self::IntoIter { self.data.into_iter() }
}

// --- IntoIterator for &Vector: yields &T ---
impl<'a, T, const N: usize> IntoIterator for &'a Vector<T, N> {
    type Item = &'a T;
    type IntoIter = Iter<'a, T>;
    #[inline] fn into_iter(self) -> Self::IntoIter { self.data.iter() }
}

// --- IntoIterator for &mut Vector: yields &mut T ---
impl<'a, T, const N: usize> IntoIterator for &'a mut Vector<T, N> {
    type Item = &'a mut T;
    type IntoIter = IterMut<'a, T>;
    #[inline] fn into_iter(self) -> Self::IntoIter { self.data.iter_mut() }
}
// Iterator Helpers
impl<T, const N: usize> Vector<T, N> {
    pub fn map<U, F: FnMut(T) -> U>(self, mut f: F) -> Vector<U, N> {
        let data = self.data.map(|x| f(x));
        Vector { data }
    }
}

// Build a Vector from any iterator; panics if length != N.
impl<T, const N: usize> FromIterator<T> for Vector<T, N> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        // Fill fixed-size array from iterator
        let mut it = iter.into_iter();
        let mut data: core::mem::MaybeUninit<[T; N]> = core::mem::MaybeUninit::uninit();
        let ptr = unsafe { &mut *data.as_mut_ptr() } as *mut [T; N] as *mut T;

        // Write elements
        for i in 0..N {
            unsafe {
                ptr.add(i).write(it.next().expect("FromIterator: not enough elements"));
            }
        }
        // Ensure no extra elements
        if it.next().is_some() {
            panic!("FromIterator: too many elements");
        }
        Vector { data: unsafe { data.assume_init() } }
    }
}



// Serialization with serde
#[cfg(feature = "serde")]
use serde::{Serialize, Deserialize, Serializer, Deserializer, ser::SerializeSeq, de::{self, SeqAccess, Visitor}};
#[cfg(feature = "serde")]
use core::{fmt, marker::PhantomData};
#[cfg(feature = "serde")]
impl<T: Serialize, const N: usize> Serialize for Vector<T, N> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(N))?;
        for elem in &self.data {
            seq.serialize_element(elem)?;
        }
        seq.end()
    }
}

#[cfg(feature = "serde")]
impl<'de, T: Float + Deserialize<'de>, const N: usize> Deserialize<'de> for Vector<T, N> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct VectorVisitor<T, const N: usize> {
            marker: PhantomData<T>,
        }

        impl<'de, T: Float + Deserialize<'de>, const N: usize> Visitor<'de> for VectorVisitor<T, N> {
            type Value = Vector<T, N>;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                write!(formatter, "an array of length {}", N)
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let mut data = [T::zero(); N];
                for i in 0..N {
                    data[i] = seq
                        .next_element()?
                        .ok_or_else(|| de::Error::invalid_length(i, &self))?;
                }
                Ok(Vector { data })
            }
        }

        deserializer.deserialize_tuple(N, VectorVisitor::<T, N> { marker: PhantomData })
    }
}


// std 
#[cfg(feature = "std")]
impl<T: Float + std::fmt::Display, const N: usize> std::fmt::Display for Vector<T, N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
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
