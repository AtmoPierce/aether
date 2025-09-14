use super::cylindrical::Cylindrical;
use super::spherical::Spherical;
use crate::attitude::Quaternion;
use crate::math::{Vector, Matrix};
use crate::reference_frame::ReferenceFrame;

use num_traits::Float;
use core::marker::PhantomData; // Reference frame tracking.
use core::slice::{Iter, IterMut};

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Cartesian<T: Float, ReferenceFrame> {
    pub data: Vector<T, 3>,
    pub _reference_frame: PhantomData<ReferenceFrame>
}

impl<T: Float, RF: ReferenceFrame> Cartesian<T, RF>{
    pub fn new(x: T, y: T, z: T) -> Self {
        Self { data: Vector { data: [x, y, z] }, _reference_frame: PhantomData }
    }
    pub fn x(&self) -> T{ self.data.data[0] }
    pub fn y(&self) -> T{ self.data.data[1] }
    pub fn z(&self) -> T{ self.data.data[2] }

    pub fn transform<To: ReferenceFrame>(&self, q_from_to: Quaternion<T>) -> Cartesian<T, To> {
        Cartesian {
            data: q_from_to.rotate_vector(self.data),
            _reference_frame: PhantomData::<To>,
        }
    }
}


// If you want a cross product method, add it to the main impl block:
impl<T: Float, RF: ReferenceFrame> Cartesian<T, RF> {
    pub fn cross(&self, rhs: &Self) -> Self {
        Self {
            data: self.data.cross(&rhs.data),
            _reference_frame: PhantomData,
        }
    }
    pub fn norm(&self)->T{
        self.data.norm()
    }
    pub fn normalize(&self) -> Self {
        Self {
            data: self.data.normalize(),
            _reference_frame: PhantomData,
        }
    }
    pub fn dot(&self, rhs: &Self) -> T {
        self.data.dot(&rhs.data)
    }
    pub fn magnitude(&self)->T{
        self.data.norm()
    }
}

use core::ops::{Add, AddAssign, Sub, SubAssign, Neg, Mul, Div};
// ----- Add -----
impl<T: Float, RF> Add for Cartesian<T, RF> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            data: self.data + rhs.data,
            _reference_frame: PhantomData,
        }
    }
}
impl<'a, T: Float, RF> Add for &'a Cartesian<T, RF>{
    type Output = Cartesian<T, RF>;
    fn add(self, rhs: Self) -> Self::Output {
        Cartesian {
            data: self.data + rhs.data,
            _reference_frame: PhantomData,
        }
    }
}

impl<'a, T: Float, RF> Add<Cartesian<T, RF>> for &'a Cartesian<T, RF>{
    type Output = Cartesian<T, RF>;
    fn add(self, rhs: Cartesian<T, RF>) -> Self::Output {
        Cartesian{
            data: self.data + rhs.data,
            _reference_frame: PhantomData,
        }
    }
}
impl<'a, T: Float, RF> Add<&'a Cartesian<T, RF>> for Cartesian<T, RF>
{
    type Output = Cartesian<T, RF>;
    fn add(self, rhs: &'a Cartesian<T, RF>) -> Self::Output {
        Cartesian{
            data: self.data + rhs.data,
            _reference_frame: PhantomData,
        }
    }
}

// ----- AddAssign -----
impl<T: Float + AddAssign, RF> AddAssign for Cartesian<T, RF> {
    fn add_assign(&mut self, rhs: Self) {
        self.data += rhs.data;
    }
}

// ----- Sub -----
impl<T: Float, RF> Sub for Cartesian<T, RF> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            data: self.data - rhs.data,
            _reference_frame: PhantomData,
        }
    }
}
impl<'a, T: Float, RF> Sub for &'a Cartesian<T, RF>{
    type Output = Cartesian<T, RF>;
    fn sub(self, rhs: Self) -> Self::Output {
        Cartesian{
            data: self.data - rhs.data,
            _reference_frame: PhantomData,
        }    
    }
}
impl<'a, T: Float, RF> Sub<Cartesian<T, RF>> for &'a Cartesian<T, RF>{
    type Output = Cartesian<T, RF>;
    fn sub(self, rhs: Cartesian<T, RF>) -> Self::Output {
        Cartesian{
            data: self.data - rhs.data,
            _reference_frame: PhantomData,
        }
    }
}
impl<'a, T: Float, RF> Sub<&'a Cartesian<T, RF>> for Cartesian<T, RF>
where
    Cartesian<T, RF>: Clone,
{
    type Output = Cartesian<T, RF>;
    fn sub(self, rhs: &'a Cartesian<T, RF>) -> Self::Output {
        Cartesian{
            data: self.data - rhs.data,
            _reference_frame: PhantomData,
        }
    }
}

// ----- SubAssign -----
impl<T: Float + SubAssign, RF> SubAssign for Cartesian<T, RF> {
    fn sub_assign(&mut self, rhs: Self) {
        self.data -= rhs.data;
    }
}

// ----- Neg -----
impl<T: Float, RF> Neg for Cartesian<T, RF> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Cartesian {
            data: -self.data,
            _reference_frame: PhantomData,
        }
    }
}
impl<'a, T: Float, RF> Neg for &'a Cartesian<T, RF>
where
    Cartesian<T, RF>: Clone,
{
    type Output = Cartesian<T, RF>;
    fn neg(self) -> Self::Output {
        Cartesian {
            data: -self.data,
            _reference_frame: PhantomData,
        }    
    }
}

// ----- Scalar Mul -----
impl<T: Float + Copy, RF> Mul<T> for Cartesian<T, RF> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Self {
            data: self.data * rhs,
            _reference_frame: PhantomData,
        }
    }
}
impl<'a, T: Float + Copy, RF> Mul<T> for &'a Cartesian<T, RF>{
    type Output = Cartesian<T, RF>;
    fn mul(self, rhs: T) -> Self::Output {
        Cartesian {
            data: self.data * rhs,
            _reference_frame: PhantomData,
        }
    }
}

// ----- Scalar Div -----
impl<T: Float + Copy, RF> Div<T> for Cartesian<T, RF> {
    type Output = Self;
    fn div(self, rhs: T) -> Self::Output {
        Cartesian {
            data: self.data / rhs,
            _reference_frame: PhantomData,
        }
    }
}
// Indexing
use core::ops::{Index, IndexMut};

impl<T: Float, ReferenceFrame> Index<usize> for Cartesian<T, ReferenceFrame> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data.data[index]
    }
}

impl<T: Float, ReferenceFrame> IndexMut<usize> for Cartesian<T, ReferenceFrame> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data.data[index]
    }
}

impl<'a, T: Float + Copy, RF> Div<T> for &'a Cartesian<T, RF>
where
    Cartesian<T, RF>: Clone,
{
    type Output = Cartesian<T, RF>;
    fn div(self, rhs: T) -> Self::Output {
        Cartesian{
            data: self.data / rhs,
            _reference_frame: PhantomData,
        }
    }
}


impl<T: Float + Copy + Default, F> Mul<Cartesian<T, F>> for Matrix<T, 3, 3> {
    type Output = Cartesian<T, F>;
    fn mul(self, rhs: Cartesian<T, F>) -> Self::Output {
        let v = self * rhs.data; // Use existing `Matrix * Vector` implementation
        Cartesian {
            data: v,
            _reference_frame: PhantomData,
        }
    }
}


// Conversions
impl<T: Float, RF: ReferenceFrame> From<&Spherical<T>> for Cartesian<T, RF> {
    fn from(p: &Spherical<T>) -> Self {
        let r = p.r();
        let phi = p.phi();     // azimuth
        let theta = p.theta(); // inclination
        let sin_theta = theta.sin();
        let x = r * sin_theta * phi.cos();
        let y = r * sin_theta * phi.sin();
        let z = r * theta.cos();
        Cartesian::new(x, y, z)
    }
}

impl<T: Float, RF: ReferenceFrame> From<&Cylindrical<T>> for Cartesian<T, RF>{
    fn from(c: &Cylindrical<T>) -> Self {
        let x = c.r()*c.theta().cos();
        let y = c.r()*c.theta().sin();
        let z = c.z(); 
        Cartesian::new(x, y, z)
    }
}

// Iterator
// --- AsRef / AsMut ---
impl<T: num_traits::Float, RF> AsRef<[T; 3]> for Cartesian<T, RF> {
    #[inline] fn as_ref(&self) -> &[T; 3] { &self.data.data }
}
impl<T: num_traits::Float, RF> AsMut<[T; 3]> for Cartesian<T, RF> {
    #[inline] fn as_mut(&mut self) -> &mut [T; 3] { &mut self.data.data }
}

// --- Convenience methods ---
impl<T: num_traits::Float, RF> Cartesian<T, RF> {
    #[inline] pub fn iter(&self) -> Iter<'_, T> { self.data.iter() }
    #[inline] pub fn iter_mut(&mut self) -> IterMut<'_, T> { self.data.iter_mut() }
}

// --- IntoIterator (by value) ---
impl<T: num_traits::Float, RF> IntoIterator for Cartesian<T, RF> {
    type Item = T;
    type IntoIter = <Vector<T, 3> as IntoIterator>::IntoIter;
    #[inline] fn into_iter(self) -> Self::IntoIter { self.data.into_iter() }
}

// --- IntoIterator for &Cartesian ---
impl<'a, T: num_traits::Float, RF> IntoIterator for &'a Cartesian<T, RF> {
    type Item = &'a T;
    type IntoIter = <&'a Vector<T, 3> as IntoIterator>::IntoIter;
    #[inline] fn into_iter(self) -> Self::IntoIter { (&self.data).into_iter() }
}

// --- IntoIterator for &mut Cartesian ---
impl<'a, T: num_traits::Float, RF> IntoIterator for &'a mut Cartesian<T, RF> {
    type Item = &'a mut T;
    type IntoIter = <&'a mut Vector<T, 3> as IntoIterator>::IntoIter;
    #[inline] fn into_iter(self) -> Self::IntoIter { (&mut self.data).into_iter() }
}
// Iterator Helpers
impl<T: num_traits::Float, RF> Cartesian<T, RF> {
    /// Map each component to a new type, preserving the reference frame.
    pub fn map<U: num_traits::Float, F>(self, f: F) -> Cartesian<U, RF>
    where
        F: FnMut(T) -> U,
    {
        Cartesian {
            data: self.data.map(f),  // uses your Vector<T,3>::map(...)
            _reference_frame: PhantomData,
        }
    }
}

// Serialization with Serde
#[cfg(feature = "serde")]
use serde::{Serialize, Deserialize};
#[cfg(feature = "serde")]
impl<T: Float + Serialize, RF: ReferenceFrame> Serialize for Cartesian<T, RF> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.data.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T: Float + Deserialize<'de>, RF: ReferenceFrame> Deserialize<'de> for Cartesian<T, RF> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let data = Vector::<T, 3>::deserialize(deserializer)?;
        Ok(Cartesian {
            data,
            _reference_frame: PhantomData,
        })
    }
}

// Print/Display
use core::{fmt, num};
impl<T, RF> fmt::Display for Cartesian<T, RF>
where
    T: Float + fmt::Display,
    RF: ReferenceFrame,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Cartesian [x,y,z]: [{}, {}, {}]",
            self.data.data[0],
            self.data.data[1],
            self.data.data[2]
        )
    }
}