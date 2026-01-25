use crate::attitude::{DirectionCosineMatrix, Quaternion};
use crate::math::Vector;
use crate::reference_frame::ReferenceFrame;
use crate::utils::angle_conversion::{ToDegrees, ToRadians};
use crate::real::Real;

use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Euler<T: Real, From: ReferenceFrame, To: ReferenceFrame> {
    pub data: Vector<T, 3>, // [roll, pitch, yaw]
    _from: PhantomData<From>,
    _to: PhantomData<To>,
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Euler<T, From, To> {
    pub fn new(roll: T, pitch: T, yaw: T) -> Self {
        Self {
            data: Vector::new([roll, pitch, yaw]),
            _from: PhantomData,
            _to: PhantomData,
        }
    }

    #[inline] pub fn roll(&self) -> T { self.data[0] }
    #[inline] pub fn pitch(&self) -> T { self.data[1] }
    #[inline] pub fn yaw(&self) -> T { self.data[2] }
}

//
// ===== Conversions =====
//

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame>
    core::convert::From<&Quaternion<T, From, To>> for Euler<T, From, To>
{
    fn from(q: &Quaternion<T, From, To>) -> Self {
        let w = q.w();
        let x = q.i();
        let y = q.j();
        let z = q.k();

        let two = T::ONE + T::ONE;
        let one = T::ONE;

        // Roll (x-axis)
        let sinr_cosp = two * (w * x + y * z);
        let cosr_cosp = one - two * (x * x + y * y);
        let roll = sinr_cosp.atan2(cosr_cosp);

        // Pitch (y-axis)
        let sinp = two * (w * y - z * x);
        let half_pi = T::FRAC_PI_2;
        let pitch = if sinp.abs() >= one {
            half_pi * sinp.signum()
        } else {
            sinp.asin()
        };

        // Yaw (z-axis)
        let siny_cosp = two * (w * z + x * y);
        let cosy_cosp = one - two * (y * y + z * z);
        let yaw = siny_cosp.atan2(cosy_cosp);

        Self::new(roll, pitch, yaw)
    }
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame>
    core::convert::From<&DirectionCosineMatrix<T, From, To>>
    for Euler<T, From, To>
{
    fn from(dcm: &DirectionCosineMatrix<T, From, To>) -> Self {
        let m = &dcm.as_matrix().data;

        // ZYX convention (yaw-pitch-roll), passive From -> To
        let yaw   = m[0][1].atan2(m[0][0]);
        let pitch = -(m[0][2]).asin();
        let roll  = m[1][2].atan2(m[2][2]);

        Self::new(roll, pitch, yaw)
    }
}

//
// ===== Arithmetic =====
//

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Add for Euler<T, From, To> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self { data: self.data + rhs.data, _from: PhantomData, _to: PhantomData }
    }
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Sub for Euler<T, From, To> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self { data: self.data - rhs.data, _from: PhantomData, _to: PhantomData }
    }
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Neg for Euler<T, From, To> {
    type Output = Self;
    fn neg(self) -> Self {
        Self { data: -self.data, _from: PhantomData, _to: PhantomData }
    }
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Mul<T>
    for Euler<T, From, To>
{
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Self { data: self.data * rhs, _from: PhantomData, _to: PhantomData }
    }
}

impl<T: Real, From: ReferenceFrame, To: ReferenceFrame> Div<T>
    for Euler<T, From, To>
{
    type Output = Self;
    fn div(self, rhs: T) -> Self {
        Self { data: self.data / rhs, _from: PhantomData, _to: PhantomData }
    }
}

//
// ===== Unit conversions =====
//

impl<T, From, To> Euler<T, From, To>
where
    T: Real + ToRadians<Output = T> + ToDegrees<Output = T> + Copy,
    From: ReferenceFrame,
    To: ReferenceFrame,
{
    pub fn from_degrees_vec(deg: Vector<T, 3>) -> Self {
        Self {
            data: deg.to_radians(),
            _from: PhantomData,
            _to: PhantomData,
        }
    }

    pub fn to_degrees_vec(&self) -> Vector<T, 3> {
        self.data.to_degrees()
    }
}

// std
#[cfg(feature = "std")]
use std::fmt;

#[cfg(feature = "std")]
impl<T, From, To> fmt::Display for Euler<T, From, To>
where
    T: Real + fmt::Display,
    From: ReferenceFrame,
    To: ReferenceFrame,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.data)
    }
}

// Tests
#[cfg(test)]
#[path = "tests/euler_tests.rs"]
mod euler_tests;
