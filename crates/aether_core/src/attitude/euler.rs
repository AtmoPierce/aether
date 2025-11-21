use crate::attitude::{DirectionCosineMatrix, Quaternion};
use crate::math::Vector;
use crate::reference_frame::ReferenceFrame;
use crate::utils::angle_conversion::{ToDegrees, ToRadians};
use core::ops::{Add, Div, Mul, Neg, Sub};
use crate::real::Real;
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Euler<T: Real> {
    pub data: Vector<T, 3>, // [roll, pitch, yaw]
}

impl<T: Real> Euler<T> {
    pub fn new(roll: T, pitch: T, yaw: T) -> Self {
        Self {
            data: Vector {
                data: [roll, pitch, yaw],
            },
        }
    }

    pub fn roll(&self) -> T {
        self.data.data[0]
    }
    pub fn pitch(&self) -> T {
        self.data.data[1]
    }
    pub fn yaw(&self) -> T {
        self.data.data[2]
    }
}

impl<T: Real> From<&Quaternion<T>> for Euler<T> {
    fn from(q: &Quaternion<T>) -> Self {
        let w = q.data.data[0];
        let x = q.data.data[1];
        let y = q.data.data[2];
        let z = q.data.data[3];
        let two = T::ONE + T::ONE;
        let one = T::ONE;

        let sinr_cosp = two * (w * x + y * z);
        let cosr_cosp = one - two * (x * x + y * y);
        let roll = sinr_cosp.atan2(cosr_cosp);

        let sinp = two * (w * y - z * x);
        let half_pi = T::FRAC_PI_2;
        let pitch = if sinp.abs() >= one {
            half_pi * sinp.signum()
        } else {
            sinp.asin()
        };

        let siny_cosp = two * (w * z + x * y);
        let cosy_cosp = one - two * (y * y + z * z);
        let yaw = siny_cosp.atan2(cosy_cosp);

        Self::new(roll, pitch, yaw)
    }
}

impl<T: Real, A: ReferenceFrame, B: ReferenceFrame> From<&DirectionCosineMatrix<T, A, B>>
    for Euler<T>
{
    fn from(dcm: &DirectionCosineMatrix<T, A, B>) -> Self {
        let m = &dcm.as_matrix().data;
        let yaw = (m[0][1] / m[0][0]).atan();
        let pitch = -(m[0][2]).asin();
        let roll = (m[1][2] / m[2][2]).atan();
        Self::new(roll, pitch, yaw)
    }
}

// Add
impl<T: Real> Add for Euler<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            data: self.data + rhs.data,
        }
    }
}

// Sub
impl<T: Real> Sub for Euler<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            data: self.data - rhs.data,
        }
    }
}

// Neg
impl<T: Real> Neg for Euler<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self { data: -self.data }
    }
}

// Scalar multiplication
impl<T: Real> Mul<T> for Euler<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Self {
            data: self.data * rhs,
        }
    }
}

// Scalar division
impl<T: Real> Div<T> for Euler<T> {
    type Output = Self;
    fn div(self, rhs: T) -> Self {
        Self {
            data: self.data / rhs,
        }
    }
}

// Unit Conversions
impl<T> Euler<T>
where
    T: Real + ToRadians<Output = T> + ToDegrees<Output = T> + Copy,
{
    /// Create Euler angles from degrees vector.
    pub fn from_degrees_vec(deg: Vector<T, 3>) -> Self {
        Self {
            data: deg.to_radians(),
        }
    }

    /// Get Euler angles as degrees vector.
    pub fn to_degrees_vec(&self) -> Vector<T, 3> {
        self.data.to_degrees()
    }
}

// std
#[cfg(feature = "std")]
use std::fmt;

#[cfg(feature = "std")]
impl<T: Real + fmt::Display> fmt::Display for Euler<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.data)
    }
}

// Tests
#[cfg(test)]
#[path = "tests/euler_tests.rs"]
mod euler_tests;
