use crate::attitude::{DirectionCosineMatrix, Quaternion};
use crate::math::Vector;
use crate::reference_frame::ReferenceFrame;
use core::fmt;
use core::ops::Mul;
use crate::real::Real;
/// A body-fixed rotation, representing orientation of body w.r.t inertial frame.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Rotation<T: Real> {
    quat: Quaternion<T>,
}

impl<T: Real> Rotation<T> {
    /// Creates a new `Rotation` from a quaternion (assumes normalized).
    pub fn from_quaternion(q: Quaternion<T>) -> Self {
        Self { quat: q }
    }

    /// Returns the internal quaternion.
    pub fn quaternion(&self) -> Quaternion<T> {
        self.quat
    }

    /// Compose two rotations: self followed by rhs.
    pub fn compose(&self, rhs: &Self) -> Self {
        Self::from_quaternion(self.quat * rhs.quat)
    }
}

impl<T: Real, A: ReferenceFrame, B: ReferenceFrame> TryFrom<&DirectionCosineMatrix<T, A, B>>
    for Rotation<T>
{
    type Error = ();

    fn try_from(dcm: &DirectionCosineMatrix<T, A, B>) -> Result<Self, Self::Error> {
        Quaternion::try_from(dcm).map(Self::from_quaternion)
    }
}

impl<T: Real> From<&Quaternion<T>> for Rotation<T> {
    fn from(q: &Quaternion<T>) -> Self {
        Self::from_quaternion(*q)
    }
}

impl<T: Real> Mul for Rotation<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        self.compose(&rhs)
    }
}

impl<T: Real> Rotation<T> {
    // !TODO Need to validate this with tests and such
    pub fn integrate(&self, omega_b: Vector<T, 3>, dt: T) -> Self {
        let delta_q = Quaternion::from_angular_velocity(omega_b, dt);
        Self::from_quaternion(delta_q * self.quat)
    }
}

impl<T: Real> Quaternion<T> {
    pub fn from_angular_velocity(omega: Vector<T, 3>, dt: T) -> Self {
        // !TODO Need to validate this with tests and such
        let mag = omega.norm();
        if mag == T::ZERO {
            return Self::identity();
        }
        let axis = omega / mag;
        let angle = mag * dt;
        let half_angle = angle * T::from_f32(0.5);
        let sin_half = half_angle.sin();
        Self::new(
            half_angle.cos(),
            axis[0] * sin_half,
            axis[1] * sin_half,
            axis[2] * sin_half,
        )
    }
}

#[cfg(feature = "std")]
impl<T: Real + fmt::Display> fmt::Display for Rotation<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Rotation: {}", self.quat)
    }
}
