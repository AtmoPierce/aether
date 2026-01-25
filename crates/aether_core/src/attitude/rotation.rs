use crate::attitude::{DirectionCosineMatrix, Quaternion};
use crate::math::Vector;
use crate::reference_frame::ReferenceFrame;
use crate::real::Real;

use core::fmt;
use core::ops::Mul;

/// A body-fixed rotation: Inertial -> Body
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rotation<T: Real, Inertial: ReferenceFrame, Body: ReferenceFrame> {
    quat: Quaternion<T, Inertial, Body>,
}

impl<T: Real, Inertial: ReferenceFrame, Body: ReferenceFrame>
    Rotation<T, Inertial, Body>
{
    /// Create from a quaternion (assumed normalized).
    pub fn from_quaternion(q: Quaternion<T, Inertial, Body>) -> Self {
        Self { quat: q }
    }

    /// Access the underlying quaternion.
    pub fn quaternion(&self) -> &Quaternion<T, Inertial, Body> {
        &self.quat
    }

    /// Compose two rotations:
    /// (Inertial -> Body1) * (Body1 -> Body2) = (Inertial -> Body2)
    pub fn compose<Next: ReferenceFrame>(
        &self,
        rhs: &Rotation<T, Body, Next>,
    ) -> Rotation<T, Inertial, Next> {
        Rotation::from_quaternion(&self.quat * &rhs.quat)
    }
}

//
// ===== Conversions =====
//

impl<T: Real, Inertial: ReferenceFrame, Body: ReferenceFrame>
    TryFrom<&DirectionCosineMatrix<T, Inertial, Body>>
    for Rotation<T, Inertial, Body>
{
    type Error = ();

    fn try_from(
        dcm: &DirectionCosineMatrix<T, Inertial, Body>,
    ) -> Result<Self, Self::Error> {
        Quaternion::try_from(dcm).map(Self::from_quaternion)
    }
}

impl<T: Real, Inertial: ReferenceFrame, Body: ReferenceFrame>
    core::convert::From<&Quaternion<T, Inertial, Body>>
    for Rotation<T, Inertial, Body>
{
    fn from(q: &Quaternion<T, Inertial, Body>) -> Self {
        Self::from_quaternion(q.normalized())
    }
}

//
// ===== Kinematics =====
//

impl<T: Real, Inertial: ReferenceFrame, Body: ReferenceFrame>
    Rotation<T, Inertial, Body>
{
    /// Integrate body angular velocity (expressed in Body frame).
    ///
    /// q_dot = 0.5 * q * w
    pub fn integrate(&self, omega_body: Vector<T, 3>, dt: T) -> Self {
        let delta_q =
            Quaternion::<T, Body, Body>::from_angular_velocity(omega_body, dt);

        Self::from_quaternion(&self.quat * &delta_q)
    }
}

//
// ===== Quaternion helper =====
//

impl<T: Real, F: ReferenceFrame>
    Quaternion<T, F, F>
{
    pub fn from_angular_velocity(omega: Vector<T, 3>, dt: T) -> Self {
        let mag = omega.norm();
        if mag == T::ZERO {
            return Self::identity();
        }

        let axis = omega / mag;
        let angle = mag * dt;
        let half = angle * T::from_f32(0.5);
        let s = half.sin();

        Self::new(
            half.cos(),
            axis[0] * s,
            axis[1] * s,
            axis[2] * s,
        )
    }
}

//
// ===== Display =====
//

#[cfg(feature = "std")]
impl<T, Inertial, Body> fmt::Display
    for Rotation<T, Inertial, Body>
where
    T: Real + fmt::Display,
    Inertial: ReferenceFrame,
    Body: ReferenceFrame,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Rotation: {}", self.quat)
    }
}
