use crate::attitude::Quaternion;
use crate::coordinate::Cartesian;
use crate::math::Vector;
use crate::real::Real;
/// Marker trait for any frame (inertial or rotating)
pub trait ReferenceFrame {}

/// Trait for fixed (non-rotating) frames
pub trait FixedFrame<T: Real>: ReferenceFrame {}

/// Trait for rotating frames with time-dependent orientation
pub trait RotatingFrame<T: Real, F: ReferenceFrame>: ReferenceFrame {
    /// Angular velocity expressed in the frame `F`
    fn angular_velocity(&self) -> Cartesian<T, F>;

    /// Reference epoch
    fn epoch(&self) -> T;

    /// Orientation of this frame at time `t`
    ///
    /// This is a passive transform: F(epoch) -> F(t)
    fn orientation_at(&self, t: T) -> Quaternion<T, F, F>;
}

/// Trait for computing the rotation from `From` to `From`
pub trait RotationBetween<T: Real, From: ReferenceFrame, To: ReferenceFrame> {
    fn rotation(t: T) -> Quaternion<T, From, To>;
}
