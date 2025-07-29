use num_traits::Float;
use crate::attitude::Quaternion;
use crate::coordinate::Cartesian;
use crate::math::Vector;

/// Marker trait for any frame (inertial or rotating)
pub trait ReferenceFrame {
    
}

/// Trait for fixed (non-rotating) frames
pub trait FixedFrame<T: Float> : ReferenceFrame {}

/// Trait for rotating frames with time-dependent orientation
pub trait RotatingFrame<T: Float, RF: ReferenceFrame> {
    /// Angular velocity in the body-fixed frame
    fn angular_velocity(&self) -> Cartesian<T, RF>;

    /// Epoch of the angular velocity
    fn epoch(&self) -> T;

    /// Orientation at time `t` relative to inertial parent
    fn orientation_at(&self, t: T) -> Quaternion<T>;
}

/// Trait for computing the rotation from `FROM` to `TO`
pub trait RotationBetween<T: Float, FROM: ReferenceFrame, TO: ReferenceFrame> {
    fn rotation(t: T) -> Quaternion<T>;
}