use crate::attitude::Quaternion;
use crate::coordinate::Cartesian;
use crate::math::Vector;
use crate::reference_frame::{ReferenceFrame, RotatingFrame};
use crate::real::Real;

use core::marker::PhantomData;

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Default, Clone, Copy)]
pub struct ITRF<T: Real> {
    pub epoch: T,
}

impl<T: Real> ReferenceFrame for ITRF<T> {}

impl<T: Real> RotatingFrame<T, ITRF<T>> for ITRF<T> {
    fn angular_velocity(&self) -> Cartesian<T, ITRF<T>> {
        Cartesian {
            data: Vector::new([T::ZERO, T::ZERO, T::from_f64(7.2921150e-5)]),
            _reference_frame: PhantomData,
        }
    }

    fn epoch(&self) -> T {
        self.epoch
    }

    fn orientation_at(&self, t: T) -> Quaternion<T, ITRF<T>, ITRF<T>> {
        // Earth's rotation rate [rad/s]
        let omega = T::from_f64(7.2921150e-5);

        let dtheta = (t - self.epoch) * omega;
        let half_theta = dtheta * T::from_f32(0.5);

        // Passive frame transform: ITRF(epoch) -> ITRF(t)
        Quaternion::new(
            half_theta.cos(),
            T::ZERO,
            T::ZERO,
            half_theta.sin(),
        )
    }
}
