use crate::reference_frame::{RotatingFrame, ReferenceFrame};
use crate::coordinate::Cartesian;
use crate::attitude::Quaternion;
use crate::math::Vector;
use num_traits::Float;
use core::marker::PhantomData;

#[derive(Debug, Default, Clone, Copy)]
pub struct ITRF<T: Float> {
    pub epoch: T,
}

impl<T: Float> ReferenceFrame for ITRF<T> {}

impl<T: Float> RotatingFrame<T, ITRF<T>> for ITRF<T> {
    fn angular_velocity(&self) -> Cartesian<T, ITRF<T>> {
        Cartesian {
            data: Vector::new([T::zero(), T::zero(), T::from(7.2921150e-5).unwrap()]),
            _reference_frame: PhantomData,
        }
    }

    fn epoch(&self) -> T {
        self.epoch
    }

    fn orientation_at(&self, t: T) -> Quaternion<T> {
        let omega = T::from(7.2921150e-5).unwrap(); // Earth's rotation rate [rad/s] (ITRF is actually more complicated and comes from ephermerides, more on that later...)
        let dtheta = omega * (t - self.epoch);
        let half_theta = dtheta * T::from(0.5).unwrap();

        Quaternion {
            data: Vector::new([
                half_theta.cos(),
                T::zero(),
                T::zero(),
                half_theta.sin(),
            ]),
        }
    }
}