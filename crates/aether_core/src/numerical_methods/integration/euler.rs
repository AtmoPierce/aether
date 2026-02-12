use crate::coordinate::Cartesian;
use crate::attitude::Euler;
use crate::math::Vector;
use crate::real::Real;
use crate::reference_frame::ReferenceFrame;
pub trait EulerIntegrate<T: Real> {
    fn integrate_euler(&self, dt: T) -> Self;
}

impl<T: Real + Copy, const N: usize> EulerIntegrate<T> for Vector<T, N> {
    fn integrate_euler(&self, dt: T) -> Self {
        *self * dt
    }
}

impl<T: Real + Copy, RF> EulerIntegrate<T> for Cartesian<T, RF> {
    fn integrate_euler(&self, dt: T) -> Self {
        Cartesian {
            data: self.data * dt,
            ..*self
        }
    }
}

impl<T: Real + Copy, From: ReferenceFrame, To: ReferenceFrame> EulerIntegrate<T>
    for Euler<T, From, To>
{
    fn integrate_euler(&self, dt: T) -> Self {
        let data = self.data * dt;
        Euler::new(data[0], data[1], data[2])
    }
}
