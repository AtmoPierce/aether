use crate::coordinate::Cartesian;
use crate::math::Vector;
use crate::real::Real;
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
