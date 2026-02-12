use crate::coordinate::Cartesian;
use crate::attitude::Euler;
use crate::math::Vector;
use core::marker::PhantomData;
use crate::real::Real;
use crate::reference_frame::ReferenceFrame;
pub trait TrapezoidalIntegrate {
    type Scalar;

    /// Integrates self and last using the trapezoidal rule and time step dt.
    fn integrate_trapezoidal(&self, last: &Self, dt: Self::Scalar) -> Self;
}

impl<T: Real + Copy, const N: usize> TrapezoidalIntegrate for Vector<T, N> {
    type Scalar = T;
    fn integrate_trapezoidal(&self, last: &Self, dt: T) -> Self {
        (*self + *last) * T::from_f32(0.5) * dt
    }
}

impl<T: Real + Copy, RF> TrapezoidalIntegrate for Cartesian<T, RF> {
    type Scalar = T;
    fn integrate_trapezoidal(&self, last: &Self, dt: T) -> Self {
        Cartesian {
            data: (self.data + last.data) * T::from_f32(0.5) * dt,
            _reference_frame: PhantomData,
        }
    }
}

impl<T: Real + Copy, From: ReferenceFrame, To: ReferenceFrame> TrapezoidalIntegrate
    for Euler<T, From, To>
{
    type Scalar = T;

    fn integrate_trapezoidal(&self, last: &Self, dt: T) -> Self {
        let data = (self.data + last.data) * T::from_f32(0.5) * dt;
        Euler::new(data[0], data[1], data[2])
    }
}
