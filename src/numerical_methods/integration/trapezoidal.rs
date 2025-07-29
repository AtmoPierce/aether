use num_traits::Float;
use crate::math::Vector;
use crate::coordinate::Cartesian;
pub trait TrapezoidalIntegrate {
    type Scalar;

    /// Integrates self and last using the trapezoidal rule and time step dt.
    fn integrate_trapezoidal(&self, last: &Self, dt: Self::Scalar) -> Self;
}

impl<T: Float + Copy + Clone, const N: usize> TrapezoidalIntegrate for Vector<T, N> {
    type Scalar = T;
    fn integrate_trapezoidal(&self, last: &Self, dt: T) -> Self {
        (*self + *last) * T::from(0.5).unwrap() * dt
    }
}

impl<T: Float + Copy + Clone, RF: Clone> TrapezoidalIntegrate for Cartesian<T, RF> {
    type Scalar = T;
    fn integrate_trapezoidal(&self, last: &Self, dt: T) -> Self {
        Cartesian {
            data: (self.data + last.data) * (T::from(0.5).unwrap()) * dt,
            ..self.clone()
        }
    }
}