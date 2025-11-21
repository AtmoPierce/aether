use crate::coordinate::Cartesian;
use crate::math::Vector;
use core::marker::PhantomData;
use crate::real::Real;
pub trait Rk4Integrate<T: Real> {
    /// RK4 step given `self` as the current state,
    /// `f` as the derivative function, and `dt` as the time step.
    fn integrate_rk4<F>(&self, f: F, dt: T) -> Self
    where
        F: Fn(&Self) -> Self;
}

impl<T: Real + Copy, const N: usize> Rk4Integrate<T> for Vector<T, N> {
    fn integrate_rk4<F>(&self, f: F, dt: T) -> Self
    where
        F: Fn(&Self) -> Self,
    {
        let k1 = f(self);
        let k2 = f(&(*self + &k1 * (dt / T::from_f32(2.0))));
        let k3 = f(&(*self + &k2 * (dt / T::from_f32(2.0))));
        let k4 = f(&(*self + &k3 * dt));
        *self + (k1 + k2 + k2 + k3 + k3 + k4) * (dt / T::from_f32(6.0))
    }
}
impl<T: Real + Copy, RF> Rk4Integrate<T> for Cartesian<T, RF> {
    fn integrate_rk4<F>(&self, f: F, dt: T) -> Self
    where
        F: Fn(&Self) -> Self,
    {
        let half_dt = dt / T::from_f32(2.0);
        let sixth_dt = dt / T::from_f32(6.0);

        let k1 = f(self);
        let k2 = f(&(self + &k1 * half_dt));
        let k3 = f(&(self + &k2 * half_dt));
        let k4 = f(&(self + &k3 * dt));

        Self {
            data: self.data
                + (k1.data
                    + k2.data * T::from_f32(2.0)
                    + k3.data * T::from_f32(2.0)
                    + k4.data)
                    * sixth_dt,
            _reference_frame: PhantomData,
        }
    }
}
