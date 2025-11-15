use core::marker::PhantomData;
use core::ops::{Add, Sub, Mul, Div, Neg};
use fpx::design::Real;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Quantity<Unit, T: Real> {
    pub(crate) value: T,
    pub(crate) _marker: PhantomData<Unit>,
}

impl<Unit, T: Real> Quantity<Unit, T> {
    #[inline]
    pub const fn new(value: T) -> Self {
        Self { value, _marker: PhantomData }
    }

    #[inline]
    pub const fn raw(self) -> T {
        self.value
    }
}

impl<Unit, T: Real> Add for Quantity<Unit, T> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.value + rhs.value)
    }
}

impl<Unit, T: Real> Sub for Quantity<Unit, T> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.value - rhs.value)
    }
}

impl<Unit, T: Real> Div<T> for Quantity<Unit, T> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Self::new(self.value / rhs)
    }
}

impl<Unit, T: Real> Mul<T> for Quantity<Unit, T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Self::new(self.value * rhs)
    }
}

impl<Unit, T: Real> Neg for Quantity<Unit, T> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.value)
    }
}

impl<Unit, T: Real> From<T> for Quantity<Unit, T> {
    #[inline]
    fn from(value: T) -> Self {
        Self::new(value)
    }
}