use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::real::Real;

pub trait ComplexField:
    Copy
    + Clone
    + core::fmt::Debug
    + PartialEq
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self>
    + MulAssign
    + Mul<Self::RealPart, Output = Self>
    + MulAssign<Self::RealPart>
    + Div<Self::RealPart, Output = Self>
    + DivAssign<Self::RealPart>
    + Neg<Output = Self>
{
    type RealPart: Real;

    fn new(re: Self::RealPart, im: Self::RealPart) -> Self;
    fn zero() -> Self;
    fn from_real(re: Self::RealPart) -> Self;
    fn real(self) -> Self::RealPart;
    fn imag(self) -> Self::RealPart;
    fn conj(self) -> Self;
    fn norm_sqr(self) -> Self::RealPart;
    fn norm(self) -> Self::RealPart;
    fn from_polar(radius: Self::RealPart, angle: Self::RealPart) -> Self;
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex<T: Real> {
    pub re: T,
    pub im: T,
}

impl<T: Real> Complex<T> {
    pub const fn new(re: T, im: T) -> Self {
        Self { re, im }
    }

    pub fn zero() -> Self {
        Self::new(T::ZERO, T::ZERO)
    }

    pub fn from_real(re: T) -> Self {
        Self::new(re, T::ZERO)
    }

    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    pub fn norm_sqr(self) -> T {
        self.re * self.re + self.im * self.im
    }

    pub fn norm(self) -> T {
        self.norm_sqr().sqrt()
    }

    pub fn from_polar(radius: T, angle: T) -> Self {
        Self::new(radius * angle.cos(), radius * angle.sin())
    }
}

impl<T: Real> ComplexField for Complex<T> {
    type RealPart = T;

    fn new(re: Self::RealPart, im: Self::RealPart) -> Self {
        Self { re, im }
    }

    fn zero() -> Self {
        Self::new(T::ZERO, T::ZERO)
    }

    fn from_real(re: Self::RealPart) -> Self {
        Self::new(re, T::ZERO)
    }

    fn real(self) -> Self::RealPart {
        self.re
    }

    fn imag(self) -> Self::RealPart {
        self.im
    }

    fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    fn norm_sqr(self) -> Self::RealPart {
        self.re * self.re + self.im * self.im
    }

    fn norm(self) -> Self::RealPart {
        self.norm_sqr().sqrt()
    }

    fn from_polar(radius: Self::RealPart, angle: Self::RealPart) -> Self {
        Self::new(radius * angle.cos(), radius * angle.sin())
    }
}

impl<T: Real> Add for Complex<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl<T: Real> AddAssign for Complex<T> {
    fn add_assign(&mut self, rhs: Self) {
        self.re = self.re + rhs.re;
        self.im = self.im + rhs.im;
    }
}

impl<T: Real> Sub for Complex<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl<T: Real> SubAssign for Complex<T> {
    fn sub_assign(&mut self, rhs: Self) {
        self.re = self.re - rhs.re;
        self.im = self.im - rhs.im;
    }
}

impl<T: Real> Neg for Complex<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.re, -self.im)
    }
}

impl<T: Real> Mul for Complex<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl<T: Real> MulAssign for Complex<T> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<T: Real> Mul<T> for Complex<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Self::new(self.re * rhs, self.im * rhs)
    }
}

impl<T: Real> MulAssign<T> for Complex<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.re = self.re * rhs;
        self.im = self.im * rhs;
    }
}

impl<T: Real> Div<T> for Complex<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Self::new(self.re / rhs, self.im / rhs)
    }
}

impl<T: Real> DivAssign<T> for Complex<T> {
    fn div_assign(&mut self, rhs: T) {
        self.re = self.re / rhs;
        self.im = self.im / rhs;
    }
}