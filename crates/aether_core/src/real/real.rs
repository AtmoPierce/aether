use core::fmt::{Debug, Display, Formatter, Result as FmtResult};
use core::ops::{Add, Sub, Mul, Div, Neg};

pub trait Real:
    RealCast
    + Debug
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    const ZERO:    Self;
    const ONE:     Self;
    const PI:      Self;
    const FRAC_PI_2: Self;
    const EPSILON: Self;
    const INFINITY: Self;
    const NEG_INFINITY: Self;

    // Basic math operations (re-exports)
    fn abs(self) -> Self;
    fn signum(self) -> Self;

    fn floor(self) -> Self;
    fn ceil(self) -> Self;
    fn round(self) -> Self;
    fn trunc(self) -> Self;
    fn fract(self) -> Self;

    fn min(self, other: Self) -> Self;
    fn max(self, other: Self) -> Self;
    fn copysign(self, sign: Self) -> Self;

    fn sqrt(self) -> Self;

    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn tan(self) -> Self;

    fn asin(self) -> Self;
    fn acos(self) -> Self;
    fn atan(self) -> Self;
    fn atan2(self, other: Self) -> Self;

    fn exp(self) -> Self;
    fn exp2(self) -> Self;
    fn ln(self) -> Self;
    fn log2(self) -> Self;
    fn log10(self) -> Self;

    fn powi(self, n: i32) -> Self;
    fn powf(self, n: Self) -> Self;

    fn to_degrees(self) -> Self;
    fn to_radians(self) -> Self;
}

#[cfg(feature = "std")]
#[path = "std_impl.rs"]
mod std_impl;

#[cfg(not(feature = "std"))]
#[path = "no_std_impl.rs"]
mod no_std_impl;

/* -------------------- f16 -------------------- */

#[cfg(feature = "f16")]
impl Real for f16 {
    const ZERO:    Self = 0.0;
    const ONE:     Self = 1.0;
    const PI:      Self = core::f16::consts::PI;
    const FRAC_PI_2: Self = core::f16::consts::FRAC_PI_2;
    const EPSILON: Self = f16::EPSILON;
    const INFINITY: Self = f16::INFINITY;
    const NEG_INFINITY: Self = f16::NEG_INFINITY;

    #[inline] fn abs(self) -> Self { self.abs() }
    #[inline] fn signum(self) -> Self { self.signum() }

    #[inline] fn floor(self) -> Self { self.floor() }
    #[inline] fn ceil(self) -> Self { self.ceil() }
    #[inline] fn round(self) -> Self { self.round() }
    #[inline] fn trunc(self) -> Self { self.trunc() }
    #[inline] fn fract(self) -> Self { self.fract() }

    #[inline] fn min(self, other: Self) -> Self { self.min(other) }
    #[inline] fn max(self, other: Self) -> Self { self.max(other) }
    #[inline] fn copysign(self, sign: Self) -> Self { self.copysign(sign) }

    #[inline] fn sqrt(self) -> Self { self.sqrt() }

    #[inline] fn sin(self) -> Self { self.sin() }
    #[inline] fn cos(self) -> Self { self.cos() }
    #[inline] fn tan(self) -> Self { self.tan() }

    #[inline] fn asin(self) -> Self { self.asin() }
    #[inline] fn acos(self) -> Self { self.acos() }
    #[inline] fn atan(self) -> Self { self.atan() }
    #[inline] fn atan2(self, other: Self) -> Self { self.atan2(other) }

    #[inline] fn exp(self) -> Self { self.exp() }
    #[inline] fn exp2(self) -> Self { self.exp2() }
    #[inline] fn ln(self) -> Self { self.ln() }
    #[inline] fn log2(self) -> Self { self.log2() }
    #[inline] fn log10(self) -> Self { self.log10() }

    #[inline] fn powi(self, n: i32) -> Self { self.powi(n) }
    #[inline] fn powf(self, n: Self) -> Self { self.powf(n) }

    #[inline] fn to_degrees(self) -> Self { self.to_degrees() }
    #[inline] fn to_radians(self) -> Self { self.to_radians() }
}

/* -------------------- f128 -------------------- */

#[cfg(feature = "f128")]
impl Real for f128 {
    const ZERO:    Self = 0.0;
    const ONE:     Self = 1.0;
    const PI:      Self = core::f128::consts::PI;
    const FRAC_PI_2: Self = core::f128::consts::FRAC_PI_2;
    const EPSILON: Self = f128::EPSILON;
    const INFINITY: Self = f128::INFINITY;
    const NEG_INFINITY: Self = f128::NEG_INFINITY;

    #[inline] fn abs(self) -> Self { self.abs() }
    #[inline] fn signum(self) -> Self { self.signum() }

    #[inline] fn floor(self) -> Self { self.floor() }
    #[inline] fn ceil(self) -> Self { self.ceil() }
    #[inline] fn round(self) -> Self { self.round() }
    #[inline] fn trunc(self) -> Self { self.trunc() }
    #[inline] fn fract(self) -> Self { self.fract() }

    #[inline] fn min(self, other: Self) -> Self { self.min(other) }
    #[inline] fn max(self, other: Self) -> Self { self.max(other) }
    #[inline] fn copysign(self, sign: Self) -> Self { self.copysign(sign) }

    #[inline] fn sqrt(self) -> Self { self.sqrt() }

    #[inline] fn sin(self) -> Self { self.sin() }
    #[inline] fn cos(self) -> Self { self.cos() }
    #[inline] fn tan(self) -> Self { self.tan() }

    #[inline] fn asin(self) -> Self { self.asin() }
    #[inline] fn acos(self) -> Self { self.acos() }
    #[inline] fn atan(self) -> Self { self.atan() }
    #[inline] fn atan2(self, other: Self) -> Self { self.atan2(other) }

    #[inline] fn exp(self) -> Self { self.exp() }
    #[inline] fn exp2(self) -> Self { self.exp2() }
    #[inline] fn ln(self) -> Self { self.ln() }
    #[inline] fn log2(self) -> Self { self.log2() }
    #[inline] fn log10(self) -> Self { self.log10() }

    #[inline] fn powi(self, n: i32) -> Self { self.powi(n) }
    #[inline] fn powf(self, n: Self) -> Self { self.powf(n) }

    #[inline] fn to_degrees(self) -> Self { self.to_degrees() }
    #[inline] fn to_radians(self) -> Self { self.to_radians() }
}

// Casting
pub trait RealCast: Copy {
    fn from_f32(x: f32) -> Self;
    fn from_f64(x: f64) -> Self;
    fn from_u32(x: u32) -> Self;
    fn from_usize(x: usize) -> Self;

    fn to_f32(self) -> f32;
    fn to_f64(self) -> f64;

    #[cfg(feature = "f16")]
    fn from_f16(x: f16) -> Self;
    #[cfg(feature = "f16")]
    fn to_f16(self) -> f16;

    #[cfg(feature = "f128")]
    fn from_f128(x: f128) -> Self;
    #[cfg(feature = "f128")]
    fn to_f128(self) -> f128;
}

macro_rules! impl_real_cast_float {
    ($t:ty) => {
        impl RealCast for $t {
            #[inline] fn from_f32(x: f32) -> Self { x as $t }
            #[inline] fn from_f64(x: f64) -> Self { x as $t }
            #[inline] fn from_u32(x: u32) -> Self { x as $t }
            #[inline] fn from_usize(x: usize) -> Self { x as $t }

            #[inline] fn to_f32(self) -> f32 { self as f32 }
            #[inline] fn to_f64(self) -> f64 { self as f64 }

            #[cfg(feature = "f16")]
            #[inline] fn from_f16(x: f16) -> Self { x as $t }
            #[cfg(feature = "f16")]
            #[inline] fn to_f16(self) -> f16 { self as f16 }

            #[cfg(feature = "f128")]
            #[inline] fn from_f128(x: f128) -> Self { x as $t }
            #[cfg(feature = "f128")]
            #[inline] fn to_f128(self) -> f128 { self as f128 }
        }
    };
}

impl_real_cast_float!(f32);
impl_real_cast_float!(f64);

#[cfg(feature = "f16")]
impl_real_cast_float!(f16);

#[cfg(feature = "f128")]
impl_real_cast_float!(f128);
