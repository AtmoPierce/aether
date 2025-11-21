use core::fmt::Debug;
use core::ops::{Add, Sub, Mul, Div, Neg};

#[cfg(feature = "f16")]
use core::f16;
#[cfg(feature = "f128")]
use core::f128;

pub trait Real:
    Copy
    + Debug
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    const ZERO: Self;
    const ONE:  Self;
    const PI:   Self;

    // Conversions
    fn from_f32(x: f32) -> Self;
    fn from_f64(x: f64) -> Self;
    fn from_u32(x: u32) -> Self;
    fn from_usize(x: usize) -> Self;

    #[cfg(feature = "f16")]
    fn from_f16(x: f16) -> Self;

    #[cfg(feature = "f128")]
    fn from_f128(x: f128) -> Self;

    // Will add as I use stuff...
    fn abs(self) -> Self;
    fn cos(self) -> Self;
    fn sin(self) -> Self;
}

/* -------------------- f32 -------------------- */

impl Real for f32 {
    const ZERO: Self = 0.0;
    const ONE:  Self = 1.0;
    const PI:   Self = core::f32::consts::PI;

    #[cfg(feature = "f16")]
    #[inline]
    fn from_f16(x: f16) -> Self { x as f32 }

    #[inline]
    fn from_f32(x: f32) -> Self { x }

    #[inline]
    fn from_f64(x: f64) -> Self { x as f32 }

    #[inline]
    fn from_u32(x: u32) -> Self { x as f32 }

    #[inline]
    fn from_usize(x: usize) -> Self { x as f32 }

    #[cfg(feature = "f128")]
    #[inline]
    fn from_f128(x: f128) -> Self { x as f32 }

    #[inline]
    fn abs(self) -> Self { f32::abs(self) }

    #[inline]
    fn cos(self) -> Self { f32::cos(self) }

    #[inline]
    fn sin(self) -> Self { f32::sin(self) }
}

/* -------------------- f64 -------------------- */

impl Real for f64 {
    const ZERO: Self = 0.0;
    const ONE:  Self = 1.0;
    const PI:   Self = core::f64::consts::PI;

    #[cfg(feature = "f16")]
    #[inline]
    fn from_f16(x: f16) -> Self { x as f64 }

    #[inline]
    fn from_f32(x: f32) -> Self { x as f64 }

    #[inline]
    fn from_f64(x: f64) -> Self { x }

    #[inline]
    fn from_u32(x: u32) -> Self { x as f64 }

    #[inline]
    fn from_usize(x: usize) -> Self { x as f64 }

    #[cfg(feature = "f128")]
    #[inline]
    fn from_f128(x: f128) -> Self { x as f64 }

    #[inline]
    fn abs(self) -> Self { f64::abs(self) }

    #[inline]
    fn cos(self) -> Self { f64::cos(self) }

    #[inline]
    fn sin(self) -> Self { f64::sin(self) }
}

/* -------------------- f16 -------------------- */

#[cfg(feature = "f16")]
impl Real for f16 {
    const ZERO: Self = 0.0;
    const ONE:  Self = 1.0;
    const PI:   Self = core::f16::consts::PI;

    #[inline]
    fn from_f16(x: f16) -> Self { x }

    #[inline]
    fn from_f32(x: f32) -> Self { x as f16 }

    #[inline]
    fn from_f64(x: f64) -> Self { x as f16 }

    #[inline]
    fn from_u32(x: u32) -> Self { x as f16 }

    #[inline]
    fn from_usize(x: usize) -> Self { x as f16 }

    #[cfg(feature = "f128")]
    #[inline]
    fn from_f128(x: f128) -> Self { x as f16 }

    #[inline]
    fn abs(self) -> Self { self.abs() }

    #[inline]
    fn cos(self) -> Self { self.cos() }

    #[inline]
    fn sin(self) -> Self { self.sin() }
}

/* -------------------- f128 -------------------- */

#[cfg(feature = "f128")]
impl Real for f128 {
    const ZERO: Self = 0.0;
    const ONE:  Self = 1.0;
    const PI:   Self = core::f128::consts::PI;

    #[cfg(feature = "f16")]
    #[inline]
    fn from_f16(x: f16) -> Self { x as f128 }

    #[inline]
    fn from_f32(x: f32) -> Self { x as f128 }

    #[inline]
    fn from_f64(x: f64) -> Self { x as f128 }

    #[inline]
    fn from_u32(x: u32) -> Self { x as f128 }

    #[inline]
    fn from_usize(x: usize) -> Self { x as f128 }

    #[inline]
    fn from_f128(x: f128) -> Self { x }

    #[inline]
    fn abs(self) -> Self { self.abs() }

    #[inline]
    fn cos(self) -> Self { self.cos() }

    #[inline]
    fn sin(self) -> Self { self.sin() }
}
