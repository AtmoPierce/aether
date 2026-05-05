use super::{Real, RealCast};

impl Real for f32 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const PI: Self = core::f32::consts::PI;
    const FRAC_PI_2: Self = core::f32::consts::FRAC_PI_2;
    const EPSILON: Self = core::f32::EPSILON;
    const INFINITY: Self = core::f32::INFINITY;
    const NEG_INFINITY: Self = core::f32::NEG_INFINITY;

    #[inline] fn abs(self) -> Self { f32::abs(self) }
    #[inline] fn signum(self) -> Self {
        if self.is_nan() { self } else if self > 0.0 { 1.0 } else if self < 0.0 { -1.0 } else { self }
    }
    #[inline] fn floor(self) -> Self {
        #[cfg(feature = "std")]
        { f32::floor(self) }
        #[cfg(not(feature = "std"))]
        { libm::floorf(self) }
    }
    #[inline] fn ceil(self) -> Self {
        #[cfg(feature = "std")]
        { f32::ceil(self) }
        #[cfg(not(feature = "std"))]
        { libm::ceilf(self) }
    }
    #[inline] fn round(self) -> Self {
        #[cfg(feature = "std")]
        { f32::round(self) }
        #[cfg(not(feature = "std"))]
        { libm::roundf(self) }
    }
    #[inline] fn trunc(self) -> Self {
        #[cfg(feature = "std")]
        { f32::trunc(self) }
        #[cfg(not(feature = "std"))]
        { libm::truncf(self) }
    }
    #[inline] fn fract(self) -> Self {
        #[cfg(feature = "std")]
        { f32::fract(self) }
        #[cfg(not(feature = "std"))]
        { self - libm::truncf(self) }
    }
    #[inline] fn min(self, other: Self) -> Self {
        if self.is_nan() { other } else if other.is_nan() || self < other { self } else { other }
    }
    #[inline] fn max(self, other: Self) -> Self {
        if self.is_nan() { other } else if other.is_nan() || self > other { self } else { other }
    }
    #[inline] fn copysign(self, sign: Self) -> Self {
        #[cfg(feature = "std")]
        { f32::copysign(self, sign) }
        #[cfg(not(feature = "std"))]
        { libm::copysignf(self, sign) }
    }
    #[inline] fn sqrt(self) -> Self {
        #[cfg(feature = "std")]
        { f32::sqrt(self) }
        #[cfg(not(feature = "std"))]
        { libm::sqrtf(self) }
    }

    #[inline] fn sin(self) -> Self { fpx_core::trig::sinf(self) }
    #[inline] fn cos(self) -> Self { fpx_core::trig::cosf(self) }
    #[inline] fn tan(self) -> Self { fpx_core::trig::tanf(self) }

    #[inline] fn asin(self) -> Self { fpx_core::trig::asinf(self) }
    #[inline] fn acos(self) -> Self { fpx_core::trig::acosf(self) }
    #[inline] fn atan(self) -> Self { fpx_core::trig::atanf(self) }
    #[inline] fn atan2(self, other: Self) -> Self {
        #[cfg(feature = "std")]
        { f32::atan2(self, other) }
        #[cfg(not(feature = "std"))]
        { libm::atan2f(self, other) }
    }

    #[inline] fn exp(self) -> Self { fpx_core::trig::expf(self) }
    #[inline] fn exp2(self) -> Self { fpx_core::trig::expf(self * core::f32::consts::LN_2) }
    #[inline] fn ln(self) -> Self { fpx_core::trig::logf(self) }
    #[inline] fn log2(self) -> Self { fpx_core::trig::logf(self) / core::f32::consts::LN_2 }
    #[inline] fn log10(self) -> Self { fpx_core::trig::logf(self) / core::f32::consts::LN_10 }

    #[inline] fn sinh(self) -> Self { fpx_core::trig::sinh(self) }
    #[inline] fn cosh(self) -> Self { fpx_core::trig::cosh(self) }
    #[inline] fn tanh(self) -> Self { fpx_core::trig::tanh(self) }

    #[inline] fn exp_m1(self) -> Self { fpx_core::trig::expm1(self) }
    #[inline] fn ln_1p(self) -> Self { fpx_core::trig::log1p(self) }

    #[inline] fn powi(self, n: i32) -> Self {
        #[cfg(feature = "std")]
        { f32::powi(self, n) }
        #[cfg(not(feature = "std"))]
        { libm::powf(self, n as f32) }
    }
    #[inline] fn powf(self, n: Self) -> Self {
        #[cfg(feature = "std")]
        { f32::powf(self, n) }
        #[cfg(not(feature = "std"))]
        { libm::powf(self, n) }
    }

    #[inline] fn to_degrees(self) -> Self { self * (180.0 / core::f32::consts::PI) }
    #[inline] fn to_radians(self) -> Self { self * (core::f32::consts::PI / 180.0) }
}

impl Real for f64 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const PI: Self = core::f64::consts::PI;
    const FRAC_PI_2: Self = core::f64::consts::FRAC_PI_2;
    const EPSILON: Self = core::f64::EPSILON;
    const INFINITY: Self = core::f64::INFINITY;
    const NEG_INFINITY: Self = core::f64::NEG_INFINITY;

    #[inline] fn abs(self) -> Self { f64::abs(self) }
    #[inline] fn signum(self) -> Self {
        if self.is_nan() { self } else if self > 0.0 { 1.0 } else if self < 0.0 { -1.0 } else { self }
    }
    #[inline] fn floor(self) -> Self {
        #[cfg(feature = "std")]
        { f64::floor(self) }
        #[cfg(not(feature = "std"))]
        { libm::floor(self) }
    }
    #[inline] fn ceil(self) -> Self {
        #[cfg(feature = "std")]
        { f64::ceil(self) }
        #[cfg(not(feature = "std"))]
        { libm::ceil(self) }
    }
    #[inline] fn round(self) -> Self {
        #[cfg(feature = "std")]
        { f64::round(self) }
        #[cfg(not(feature = "std"))]
        { libm::round(self) }
    }
    #[inline] fn trunc(self) -> Self {
        #[cfg(feature = "std")]
        { f64::trunc(self) }
        #[cfg(not(feature = "std"))]
        { libm::trunc(self) }
    }
    #[inline] fn fract(self) -> Self {
        #[cfg(feature = "std")]
        { f64::fract(self) }
        #[cfg(not(feature = "std"))]
        { self - libm::trunc(self) }
    }
    #[inline] fn min(self, other: Self) -> Self {
        if self.is_nan() { other } else if other.is_nan() || self < other { self } else { other }
    }
    #[inline] fn max(self, other: Self) -> Self {
        if self.is_nan() { other } else if other.is_nan() || self > other { self } else { other }
    }
    #[inline] fn copysign(self, sign: Self) -> Self {
        #[cfg(feature = "std")]
        { f64::copysign(self, sign) }
        #[cfg(not(feature = "std"))]
        { libm::copysign(self, sign) }
    }
    #[inline] fn sqrt(self) -> Self {
        #[cfg(feature = "std")]
        { f64::sqrt(self) }
        #[cfg(not(feature = "std"))]
        { libm::sqrt(self) }
    }

    #[inline] fn sin(self) -> Self { fpx_core::trig::sin(self) }
    #[inline] fn cos(self) -> Self { fpx_core::trig::cos(self) }
    #[inline] fn tan(self) -> Self { fpx_core::trig::tan(self) }

    #[inline] fn asin(self) -> Self { fpx_core::trig::asin(self) }
    #[inline] fn acos(self) -> Self { fpx_core::trig::acos(self) }
    #[inline] fn atan(self) -> Self { fpx_core::trig::atan(self) }
    #[inline] fn atan2(self, other: Self) -> Self {
        #[cfg(feature = "std")]
        { f64::atan2(self, other) }
        #[cfg(not(feature = "std"))]
        { libm::atan2(self, other) }
    }

    #[inline] fn exp(self) -> Self { fpx_core::trig::exp(self) }
    #[inline] fn exp2(self) -> Self { fpx_core::trig::exp(self * core::f64::consts::LN_2) }
    #[inline] fn ln(self) -> Self { fpx_core::trig::log(self) }
    #[inline] fn log2(self) -> Self { fpx_core::trig::log(self) / core::f64::consts::LN_2 }
    #[inline] fn log10(self) -> Self { fpx_core::trig::log(self) / core::f64::consts::LN_10 }

    #[inline] fn sinh(self) -> Self { fpx_core::trig::sinh(self) }
    #[inline] fn cosh(self) -> Self { fpx_core::trig::cosh(self) }
    #[inline] fn tanh(self) -> Self { fpx_core::trig::tanh(self) }

    #[inline] fn exp_m1(self) -> Self { fpx_core::trig::expm1(self) }
    #[inline] fn ln_1p(self) -> Self { fpx_core::trig::log1p(self) }

    #[inline] fn powi(self, n: i32) -> Self {
        #[cfg(feature = "std")]
        { f64::powi(self, n) }
        #[cfg(not(feature = "std"))]
        { libm::pow(self, n as f64) }
    }
    #[inline] fn powf(self, n: Self) -> Self {
        #[cfg(feature = "std")]
        { f64::powf(self, n) }
        #[cfg(not(feature = "std"))]
        { libm::pow(self, n) }
    }

    #[inline] fn to_degrees(self) -> Self { self * (180.0 / core::f64::consts::PI) }
    #[inline] fn to_radians(self) -> Self { self * (core::f64::consts::PI / 180.0) }
}

#[cfg(feature = "f16")]
impl Real for f16 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const PI: Self = core::f16::consts::PI;
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

    #[inline] fn sin(self) -> Self { fpx_core::trig::sins(self) }
    #[inline] fn cos(self) -> Self { fpx_core::trig::coss(self) }
    #[inline] fn tan(self) -> Self { fpx_core::trig::tans(self) }

    #[inline] fn asin(self) -> Self { fpx_core::trig::asins(self) }
    #[inline] fn acos(self) -> Self { fpx_core::trig::acoss(self) }
    #[inline] fn atan(self) -> Self { fpx_core::trig::atans(self) }
    #[inline] fn atan2(self, other: Self) -> Self { self.atan2(other) }

    #[inline] fn exp(self) -> Self { fpx_core::trig::exps(self) }
    #[inline] fn exp2(self) -> Self { fpx_core::trig::exps(self * f16::from_f32(core::f32::consts::LN_2)) }
    #[inline] fn ln(self) -> Self { fpx_core::trig::logs(self) }
    #[inline] fn log2(self) -> Self { fpx_core::trig::logs(self) / f16::from_f32(core::f32::consts::LN_2) }
    #[inline] fn log10(self) -> Self { fpx_core::trig::logs(self) / f16::from_f32(core::f32::consts::LN_10) }

    #[inline] fn sinh(self) -> Self { fpx_core::trig::sinhs(self) }
    #[inline] fn cosh(self) -> Self { fpx_core::trig::coshs(self) }
    #[inline] fn tanh(self) -> Self { fpx_core::trig::tanhs(self) }

    #[inline] fn exp_m1(self) -> Self { fpx_core::trig::expm1s(self) }
    #[inline] fn ln_1p(self) -> Self { fpx_core::trig::log1ps(self) }
    
    
    #[inline] fn powi(self, n: i32) -> Self { self.powi(n) }
    #[inline] fn powf(self, n: Self) -> Self { self.powf(n) }

    #[inline] fn to_degrees(self) -> Self { self.to_degrees() }
    #[inline] fn to_radians(self) -> Self { self.to_radians() }
}

#[cfg(feature = "f128")]
impl Real for f128 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const PI: Self = core::f128::consts::PI;
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

    #[inline] fn sin(self) -> Self { fpx_core::trig::sinl(self) }
    #[inline] fn cos(self) -> Self { fpx_core::trig::cosl(self) }
    #[inline] fn tan(self) -> Self { fpx_core::trig::tanl(self) }

    #[inline] fn asin(self) -> Self { fpx_core::trig::asinl(self) }
    #[inline] fn acos(self) -> Self { fpx_core::trig::acosl(self) }
    #[inline] fn atan(self) -> Self { fpx_core::trig::atanl(self) }
    #[inline] fn atan2(self, other: Self) -> Self { self.atan2(other) }

    #[inline] fn exp(self) -> Self { fpx_core::trig::expl(self) }
    #[inline] fn exp2(self) -> Self { fpx_core::trig::expl(self * core::f128::consts::LN_2) }
    #[inline] fn ln(self) -> Self { fpx_core::trig::logl(self) }
    #[inline] fn log2(self) -> Self { fpx_core::trig::logl(self) / core::f128::consts::LN_2 }
    #[inline] fn log10(self) -> Self { fpx_core::trig::logl(self) / core::f128::consts::LN_10 }

    #[inline] fn sinh(self) -> Self { fpx_core::trig::sinhl(self) }
    #[inline] fn cosh(self) -> Self { fpx_core::trig::coshl(self) }
    #[inline] fn tanh(self) -> Self { fpx_core::trig::tanhl(self) }

    #[inline] fn exp_m1(self) -> Self { fpx_core::trig::expm1l(self) }
    #[inline] fn ln_1p(self) -> Self { fpx_core::trig::log1pl(self) }

    #[inline] fn powi(self, n: i32) -> Self { self.powi(n) }
    #[inline] fn powf(self, n: Self) -> Self { self.powf(n) }

    #[inline] fn to_degrees(self) -> Self { self.to_degrees() }
    #[inline] fn to_radians(self) -> Self { self.to_radians() }
}