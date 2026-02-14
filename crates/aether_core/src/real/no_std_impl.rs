use super::Real;

impl Real for f32 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const PI: Self = core::f32::consts::PI;
    const FRAC_PI_2: Self = core::f32::consts::FRAC_PI_2;
    const EPSILON: Self = core::f32::EPSILON;
    const INFINITY: Self = core::f32::INFINITY;
    const NEG_INFINITY: Self = core::f32::NEG_INFINITY;

    #[inline] fn abs(self) -> Self { libm::fabsf(self) }
    #[inline] fn signum(self) -> Self {
        if self.is_nan() { self } else if self > 0.0 { 1.0 } else if self < 0.0 { -1.0 } else { self }
    }
    #[inline] fn floor(self) -> Self { libm::floorf(self) }
    #[inline] fn ceil(self) -> Self { libm::ceilf(self) }
    #[inline] fn round(self) -> Self { libm::roundf(self) }
    #[inline] fn trunc(self) -> Self { libm::truncf(self) }
    #[inline] fn fract(self) -> Self { self - libm::truncf(self) }
    #[inline] fn min(self, other: Self) -> Self {
        if self.is_nan() { other } else if other.is_nan() || self < other { self } else { other }
    }
    #[inline] fn max(self, other: Self) -> Self {
        if self.is_nan() { other } else if other.is_nan() || self > other { self } else { other }
    }
    #[inline] fn copysign(self, sign: Self) -> Self { libm::copysignf(self, sign) }
    #[inline] fn sqrt(self) -> Self { libm::sqrtf(self) }
    #[inline] fn sin(self) -> Self { libm::sinf(self) }
    #[inline] fn cos(self) -> Self { libm::cosf(self) }
    #[inline] fn tan(self) -> Self { libm::tanf(self) }
    #[inline] fn asin(self) -> Self { libm::asinf(self) }
    #[inline] fn acos(self) -> Self { libm::acosf(self) }
    #[inline] fn atan(self) -> Self { libm::atanf(self) }
    #[inline] fn atan2(self, other: Self) -> Self { libm::atan2f(self, other) }
    #[inline] fn exp(self) -> Self { libm::expf(self) }
    #[inline] fn exp2(self) -> Self { libm::exp2f(self) }
    #[inline] fn ln(self) -> Self { libm::logf(self) }
    #[inline] fn log2(self) -> Self { libm::log2f(self) }
    #[inline] fn log10(self) -> Self { libm::log10f(self) }
    #[inline] fn powi(self, n: i32) -> Self { libm::powf(self, n as f32) }
    #[inline] fn powf(self, n: Self) -> Self { libm::powf(self, n) }
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

    #[inline] fn abs(self) -> Self { libm::fabs(self) }
    #[inline] fn signum(self) -> Self {
        if self.is_nan() { self } else if self > 0.0 { 1.0 } else if self < 0.0 { -1.0 } else { self }
    }
    #[inline] fn floor(self) -> Self { libm::floor(self) }
    #[inline] fn ceil(self) -> Self { libm::ceil(self) }
    #[inline] fn round(self) -> Self { libm::round(self) }
    #[inline] fn trunc(self) -> Self { libm::trunc(self) }
    #[inline] fn fract(self) -> Self { self - libm::trunc(self) }
    #[inline] fn min(self, other: Self) -> Self {
        if self.is_nan() { other } else if other.is_nan() || self < other { self } else { other }
    }
    #[inline] fn max(self, other: Self) -> Self {
        if self.is_nan() { other } else if other.is_nan() || self > other { self } else { other }
    }
    #[inline] fn copysign(self, sign: Self) -> Self { libm::copysign(self, sign) }
    #[inline] fn sqrt(self) -> Self { libm::sqrt(self) }
    #[inline] fn sin(self) -> Self { libm::sin(self) }
    #[inline] fn cos(self) -> Self { libm::cos(self) }
    #[inline] fn tan(self) -> Self { libm::tan(self) }
    #[inline] fn asin(self) -> Self { libm::asin(self) }
    #[inline] fn acos(self) -> Self { libm::acos(self) }
    #[inline] fn atan(self) -> Self { libm::atan(self) }
    #[inline] fn atan2(self, other: Self) -> Self { libm::atan2(self, other) }
    #[inline] fn exp(self) -> Self { libm::exp(self) }
    #[inline] fn exp2(self) -> Self { libm::exp2(self) }
    #[inline] fn ln(self) -> Self { libm::log(self) }
    #[inline] fn log2(self) -> Self { libm::log2(self) }
    #[inline] fn log10(self) -> Self { libm::log10(self) }
    #[inline] fn powi(self, n: i32) -> Self { libm::pow(self, n as f64) }
    #[inline] fn powf(self, n: Self) -> Self { libm::pow(self, n) }
    #[inline] fn to_degrees(self) -> Self { self * (180.0 / core::f64::consts::PI) }
    #[inline] fn to_radians(self) -> Self { self * (core::f64::consts::PI / 180.0) }
}
