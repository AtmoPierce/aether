use super::Real;

impl Real for f32 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const PI: Self = core::f32::consts::PI;
    const FRAC_PI_2: Self = core::f32::consts::FRAC_PI_2;
    const EPSILON: Self = core::f32::EPSILON;
    const INFINITY: Self = core::f32::INFINITY;
    const NEG_INFINITY: Self = core::f32::NEG_INFINITY;

    #[inline] fn abs(self) -> Self { f32::abs(self) }
    #[inline] fn signum(self) -> Self { f32::signum(self) }
    #[inline] fn floor(self) -> Self { f32::floor(self) }
    #[inline] fn ceil(self) -> Self { f32::ceil(self) }
    #[inline] fn round(self) -> Self { f32::round(self) }
    #[inline] fn trunc(self) -> Self { f32::trunc(self) }
    #[inline] fn fract(self) -> Self { f32::fract(self) }
    #[inline] fn min(self, other: Self) -> Self { f32::min(self, other) }
    #[inline] fn max(self, other: Self) -> Self { f32::max(self, other) }
    #[inline] fn copysign(self, sign: Self) -> Self { f32::copysign(self, sign) }
    #[inline] fn sqrt(self) -> Self { f32::sqrt(self) }
    #[inline] fn sin(self) -> Self { f32::sin(self) }
    #[inline] fn cos(self) -> Self { f32::cos(self) }
    #[inline] fn tan(self) -> Self { f32::tan(self) }
    #[inline] fn asin(self) -> Self { f32::asin(self) }
    #[inline] fn acos(self) -> Self { f32::acos(self) }
    #[inline] fn atan(self) -> Self { f32::atan(self) }
    #[inline] fn atan2(self, other: Self) -> Self { f32::atan2(self, other) }
    #[inline] fn exp(self) -> Self { f32::exp(self) }
    #[inline] fn exp2(self) -> Self { f32::exp2(self) }
    #[inline] fn ln(self) -> Self { f32::ln(self) }
    #[inline] fn log2(self) -> Self { f32::log2(self) }
    #[inline] fn log10(self) -> Self { f32::log10(self) }
    #[inline] fn powi(self, n: i32) -> Self { f32::powi(self, n) }
    #[inline] fn powf(self, n: Self) -> Self { f32::powf(self, n) }
    #[inline] fn to_degrees(self) -> Self { f32::to_degrees(self) }
    #[inline] fn to_radians(self) -> Self { f32::to_radians(self) }
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
    #[inline] fn signum(self) -> Self { f64::signum(self) }
    #[inline] fn floor(self) -> Self { f64::floor(self) }
    #[inline] fn ceil(self) -> Self { f64::ceil(self) }
    #[inline] fn round(self) -> Self { f64::round(self) }
    #[inline] fn trunc(self) -> Self { f64::trunc(self) }
    #[inline] fn fract(self) -> Self { f64::fract(self) }
    #[inline] fn min(self, other: Self) -> Self { f64::min(self, other) }
    #[inline] fn max(self, other: Self) -> Self { f64::max(self, other) }
    #[inline] fn copysign(self, sign: Self) -> Self { f64::copysign(self, sign) }
    #[inline] fn sqrt(self) -> Self { f64::sqrt(self) }
    #[inline] fn sin(self) -> Self { f64::sin(self) }
    #[inline] fn cos(self) -> Self { f64::cos(self) }
    #[inline] fn tan(self) -> Self { f64::tan(self) }
    #[inline] fn asin(self) -> Self { f64::asin(self) }
    #[inline] fn acos(self) -> Self { f64::acos(self) }
    #[inline] fn atan(self) -> Self { f64::atan(self) }
    #[inline] fn atan2(self, other: Self) -> Self { f64::atan2(self, other) }
    #[inline] fn exp(self) -> Self { f64::exp(self) }
    #[inline] fn exp2(self) -> Self { f64::exp2(self) }
    #[inline] fn ln(self) -> Self { f64::ln(self) }
    #[inline] fn log2(self) -> Self { f64::log2(self) }
    #[inline] fn log10(self) -> Self { f64::log10(self) }
    #[inline] fn powi(self, n: i32) -> Self { f64::powi(self, n) }
    #[inline] fn powf(self, n: Self) -> Self { f64::powf(self, n) }
    #[inline] fn to_degrees(self) -> Self { f64::to_degrees(self) }
    #[inline] fn to_radians(self) -> Self { f64::to_radians(self) }
}
