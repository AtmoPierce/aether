use super::{Real, RealCast};

impl Real for f32 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const PI: Self = core::f32::consts::PI;
    const FRAC_PI_2: Self = core::f32::consts::FRAC_PI_2;
    const EPSILON: Self = core::f32::EPSILON;
    const INFINITY: Self = core::f32::INFINITY;
    const NEG_INFINITY: Self = core::f32::NEG_INFINITY;

    #[inline]
    fn abs(self) -> Self {
        f32::abs(self)
    }
    #[inline]
    fn signum(self) -> Self {
        if self.is_nan() {
            self
        } else if self > 0.0 {
            1.0
        } else if self < 0.0 {
            -1.0
        } else {
            self
        }
    }
    #[inline]
    fn floor(self) -> Self {
        #[cfg(feature = "std")]
        {
            f32::floor(self)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::floorf(self)
        }
    }
    #[inline]
    fn ceil(self) -> Self {
        #[cfg(feature = "std")]
        {
            f32::ceil(self)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::ceilf(self)
        }
    }
    #[inline]
    fn round(self) -> Self {
        #[cfg(feature = "std")]
        {
            f32::round(self)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::roundf(self)
        }
    }
    #[inline]
    fn trunc(self) -> Self {
        #[cfg(feature = "std")]
        {
            f32::trunc(self)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::truncf(self)
        }
    }
    #[inline]
    fn fract(self) -> Self {
        #[cfg(feature = "std")]
        {
            f32::fract(self)
        }
        #[cfg(not(feature = "std"))]
        {
            self - libm::truncf(self)
        }
    }
    #[inline]
    fn min(self, other: Self) -> Self {
        if self.is_nan() {
            other
        } else if other.is_nan() || self < other {
            self
        } else {
            other
        }
    }
    #[inline]
    fn max(self, other: Self) -> Self {
        if self.is_nan() {
            other
        } else if other.is_nan() || self > other {
            self
        } else {
            other
        }
    }
    #[inline]
    fn copysign(self, sign: Self) -> Self {
        #[cfg(feature = "std")]
        {
            f32::copysign(self, sign)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::copysignf(self, sign)
        }
    }
    #[inline]
    fn sqrt(self) -> Self {
        #[cfg(feature = "std")]
        {
            f32::sqrt(self)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::sqrtf(self)
        }
    }

    #[inline]
    fn sin(self) -> Self {
        fpx_core::trig::sinf(self)
    }
    #[inline]
    fn cos(self) -> Self {
        fpx_core::trig::cosf(self)
    }
    #[inline]
    fn tan(self) -> Self {
        fpx_core::trig::tanf(self)
    }

    #[inline]
    fn asin(self) -> Self {
        fpx_core::trig::asinf(self)
    }
    #[inline]
    fn acos(self) -> Self {
        fpx_core::trig::acosf(self)
    }
    #[inline]
    fn atan(self) -> Self {
        fpx_core::trig::atanf(self)
    }
    #[inline]
    fn atan2(self, other: Self) -> Self {
        #[cfg(feature = "std")]
        {
            f32::atan2(self, other)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::atan2f(self, other)
        }
    }

    #[inline]
    fn exp(self) -> Self {
        fpx_core::trig::expf(self)
    }
    #[inline]
    fn exp2(self) -> Self {
        fpx_core::trig::expf(self * core::f32::consts::LN_2)
    }
    #[inline]
    fn ln(self) -> Self {
        fpx_core::trig::logf(self)
    }
    #[inline]
    fn log2(self) -> Self {
        fpx_core::trig::logf(self) / core::f32::consts::LN_2
    }
    #[inline]
    fn log10(self) -> Self {
        fpx_core::trig::logf(self) / core::f32::consts::LN_10
    }

    #[inline]
    fn sinh(self) -> Self {
        fpx_core::trig::sinhf(self)
    }
    #[inline]
    fn cosh(self) -> Self {
        fpx_core::trig::coshf(self)
    }
    #[inline]
    fn tanh(self) -> Self {
        fpx_core::trig::tanhf(self)
    }

    #[inline]
    fn exp_m1(self) -> Self {
        fpx_core::trig::expm1f(self)
    }
    #[inline]
    fn ln_1p(self) -> Self {
        fpx_core::trig::log1pf(self)
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        #[cfg(feature = "std")]
        {
            f32::powi(self, n)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::powf(self, n as f32)
        }
    }
    #[inline]
    fn powf(self, n: Self) -> Self {
        #[cfg(feature = "std")]
        {
            f32::powf(self, n)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::powf(self, n)
        }
    }

    #[inline]
    fn to_degrees(self) -> Self {
        self * (180.0 / core::f32::consts::PI)
    }
    #[inline]
    fn to_radians(self) -> Self {
        self * (core::f32::consts::PI / 180.0)
    }
}

impl Real for f64 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const PI: Self = core::f64::consts::PI;
    const FRAC_PI_2: Self = core::f64::consts::FRAC_PI_2;
    const EPSILON: Self = core::f64::EPSILON;
    const INFINITY: Self = core::f64::INFINITY;
    const NEG_INFINITY: Self = core::f64::NEG_INFINITY;

    #[inline]
    fn abs(self) -> Self {
        f64::abs(self)
    }
    #[inline]
    fn signum(self) -> Self {
        if self.is_nan() {
            self
        } else if self > 0.0 {
            1.0
        } else if self < 0.0 {
            -1.0
        } else {
            self
        }
    }
    #[inline]
    fn floor(self) -> Self {
        #[cfg(feature = "std")]
        {
            f64::floor(self)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::floor(self)
        }
    }
    #[inline]
    fn ceil(self) -> Self {
        #[cfg(feature = "std")]
        {
            f64::ceil(self)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::ceil(self)
        }
    }
    #[inline]
    fn round(self) -> Self {
        #[cfg(feature = "std")]
        {
            f64::round(self)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::round(self)
        }
    }
    #[inline]
    fn trunc(self) -> Self {
        #[cfg(feature = "std")]
        {
            f64::trunc(self)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::trunc(self)
        }
    }
    #[inline]
    fn fract(self) -> Self {
        #[cfg(feature = "std")]
        {
            f64::fract(self)
        }
        #[cfg(not(feature = "std"))]
        {
            self - libm::trunc(self)
        }
    }
    #[inline]
    fn min(self, other: Self) -> Self {
        if self.is_nan() {
            other
        } else if other.is_nan() || self < other {
            self
        } else {
            other
        }
    }
    #[inline]
    fn max(self, other: Self) -> Self {
        if self.is_nan() {
            other
        } else if other.is_nan() || self > other {
            self
        } else {
            other
        }
    }
    #[inline]
    fn copysign(self, sign: Self) -> Self {
        #[cfg(feature = "std")]
        {
            f64::copysign(self, sign)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::copysign(self, sign)
        }
    }
    #[inline]
    fn sqrt(self) -> Self {
        #[cfg(feature = "std")]
        {
            f64::sqrt(self)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::sqrt(self)
        }
    }

    #[inline]
    fn sin(self) -> Self {
        fpx_core::trig::sin(self)
    }
    #[inline]
    fn cos(self) -> Self {
        fpx_core::trig::cos(self)
    }
    #[inline]
    fn tan(self) -> Self {
        fpx_core::trig::tan(self)
    }

    #[inline]
    fn asin(self) -> Self {
        fpx_core::trig::asin(self)
    }
    #[inline]
    fn acos(self) -> Self {
        fpx_core::trig::acos(self)
    }
    #[inline]
    fn atan(self) -> Self {
        fpx_core::trig::atan(self)
    }
    #[inline]
    fn atan2(self, other: Self) -> Self {
        #[cfg(feature = "std")]
        {
            f64::atan2(self, other)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::atan2(self, other)
        }
    }

    #[inline]
    fn exp(self) -> Self {
        fpx_core::trig::exp(self)
    }
    #[inline]
    fn exp2(self) -> Self {
        fpx_core::trig::exp(self * core::f64::consts::LN_2)
    }
    #[inline]
    fn ln(self) -> Self {
        fpx_core::trig::log(self)
    }
    #[inline]
    fn log2(self) -> Self {
        fpx_core::trig::log(self) / core::f64::consts::LN_2
    }
    #[inline]
    fn log10(self) -> Self {
        fpx_core::trig::log(self) / core::f64::consts::LN_10
    }

    #[inline]
    fn sinh(self) -> Self {
        fpx_core::trig::sinh(self)
    }
    #[inline]
    fn cosh(self) -> Self {
        fpx_core::trig::cosh(self)
    }
    #[inline]
    fn tanh(self) -> Self {
        fpx_core::trig::tanh(self)
    }

    #[inline]
    fn exp_m1(self) -> Self {
        fpx_core::trig::expm1(self)
    }
    #[inline]
    fn ln_1p(self) -> Self {
        fpx_core::trig::log1p(self)
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        #[cfg(feature = "std")]
        {
            f64::powi(self, n)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::pow(self, n as f64)
        }
    }
    #[inline]
    fn powf(self, n: Self) -> Self {
        #[cfg(feature = "std")]
        {
            f64::powf(self, n)
        }
        #[cfg(not(feature = "std"))]
        {
            libm::pow(self, n)
        }
    }

    #[inline]
    fn to_degrees(self) -> Self {
        self * (180.0 / core::f64::consts::PI)
    }
    #[inline]
    fn to_radians(self) -> Self {
        self * (core::f64::consts::PI / 180.0)
    }
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

    #[inline]
    fn abs(self) -> Self {
        self.abs()
    }
    #[inline]
    fn signum(self) -> Self {
        self.signum()
    }

    #[inline]
    fn floor(self) -> Self {
        self.floor()
    }
    #[inline]
    fn ceil(self) -> Self {
        self.ceil()
    }
    #[inline]
    fn round(self) -> Self {
        self.round()
    }
    #[inline]
    fn trunc(self) -> Self {
        self.trunc()
    }
    #[inline]
    fn fract(self) -> Self {
        self.fract()
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        self.min(other)
    }
    #[inline]
    fn max(self, other: Self) -> Self {
        self.max(other)
    }
    #[inline]
    fn copysign(self, sign: Self) -> Self {
        self.copysign(sign)
    }

    #[inline]
    fn sqrt(self) -> Self {
        self.sqrt()
    }

    #[inline]
    fn sin(self) -> Self {
        fpx_core::trig::sins(self)
    }
    #[inline]
    fn cos(self) -> Self {
        fpx_core::trig::coss(self)
    }
    #[inline]
    fn tan(self) -> Self {
        fpx_core::trig::tans(self)
    }

    #[inline]
    fn asin(self) -> Self {
        fpx_core::trig::asins(self)
    }
    #[inline]
    fn acos(self) -> Self {
        fpx_core::trig::acoss(self)
    }
    #[inline]
    fn atan(self) -> Self {
        fpx_core::trig::atans(self)
    }
    #[inline]
    fn atan2(self, other: Self) -> Self {
        self.atan2(other)
    }

    #[inline]
    fn exp(self) -> Self {
        fpx_core::trig::exps(self)
    }
    #[inline]
    fn exp2(self) -> Self {
        fpx_core::trig::exps(self * f16::from_f32(core::f32::consts::LN_2))
    }
    #[inline]
    fn ln(self) -> Self {
        fpx_core::trig::logs(self)
    }
    #[inline]
    fn log2(self) -> Self {
        fpx_core::trig::logs(self) / f16::from_f32(core::f32::consts::LN_2)
    }
    #[inline]
    fn log10(self) -> Self {
        fpx_core::trig::logs(self) / f16::from_f32(core::f32::consts::LN_10)
    }

    #[inline]
    fn sinh(self) -> Self {
        fpx_core::trig::sinhs(self)
    }
    #[inline]
    fn cosh(self) -> Self {
        fpx_core::trig::coshs(self)
    }
    #[inline]
    fn tanh(self) -> Self {
        fpx_core::trig::tanhs(self)
    }

    #[inline]
    fn exp_m1(self) -> Self {
        fpx_core::trig::expm1s(self)
    }
    #[inline]
    fn ln_1p(self) -> Self {
        fpx_core::trig::log1ps(self)
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        self.powi(n)
    }
    #[inline]
    fn powf(self, n: Self) -> Self {
        self.powf(n)
    }

    #[inline]
    fn to_degrees(self) -> Self {
        self.to_degrees()
    }
    #[inline]
    fn to_radians(self) -> Self {
        self.to_radians()
    }
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

    #[inline]
    fn abs(self) -> Self {
        self.abs()
    }
    #[inline]
    fn signum(self) -> Self {
        self.signum()
    }

    #[inline]
    fn floor(self) -> Self {
        self.floor()
    }
    #[inline]
    fn ceil(self) -> Self {
        self.ceil()
    }
    #[inline]
    fn round(self) -> Self {
        self.round()
    }
    #[inline]
    fn trunc(self) -> Self {
        self.trunc()
    }
    #[inline]
    fn fract(self) -> Self {
        self.fract()
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        self.min(other)
    }
    #[inline]
    fn max(self, other: Self) -> Self {
        self.max(other)
    }
    #[inline]
    fn copysign(self, sign: Self) -> Self {
        self.copysign(sign)
    }

    #[inline]
    fn sqrt(self) -> Self {
        self.sqrt()
    }

    #[inline]
    fn sin(self) -> Self {
        fpx_core::trig::sinl(self)
    }
    #[inline]
    fn cos(self) -> Self {
        fpx_core::trig::cosl(self)
    }
    #[inline]
    fn tan(self) -> Self {
        fpx_core::trig::tanl(self)
    }

    #[inline]
    fn asin(self) -> Self {
        fpx_core::trig::asinl(self)
    }
    #[inline]
    fn acos(self) -> Self {
        fpx_core::trig::acosl(self)
    }
    #[inline]
    fn atan(self) -> Self {
        fpx_core::trig::atanl(self)
    }
    #[inline]
    fn atan2(self, other: Self) -> Self {
        self.atan2(other)
    }

    #[inline]
    fn exp(self) -> Self {
        fpx_core::trig::expl(self)
    }
    #[inline]
    fn exp2(self) -> Self {
        fpx_core::trig::expl(self * core::f128::consts::LN_2)
    }
    #[inline]
    fn ln(self) -> Self {
        fpx_core::trig::logl(self)
    }
    #[inline]
    fn log2(self) -> Self {
        fpx_core::trig::logl(self) / core::f128::consts::LN_2
    }
    #[inline]
    fn log10(self) -> Self {
        fpx_core::trig::logl(self) / core::f128::consts::LN_10
    }

    #[inline]
    fn sinh(self) -> Self {
        fpx_core::trig::sinhl(self)
    }
    #[inline]
    fn cosh(self) -> Self {
        fpx_core::trig::coshl(self)
    }
    #[inline]
    fn tanh(self) -> Self {
        fpx_core::trig::tanhl(self)
    }

    #[inline]
    fn exp_m1(self) -> Self {
        fpx_core::trig::expm1l(self)
    }
    #[inline]
    fn ln_1p(self) -> Self {
        fpx_core::trig::log1pl(self)
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        self.powi(n)
    }
    #[inline]
    fn powf(self, n: Self) -> Self {
        self.powf(n)
    }

    #[inline]
    fn to_degrees(self) -> Self {
        self.to_degrees()
    }
    #[inline]
    fn to_radians(self) -> Self {
        self.to_radians()
    }
}

#[cfg(feature = "mx")]
trait MxRealFunctions: RealCast {
    const ZERO: Self;
    const ONE: Self;
    const PI: Self;
    const FRAC_PI_2: Self;
    const EPSILON: Self;
    const INFINITY: Self;
    const NEG_INFINITY: Self;

    fn fpx_sin(self) -> Self;
    fn fpx_cos(self) -> Self;
    fn fpx_tan(self) -> Self;
    fn fpx_asin(self) -> Self;
    fn fpx_acos(self) -> Self;
    fn fpx_atan(self) -> Self;
    fn fpx_exp(self) -> Self;
    fn fpx_log(self) -> Self;
    fn fpx_sinh(self) -> Self;
    fn fpx_cosh(self) -> Self;
    fn fpx_tanh(self) -> Self;
    fn fpx_expm1(self) -> Self;
    fn fpx_log1p(self) -> Self;
}

#[cfg(feature = "mx")]
impl<T> Real for T
where
    T: MxRealFunctions
        + core::fmt::Debug
        + PartialEq
        + PartialOrd
        + core::ops::Add<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Mul<Output = T>
        + core::ops::Div<Output = T>
        + core::ops::Neg<Output = T>,
{
    const ZERO: Self = <Self as MxRealFunctions>::ZERO;
    const ONE: Self = <Self as MxRealFunctions>::ONE;
    const PI: Self = <Self as MxRealFunctions>::PI;
    const FRAC_PI_2: Self = <Self as MxRealFunctions>::FRAC_PI_2;
    const EPSILON: Self = <Self as MxRealFunctions>::EPSILON;
    const INFINITY: Self = <Self as MxRealFunctions>::INFINITY;
    const NEG_INFINITY: Self = <Self as MxRealFunctions>::NEG_INFINITY;

    #[inline]
    fn abs(self) -> Self {
        Self::from_f32(self.to_f32().abs())
    }

    #[inline]
    fn signum(self) -> Self {
        let value = self.to_f32();
        if value.is_nan() {
            self
        } else if value > 0.0 {
            Self::ONE
        } else if value < 0.0 {
            -Self::ONE
        } else {
            self
        }
    }

    #[inline]
    fn floor(self) -> Self {
        Self::from_f32(self.to_f32().floor())
    }

    #[inline]
    fn ceil(self) -> Self {
        Self::from_f32(self.to_f32().ceil())
    }

    #[inline]
    fn round(self) -> Self {
        Self::from_f32(self.to_f32().round())
    }

    #[inline]
    fn trunc(self) -> Self {
        Self::from_f32(self.to_f32().trunc())
    }

    #[inline]
    fn fract(self) -> Self {
        Self::from_f32(self.to_f32().fract())
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        Self::from_f32(self.to_f32().min(other.to_f32()))
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        Self::from_f32(self.to_f32().max(other.to_f32()))
    }

    #[inline]
    fn copysign(self, sign: Self) -> Self {
        Self::from_f32(self.to_f32().copysign(sign.to_f32()))
    }

    #[inline]
    fn sqrt(self) -> Self {
        Self::from_f32(self.to_f32().sqrt())
    }

    #[inline]
    fn sin(self) -> Self {
        self.fpx_sin()
    }

    #[inline]
    fn cos(self) -> Self {
        self.fpx_cos()
    }

    #[inline]
    fn tan(self) -> Self {
        self.fpx_tan()
    }

    #[inline]
    fn asin(self) -> Self {
        self.fpx_asin()
    }

    #[inline]
    fn acos(self) -> Self {
        self.fpx_acos()
    }

    #[inline]
    fn atan(self) -> Self {
        self.fpx_atan()
    }

    #[inline]
    fn atan2(self, other: Self) -> Self {
        Self::from_f32(self.to_f32().atan2(other.to_f32()))
    }

    #[inline]
    fn exp(self) -> Self {
        self.fpx_exp()
    }

    #[inline]
    fn exp2(self) -> Self {
        Self::from_f32(self.to_f32() * core::f32::consts::LN_2).fpx_exp()
    }

    #[inline]
    fn ln(self) -> Self {
        self.fpx_log()
    }

    #[inline]
    fn log2(self) -> Self {
        self.fpx_log() / Self::from_f32(core::f32::consts::LN_2)
    }

    #[inline]
    fn log10(self) -> Self {
        self.fpx_log() / Self::from_f32(core::f32::consts::LN_10)
    }

    #[inline]
    fn sinh(self) -> Self {
        self.fpx_sinh()
    }

    #[inline]
    fn cosh(self) -> Self {
        self.fpx_cosh()
    }

    #[inline]
    fn tanh(self) -> Self {
        self.fpx_tanh()
    }

    #[inline]
    fn exp_m1(self) -> Self {
        self.fpx_expm1()
    }

    #[inline]
    fn ln_1p(self) -> Self {
        self.fpx_log1p()
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        Self::from_f32(self.to_f32().powi(n))
    }

    #[inline]
    fn powf(self, n: Self) -> Self {
        Self::from_f32(self.to_f32().powf(n.to_f32()))
    }

    #[inline]
    fn to_degrees(self) -> Self {
        Self::from_f32(self.to_f32().to_degrees())
    }

    #[inline]
    fn to_radians(self) -> Self {
        Self::from_f32(self.to_f32().to_radians())
    }
}

#[cfg(feature = "mx")]
macro_rules! impl_mx_real_functions {
    (
        $type:ty,
        $zero:expr, $one:expr, $pi:expr, $frac_pi_2:expr, $epsilon:expr,
        $infinity:expr, $neg_infinity:expr,
        $sin:ident, $cos:ident, $tan:ident,
        $asin:ident, $acos:ident, $atan:ident,
        $exp:ident, $log:ident,
        $sinh:ident, $cosh:ident, $tanh:ident,
        $expm1:ident, $log1p:ident
    ) => {
        impl MxRealFunctions for $type {
            const ZERO: Self = $zero;
            const ONE: Self = $one;
            const PI: Self = $pi;
            const FRAC_PI_2: Self = $frac_pi_2;
            const EPSILON: Self = $epsilon;
            const INFINITY: Self = $infinity;
            const NEG_INFINITY: Self = $neg_infinity;

            #[inline]
            fn fpx_sin(self) -> Self {
                fpx_core::trig::$sin(self)
            }
            #[inline]
            fn fpx_cos(self) -> Self {
                fpx_core::trig::$cos(self)
            }
            #[inline]
            fn fpx_tan(self) -> Self {
                fpx_core::trig::$tan(self)
            }
            #[inline]
            fn fpx_asin(self) -> Self {
                fpx_core::trig::$asin(self)
            }
            #[inline]
            fn fpx_acos(self) -> Self {
                fpx_core::trig::$acos(self)
            }
            #[inline]
            fn fpx_atan(self) -> Self {
                fpx_core::trig::$atan(self)
            }
            #[inline]
            fn fpx_exp(self) -> Self {
                fpx_core::trig::$exp(self)
            }
            #[inline]
            fn fpx_log(self) -> Self {
                fpx_core::trig::$log(self)
            }
            #[inline]
            fn fpx_sinh(self) -> Self {
                fpx_core::trig::$sinh(self)
            }
            #[inline]
            fn fpx_cosh(self) -> Self {
                fpx_core::trig::$cosh(self)
            }
            #[inline]
            fn fpx_tanh(self) -> Self {
                fpx_core::trig::$tanh(self)
            }
            #[inline]
            fn fpx_expm1(self) -> Self {
                fpx_core::trig::$expm1(self)
            }
            #[inline]
            fn fpx_log1p(self) -> Self {
                fpx_core::trig::$log1p(self)
            }
        }
    };
}

#[cfg(feature = "mx")]
impl_mx_real_functions!(
    fpx_core::ocp::fp4::Fp4E2M1,
    fpx_core::ocp::fp4::Fp4E2M1::from_bits(0x0),
    fpx_core::ocp::fp4::Fp4E2M1::from_bits(0x2),
    fpx_core::ocp::fp4::Fp4E2M1::from_bits(0x5),
    fpx_core::ocp::fp4::Fp4E2M1::from_bits(0x3),
    fpx_core::ocp::fp4::Fp4E2M1::from_bits(0x1),
    fpx_core::ocp::fp4::Fp4E2M1::from_bits(0x7),
    fpx_core::ocp::fp4::Fp4E2M1::from_bits(0xf),
    sin_fp4_e2m1,
    cos_fp4_e2m1,
    tan_fp4_e2m1,
    asin_fp4_e2m1,
    acos_fp4_e2m1,
    atan_fp4_e2m1,
    exp_fp4_e2m1,
    log_fp4_e2m1,
    sinh_fp4_e2m1,
    cosh_fp4_e2m1,
    tanh_fp4_e2m1,
    expm1_fp4_e2m1,
    log1p_fp4_e2m1
);

#[cfg(feature = "mx")]
impl_mx_real_functions!(
    fpx_core::ocp::fp6::Fp6E3M2,
    fpx_core::ocp::fp6::Fp6E3M2::from_bits(0x00),
    fpx_core::ocp::fp6::Fp6E3M2::from_bits(0x0c),
    fpx_core::ocp::fp6::Fp6E3M2::from_bits(0x12),
    fpx_core::ocp::fp6::Fp6E3M2::from_bits(0x0e),
    fpx_core::ocp::fp6::Fp6E3M2::from_bits(0x04),
    fpx_core::ocp::fp6::Fp6E3M2::from_bits(0x1f),
    fpx_core::ocp::fp6::Fp6E3M2::from_bits(0x3f),
    sin_fp6_e3m2,
    cos_fp6_e3m2,
    tan_fp6_e3m2,
    asin_fp6_e3m2,
    acos_fp6_e3m2,
    atan_fp6_e3m2,
    exp_fp6_e3m2,
    log_fp6_e3m2,
    sinh_fp6_e3m2,
    cosh_fp6_e3m2,
    tanh_fp6_e3m2,
    expm1_fp6_e3m2,
    log1p_fp6_e3m2
);

#[cfg(feature = "mx")]
impl_mx_real_functions!(
    fpx_core::ocp::fp8::Fp8E4M3,
    fpx_core::ocp::fp8::Fp8E4M3::from_bits(0x00),
    fpx_core::ocp::fp8::Fp8E4M3::from_bits(0x38),
    fpx_core::ocp::fp8::Fp8E4M3::from_bits(0x45),
    fpx_core::ocp::fp8::Fp8E4M3::from_bits(0x3d),
    fpx_core::ocp::fp8::Fp8E4M3::from_bits(0x20),
    fpx_core::ocp::fp8::Fp8E4M3::from_bits(0x7e),
    fpx_core::ocp::fp8::Fp8E4M3::from_bits(0xfe),
    sin_fp8_e4m3,
    cos_fp8_e4m3,
    tan_fp8_e4m3,
    asin_fp8_e4m3,
    acos_fp8_e4m3,
    atan_fp8_e4m3,
    exp_fp8_e4m3,
    log_fp8_e4m3,
    sinh_fp8_e4m3,
    cosh_fp8_e4m3,
    tanh_fp8_e4m3,
    expm1_fp8_e4m3,
    log1p_fp8_e4m3
);

#[cfg(feature = "mx")]
impl_mx_real_functions!(
    fpx_core::ocp::fp8::Fp8E5M2,
    fpx_core::ocp::fp8::Fp8E5M2::from_bits(0x00),
    fpx_core::ocp::fp8::Fp8E5M2::from_bits(0x3c),
    fpx_core::ocp::fp8::Fp8E5M2::from_bits(0x42),
    fpx_core::ocp::fp8::Fp8E5M2::from_bits(0x3e),
    fpx_core::ocp::fp8::Fp8E5M2::from_bits(0x34),
    fpx_core::ocp::fp8::Fp8E5M2::POS_INF,
    fpx_core::ocp::fp8::Fp8E5M2::NEG_INF,
    sin_fp8_e5m2,
    cos_fp8_e5m2,
    tan_fp8_e5m2,
    asin_fp8_e5m2,
    acos_fp8_e5m2,
    atan_fp8_e5m2,
    exp_fp8_e5m2,
    log_fp8_e5m2,
    sinh_fp8_e5m2,
    cosh_fp8_e5m2,
    tanh_fp8_e5m2,
    expm1_fp8_e5m2,
    log1p_fp8_e5m2
);
