#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct E8M0 {
    /// Exponent in [-127, 127], with -128 as NaN sentinel.
    value: i8,
}

impl E8M0 {
    pub const NAN: Self = Self { value: i8::MIN }; // -128

    /// Construct from an exponent in [-127, 127].
    #[inline]
    pub fn new(exp: i8) -> Self {
        assert!(
            (-127..=127).contains(&exp),
            "E8M0 exponent out of range (-127..=127)"
        );
        Self { value: exp }
    }

    /// Construct from raw bits, treating -128 as NaN.
    #[inline]
    pub fn from_bits(bits: i8) -> Self {
        if bits == i8::MIN {
            Self::NAN
        } else {
            Self { value: bits }
        }
    }

    /// Raw stored value (including -128 for NaN, i.e. no nan check).
    #[inline]
    pub const fn bits(self) -> i8 {
        self.value
    }

    // Self-explanatory
    #[inline]
    pub const fn is_nan(self) -> bool {
        self.value == i8::MIN
    }

    /// Get the exponent as i8 if not NaN.
    #[inline]
    pub const fn exponent(self) -> Option<i8> {
        if self.is_nan() {
            None
        } else {
            Some(self.value)
        }
    }

    /// Power-of-two scale factor as f32, or NaN if E8M0 is NaN.
    #[inline]
    pub fn scale_f32(self) -> f32 {
        match self.exponent() {
            Some(e) => 2f32.powi(e as i32),
            None => f32::NAN,
        }
    }
}

// ===============================
// Display impls
// ===============================
use core::fmt;
impl fmt::Display for E8M0 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let bits = self.bits();
        if self.is_nan() {
            write!(f, "E8M0(NaN, bits={})", bits)
        } else if let Some(e) = self.exponent() {
            let s = self.scale_f32();
            write!(f, "E8M0(exp={}, scale={:.8e}, bits={})", e, s, bits)
        } else {
            // shouldn't happen, but be defensive
            write!(f, "E8M0(Invalid, bits={})", bits)
        }
    }
}
