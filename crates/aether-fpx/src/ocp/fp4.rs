use crate::ocp::exp2i::exp2i;
use crate::ocp::rep::UnpackedE2;

const FP4_SIGN_MASK: u8 = 0b1000;
const FP4_INDEX_MASK: u8 = 0b0111;

/// Magnitude table for FP4 "E2M1-like":
///
/// index -> |value|
/// 0 -> 0.0
/// 1 -> 0.5
/// 2 -> 1.0
/// 3 -> 1.5
/// 4 -> 2.0
/// 5 -> 3.0
/// 6 -> 4.0
/// 7 -> 6.0
const FP4_MAGS: [f32; 8] = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0];

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(transparent)]
pub struct Fp4E2M1(u8);

impl Fp4E2M1 {
    #[inline]
    pub const fn from_bits(bits: u8) -> Self {
        Self(bits & 0x0F)
    }

    #[inline]
    pub const fn to_bits(self) -> u8 {
        self.0 & 0x0F
    }

    /// Max finite representable magnitude as f32 (|6.0|).
    #[inline]
    pub fn max_finite_f32() -> f32 {
        6.0
    }

    /// Decode into `UnpackedE2` for integer-ish arithmetic.
    ///
    /// Mapping (positive magnitudes):
    ///   0.5 -> mant=2, exp=-2
    ///   1.0 -> mant=2, exp=-1
    ///   1.5 -> mant=3, exp=-1
    ///   2.0 -> mant=2, exp= 0
    ///   3.0 -> mant=3, exp= 0
    ///   4.0 -> mant=2, exp= 1
    ///   6.0 -> mant=3, exp= 1
    ///
    /// So that value = sign * mant * 2^exp.
    #[inline]
    fn unpack(self) -> UnpackedE2 {
        let bits = self.to_bits();
        let sign_bit = (bits & FP4_SIGN_MASK) >> 3;
        let idx = bits & FP4_INDEX_MASK;

        // zero (+0 or -0 both treated as is_zero)
        if idx == 0 {
            return UnpackedE2 {
                sign: 1,
                exp: 0,
                mant: 0,
                is_zero: true,
            };
        }

        let sign = if sign_bit == 0 { 1 } else { -1 };

        let (mant, exp) = match idx {
            1 => (2u8, -2i8), // 0.5
            2 => (2u8, -1i8), // 1.0
            3 => (3u8, -1i8), // 1.5
            4 => (2u8, 0i8),  // 2.0
            5 => (3u8, 0i8),  // 3.0
            6 => (2u8, 1i8),  // 4.0
            7 => (3u8, 1i8),  // 6.0
            _ => unreachable!("idx is 0..7 after masking"),
        };

        UnpackedE2 {
            sign,
            exp,
            mant,
            is_zero: false,
        }
    }

    /// Pack an `UnpackedE2` back to FP4 by rounding to nearest representable
    /// magnitude in `FP4_MAGS` and applying the sign.
    #[inline]
    fn pack(u: UnpackedE2) -> Self {
        // Zero / underflow
        if u.is_zero || u.mant == 0 {
            return Fp4E2M1::from_bits(0x0); // +0.0
        }

        let sign_bit = if u.sign < 0 { 1u8 } else { 0u8 };

        // Compute approximate positive magnitude
        let mag = (u.mant as f32) * exp2i(u.exp as i32).abs();

        // Find nearest entry in FP4_MAGS (including 0.0 for possible underflow).
        let mut best_idx: u8 = 0;
        let mut best_err: f32 = f32::INFINITY;

        for idx in 0u8..8 {
            let target = FP4_MAGS[idx as usize];
            let err = (mag - target).abs();
            if err < best_err {
                best_err = err;
                best_idx = idx;
            }
        }

        // Construct bits: sign in bit 3, index in bits 2..0
        let bits = (sign_bit << 3) | (best_idx & 0x07);
        Fp4E2M1::from_bits(bits)
    }

    /// Quantize from f32 -> Fp4E2M1 by rounding to nearest magnitude, with sign.
    pub fn quantize(v: f32) -> Self {
        if v == 0.0 || !v.is_finite() {
            // No NaN/Inf encoding; treat non-finite as 0 for FP4.
            return Fp4E2M1::from_bits(0x0);
        }

        let sign_bit = if v.is_sign_negative() { 1u8 } else { 0u8 };
        let mag = v.abs();

        // Nearest neighbor among magnitudes 0.5..6.0 (indices 1..7).
        let mut best_idx: u8 = 1;
        let mut best_err: f32 = f32::INFINITY;

        for idx in 1u8..8 {
            let target = FP4_MAGS[idx as usize];
            let err = (mag - target).abs();
            if err < best_err {
                best_err = err;
                best_idx = idx;
            }
        }

        let bits = (sign_bit << 3) | (best_idx & 0x07);
        Fp4E2M1::from_bits(bits)
    }
}

/// Decode bits -> f32 according to the FP4 table.
///
/// Uses FP4_MAGS and a sign bit. 0x0 -> +0.0, 0x8 -> -0.0.
#[inline]
fn fp4e2m1_decode(bits: u8) -> f32 {
    let bits = bits & 0x0F;
    let sign_bit = (bits & FP4_SIGN_MASK) >> 3;
    let idx = bits & FP4_INDEX_MASK;

    let mag = FP4_MAGS[idx as usize];

    if mag == 0.0 {
        // Preserve sign of zero
        if sign_bit == 0 {
            0.0
        } else {
            -0.0
        }
    } else {
        let sign = if sign_bit == 0 { 1.0f32 } else { -1.0f32 };
        sign * mag
    }
}

impl From<Fp4E2M1> for f32 {
    #[inline]
    fn from(x: Fp4E2M1) -> Self {
        fp4e2m1_decode(x.to_bits())
    }
}

impl From<f32> for Fp4E2M1 {
    #[inline]
    fn from(v: f32) -> Self {
        Self::quantize(v)
    }
}

/////////////////////////////////////
///         Operations
/////////////////////////////////////
use core::cmp::Ordering;
use core::ops::{Add, Div, Mul, Neg, Sub};

impl Neg for Fp4E2M1 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let u = self.unpack();
        if u.is_zero {
            // Preserve signed zero through decode/encode
            let bits = self.to_bits() ^ FP4_SIGN_MASK;
            Fp4E2M1::from_bits(bits)
        } else {
            Fp4E2M1::pack(UnpackedE2 { sign: -u.sign, ..u })
        }
    }
}

impl Add for Fp4E2M1 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let a = self.unpack();
        let b = rhs.unpack();

        if a.is_zero {
            return rhs;
        }
        if b.is_zero {
            return self;
        }

        // Align exponents
        let mut exp = a.exp;
        let (mut ma, mut mb) = (a.mant as i16, b.mant as i16);
        let (sa, sb) = (a.sign as i16, b.sign as i16);

        if a.exp > b.exp {
            let shift = (a.exp - b.exp) as u32;
            if shift != 0 {
                mb = (mb + (1 << (shift - 1)) as i16) >> shift;
            }
        } else if a.exp < b.exp {
            let shift = (b.exp - a.exp) as u32;
            if shift != 0 {
                ma = (ma + (1 << (shift - 1)) as i16) >> shift;
            }
            exp = b.exp;
        }

        let va = sa * ma;
        let vb = sb * mb;
        let sum = va + vb;

        if sum == 0 {
            return Fp4E2M1::from_bits(0x0);
        }

        let sign = if sum < 0 { -1 } else { 1 };
        let mant = sum.abs() as u8;

        Fp4E2M1::pack(UnpackedE2 {
            sign,
            exp,
            mant,
            is_zero: false,
        })
    }
}

impl Sub for Fp4E2M1 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl Mul for Fp4E2M1 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let a = self.unpack();
        let b = rhs.unpack();

        if a.is_zero || b.is_zero {
            return Fp4E2M1::from_bits(0x0);
        }

        let sign = a.sign * b.sign;
        let exp = a.exp + b.exp;
        let mant = (a.mant as i16 * b.mant as i16) as u8;

        Fp4E2M1::pack(UnpackedE2 {
            sign,
            exp,
            mant,
            is_zero: false,
        })
    }
}

impl Div for Fp4E2M1 {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let a = self.unpack();
        let b = rhs.unpack();

        // 0 / 0 or x / 0 -> treat as 0; decode path has no NaN encoding.
        if b.is_zero {
            return Fp4E2M1::from_bits(0x0);
        }

        if a.is_zero {
            return Fp4E2M1::from_bits(0x0);
        }

        let sign = a.sign * b.sign;
        // `mant_i` is scaled by 2 to retain one fraction bit. Compensate in
        // the integer exponent so the represented quotient keeps its scale.
        let exp = a.exp - b.exp - 1;

        // Approx mant_res ~= 2 * Ma / Mb
        let num = (a.mant as i16) * 2; // 2 * Ma
        let den = b.mant as i16;
        let mant_i = (num + den / 2) / den; // round-to-nearest

        if mant_i <= 0 {
            return Fp4E2M1::from_bits(0x0);
        }

        let mant = mant_i as u8;

        Fp4E2M1::pack(UnpackedE2 {
            sign,
            exp,
            mant,
            is_zero: false,
        })
    }
}

impl PartialOrd for Fp4E2M1 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let a = self.unpack();
        let b = other.unpack();

        if a.is_zero && b.is_zero {
            return Some(Ordering::Equal);
        }
        if a.is_zero {
            return Some(if b.sign < 0 {
                Ordering::Greater
            } else {
                Ordering::Less
            });
        }
        if b.is_zero {
            return Some(if a.sign < 0 {
                Ordering::Less
            } else {
                Ordering::Greater
            });
        }

        if a.sign != b.sign {
            return if a.sign < b.sign {
                Some(Ordering::Less)
            } else {
                Some(Ordering::Greater)
            };
        }

        if a.exp != b.exp {
            return if a.exp < b.exp {
                Some(if a.sign > 0 {
                    Ordering::Less
                } else {
                    Ordering::Greater
                })
            } else {
                Some(if a.sign > 0 {
                    Ordering::Greater
                } else {
                    Ordering::Less
                })
            };
        }

        if a.mant != b.mant {
            return if a.mant < b.mant {
                Some(if a.sign > 0 {
                    Ordering::Less
                } else {
                    Ordering::Greater
                })
            } else {
                Some(if a.sign > 0 {
                    Ordering::Greater
                } else {
                    Ordering::Less
                })
            };
        }

        Some(Ordering::Equal)
    }
}

impl core::fmt::Display for Fp4E2M1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let bits = self.to_bits();
        let v: f32 = (*self).into();
        write!(f, "Fp4E2M1({:.8e}, bits=0x{:01X})", v, bits)
    }
}
