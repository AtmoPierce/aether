use crate::core::exp2i::*;
use crate::core::rep::UnpackedE3;
use core::cmp::Ordering;
use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};

/// FP6 E3M2 (1 sign, 3 exponent, 2 mantissa)
///
/// Layout: s e e e m m (6 bits)
/// Bias = 3
///
/// For e = 0 (subnormals / zero):
///   value = sign * (m / 4) * 2^(1 - bias) = sign * (m/4) * 2^-2
///
/// For 1 <= e <= 7 (normals):
///   value = sign * (1 + m/4) * 2^(e - bias)
///
/// No encodings are reserved for Inf or NaN; all patterns are finite.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(transparent)]
pub struct Fp6E3M2(u8);

// ----------------------------------------------------
// Internal unpacked representation (finite only)
// value = sign * mant * 2^exp, mant > 0 for non-zero
// ----------------------------------------------------

impl Fp6E3M2 {
    /// Max finite positive: S 111 11_2 = + 2^4 × 1.75 = +28.0
    pub const MAX_FINITE: Self = Self(0b0_111_11);

    /// Min normal positive: S 001 00_2 = + 2^-2 × 1.0 = +0.25
    pub const MIN_NORMAL_POS: Self = Self(0b0_001_00);

    /// Max subnormal positive: S 000 11_2 = + 2^-2 × 0.75 = +0.1875
    pub const MAX_SUBNORMAL_POS: Self = Self(0b0_000_11);

    /// Min subnormal positive: S 000 01_2 = + 2^-2 × 0.25 = +0.0625
    pub const MIN_SUBNORMAL_POS: Self = Self(0b0_000_01);

    #[inline]
    pub const fn from_bits(bits: u8) -> Self {
        Self(bits & 0x3F)
    }

    #[inline]
    pub const fn to_bits(self) -> u8 {
        self.0 & 0x3F
    }

    /// Saturating, round-to-nearest-ties-to-even conversion from f32.
    ///
    /// - Values beyond FP6 range clamp to +-MAX_FINITE.
    /// - Magnitudes below min subnormal after rounding convert to +-0.0.
    /// - NaN is mapped to 0.0 (implementation-defined per spec).
    #[inline]
    pub fn quantize(v: f32) -> Self {
        if v.is_nan() {
            return Self(0);
        }

        let mut best_bits: u8 = 0;
        let mut best_err: f32 = f32::INFINITY;

        for bits in 0u8..64 {
            let decoded = fp6e3m2_decode(bits);
            let err = (decoded - v).abs();

            if err < best_err {
                best_err = err;
                best_bits = bits;
            } else if err == best_err {
                // prefer LSB=0
                let cur_lsb = best_bits & 0x1;
                let new_lsb = bits & 0x1;
                if new_lsb == 0 && cur_lsb == 1 {
                    best_bits = bits;
                }
            }
        }

        Self(best_bits)
    }

    // -------------------------
    // unpack / pack (integer)
    // -------------------------

    /// Unpack bits into integer representation:
    /// value = sign * mant * 2^exp
    #[inline]
    fn unpack(self) -> UnpackedE3 {
        let bits = self.to_bits();
        let s = (bits >> 5) & 0x1;
        let e = (bits >> FP6_E3M2_MAN_BITS) & 0x07;
        let m = bits & 0x03;

        let sign = if s == 0 { 1 } else { -1 };

        // zero
        if e == 0 && m == 0 {
            return UnpackedE3 {
                sign: 1,
                exp: 0,
                mant: 0,
                is_zero: true,
                is_subnormal: false,
            };
        }

        if e == 0 {
            // subnormals:
            // value = sign * (m/4) * 2^(1-bias)
            //       = sign * m * 2^(1-bias-2)
            let exp = (1 - FP6_E3M2_BIAS - FP6_E3M2_MAN_BITS as i32) as i16; // -4
            let mant = m as u16;
            return UnpackedE3 {
                sign,
                exp,
                mant,
                is_zero: false,
                is_subnormal: true,
            };
        }

        // normals:
        // value = sign * (1 + m/4) * 2^(e-bias)
        //       = sign * (4+m) * 2^(e-bias-2)
        let exp = (e as i32 - FP6_E3M2_BIAS - FP6_E3M2_MAN_BITS as i32) as i16;
        let mant = (4 + m) as u16;

        UnpackedE3 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
        }
    }

    /// Pack integer representation back into Fp6E3M2,
    /// saturate in [MIN_SUBNORMAL, MAX_FINITE].
    fn pack_unpacked(u: UnpackedE3) -> Self {
        if u.is_zero || u.mant == 0 {
            return Fp6E3M2::from_bits(0);
        }

        let mut sign = u.sign;
        let mut exp = u.exp;
        let mut mant = u.mant as i32;

        if mant < 0 {
            mant = -mant;
            sign = -sign;
        }

        if mant == 0 {
            return Fp6E3M2::from_bits(0);
        }

        // Normalize mant so it fits in [1,7]
        while mant > 7 {
            // rounded shift right by 1
            mant = if mant & 1 != 0 {
                (mant + 1) >> 1
            } else {
                mant >> 1
            };
            exp += 1;
        }

        // normal (mant >= 4)
        // -4 (subnormal floor).
        while mant < 4 && mant > 0 && exp > -4 {
            mant <<= 1;
            exp -= 1;
        }

        // Exponent range for representable values:
        // from the table for MX definition exp in [-4, 2].
        if exp < -4 {
            // underflow to zero
            return Fp6E3M2::from_bits(0);
        }
        if exp > 2 {
            // overflow to max finite
            let s_bit = if sign < 0 { 1u8 } else { 0u8 };
            let bits = (s_bit << 5) | (0b111 << 2) | 0b11; // (+/-) 0_111_11
            return Fp6E3M2::from_bits(bits);
        }

        if mant <= 0 {
            return Fp6E3M2::from_bits(0);
        }
        if mant > 7 {
            mant = 7;
        }

        let s_bit = if sign < 0 { 1u8 } else { 0u8 };

        if mant >= 4 {
            // normal
            let e = (exp + FP6_E3M2_BIAS as i16 + FP6_E3M2_MAN_BITS as i16) as u8; // 1..7
            let m_bits = (mant as u8).saturating_sub(4) & 0x03;
            let bits = (s_bit << 5) | ((e & 0x07) << 2) | m_bits;
            return Fp6E3M2::from_bits(bits);
        } else {
            // subnormal: exp must be -4 here, mant in 1..3
            if mant < 1 {
                return Fp6E3M2::from_bits(0);
            }
            let m_bits = (mant as u8).min(3);
            let bits = (s_bit << 5) | m_bits;
            return Fp6E3M2::from_bits(bits);
        }
    }
}

// -------- spec constants --------

const FP6_E3M2_EXP_BITS: u32 = 3;
const FP6_E3M2_MAN_BITS: u32 = 2;
const FP6_E3M2_BIAS: i32 = 3;

/// Decode FP6 E3M2 -> f32 (per spec table).
#[inline]
fn fp6e3m2_decode(bits: u8) -> f32 {
    let bits = bits & 0x3F;

    let s = (bits >> 5) & 0x1;
    let e = (bits >> FP6_E3M2_MAN_BITS) & 0x07; // 3 exponent bits
    let m = bits & 0x03; // 2 mantissa bits

    let sign = if s == 0 { 1.0f32 } else { -1.0f32 };

    if e == 0 {
        // zero or subnormal
        if m == 0 {
            return sign * 0.0;
        }

        // subnormal: value = sign * (m/4) * 2^(1 - bias) = sign * (m/4)*2^-2
        let exp = 1 - FP6_E3M2_BIAS; // = -2
        let frac = (m as f32) / (1u32 << FP6_E3M2_MAN_BITS) as f32;
        return sign * exp2i(exp) * frac;
    }

    // normal: value = sign * (1 + m/4) * 2^(e - bias)
    let exp = (e as i32 - FP6_E3M2_BIAS) as i32;
    let frac = 1.0 + (m as f32) / (1u32 << FP6_E3M2_MAN_BITS) as f32;

    sign * exp2i(exp) * frac
}

// -------- conversions --------

impl From<Fp6E3M2> for f32 {
    #[inline]
    fn from(x: Fp6E3M2) -> Self {
        fp6e3m2_decode(x.to_bits())
    }
}

impl From<f32> for Fp6E3M2 {
    #[inline]
    fn from(v: f32) -> Self {
        Fp6E3M2::quantize(v)
    }
}

// -------- Display --------

impl fmt::Display for Fp6E3M2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let bits = self.to_bits();
        let val: f32 = (*self).into();
        write!(f, "Fp6E3M2({:.8e}, bits=0x{:02X})", val, bits)
    }
}

// -------- Core ops via UnpackedE3 (finite-only) --------

impl Neg for Fp6E3M2 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let u = self.unpack();
        if u.is_zero {
            self
        } else {
            Fp6E3M2::pack_unpacked(UnpackedE3 { sign: -u.sign, ..u })
        }
    }
}

impl Add for Fp6E3M2 {
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

        // align exponents
        let mut exp = a.exp;
        let (mut ma, mut mb) = (a.mant as i32, b.mant as i32);
        let (sa, sb) = (a.sign as i32, b.sign as i32);

        if a.exp > b.exp {
            let shift = (a.exp - b.exp) as u32;
            if shift != 0 {
                mb = (mb + (1 << (shift - 1)) as i32) >> shift;
            }
        } else if a.exp < b.exp {
            let shift = (b.exp - a.exp) as u32;
            if shift != 0 {
                ma = (ma + (1 << (shift - 1)) as i32) >> shift;
            }
            exp = b.exp;
        }

        let va = sa * ma;
        let vb = sb * mb;
        let sum = va + vb;

        if sum == 0 {
            return Fp6E3M2::from_bits(0);
        }

        let sign = if sum < 0 { -1 } else { 1 };
        let mant = sum.abs() as u16;

        Fp6E3M2::pack_unpacked(UnpackedE3 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
        })
    }
}

impl Sub for Fp6E3M2 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl Mul for Fp6E3M2 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let a = self.unpack();
        let b = rhs.unpack();

        if a.is_zero || b.is_zero {
            return Fp6E3M2::from_bits(0);
        }

        let sign = a.sign * b.sign;
        let exp = a.exp + b.exp;
        let mant = (a.mant as i32 * b.mant as i32) as u16;

        Fp6E3M2::pack_unpacked(UnpackedE3 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
        })
    }
}

impl Div for Fp6E3M2 {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let a = self.unpack();
        let b = rhs.unpack();

        if a.is_zero && b.is_zero {
            // 0/0: implementation-defined; just return 0 -> if cast and divided later this would be caught
            return Fp6E3M2::from_bits(0);
        }

        if b.is_zero {
            // finite / 0: saturate to max magnitude
            let sign = a.sign * b.sign.max(1); // if b sign is 0 for zero, treat as +1
            let s_bit = if sign < 0 { 1u8 } else { 0u8 };
            let bits = (s_bit << 5) | (0b111 << 2) | 0b11;
            return Fp6E3M2::from_bits(bits);
        }

        if a.is_zero {
            // 0 / finite = 0
            return Fp6E3M2::from_bits(0);
        }

        let sign = a.sign * b.sign;
        let exp = a.exp - b.exp;

        // mant_res ~= (a.mant << 2) / b.mant for a bit of extra resolution
        let num = (a.mant as i32) << 2;
        let den = b.mant as i32;
        let mant_i = if num >= 0 {
            (num + den / 2) / den
        } else {
            (num - den / 2) / den
        };

        if mant_i <= 0 {
            return Fp6E3M2::from_bits(0);
        }

        let mant = mant_i as u16;

        Fp6E3M2::pack_unpacked(UnpackedE3 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
        })
    }
}

impl PartialOrd for Fp6E3M2 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let a = self.unpack();
        let b = other.unpack();

        // no NaNs in encoding -> always comparable
        if a.is_zero && b.is_zero {
            return Some(Ordering::Equal);
        }

        if a.sign != b.sign {
            return Some(if a.sign < b.sign {
                Ordering::Less
            } else {
                Ordering::Greater
            });
        }

        // same sign
        if a.exp != b.exp {
            if a.exp < b.exp {
                return Some(if a.sign > 0 {
                    Ordering::Less
                } else {
                    Ordering::Greater
                });
            } else {
                return Some(if a.sign > 0 {
                    Ordering::Greater
                } else {
                    Ordering::Less
                });
            }
        }

        if a.mant != b.mant {
            if a.mant < b.mant {
                return Some(if a.sign > 0 {
                    Ordering::Less
                } else {
                    Ordering::Greater
                });
            } else {
                return Some(if a.sign > 0 {
                    Ordering::Greater
                } else {
                    Ordering::Less
                });
            }
        }

        Some(Ordering::Equal)
    }
}
