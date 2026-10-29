use core::cmp::Ordering;
use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::ocp::exp2i::*;
use crate::ocp::rep::{UnpackedE4, UnpackedE5};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(transparent)]
pub struct Fp8E5M2(u8);

impl Fp8E5M2 {
    // ---- constants ----

    /// Canonical NaN: s=0, e=31, m=3 => 0b0_11111_11 = 0x7F
    pub const CANONICAL_NAN: Self = Self(0x7F);

    /// +Infinity: s=0, e=31, m=0 => 0b0_11111_00 = 0x7C
    pub const POS_INF: Self = Self(0x7C);

    /// -Infinity: s=1, e=31, m=0 => 0b1_11111_00 = 0xFC
    pub const NEG_INF: Self = Self(0xFC);

    /// Max finite positive: s=0, e=30, m=3 => 0b0_11110_11 = 0x7B
    pub const MAX_FINITE: Self = Self(0x7B);

    /// Smallest positive normal: s=0, e=1, m=0 => 0b0_00001_00 = 0x04
    pub const MIN_NORMAL_POS: Self = Self(0x04);

    /// Smallest positive subnormal: s=0, e=0, m=1 => 0b0_00000_01 = 0x01
    pub const MIN_SUBNORMAL_POS: Self = Self(0x01);

    #[inline]
    pub const fn from_bits(bits: u8) -> Self {
        Self(bits)
    }

    #[inline]
    pub const fn to_bits(self) -> u8 {
        self.0
    }

    /// Saturates to nearest-neighbor from f32 to FP8 E5M2.
    #[inline]
    pub fn quantize(v: f32) -> Self {
        if v.is_nan() {
            return Self::CANONICAL_NAN;
        }

        if v.is_infinite() {
            return if v.is_sign_positive() {
                Self::POS_INF
            } else {
                Self::NEG_INF
            };
        }

        let mut best_bits: u8 = 0;
        let mut best_err: f32 = f32::INFINITY;

        for bits in 0u8..=255 {
            let decoded = fp8e5m2_decode(bits);

            // skip NaNs (and keep Inf reserved)
            if decoded.is_nan() || decoded.is_infinite() {
                continue;
            }

            let err = (decoded - v).abs();
            if err < best_err {
                best_err = err;
                best_bits = bits;
            } else if err == best_err {
                // tie-break: prefer even mantissa LSB
                let cur_lsb = best_bits & 0x1;
                let new_lsb = bits & 0x1;
                if new_lsb == 0 && cur_lsb == 1 {
                    best_bits = bits;
                }
            }
        }

        Self(best_bits)
    }

    /// Unpack into integer representation:
    ///
    /// For normals (1 <= e <= 30):
    ///   value = sign * (1 + m/4) * 2^(e - BIAS)
    ///         = sign * (4 + m) * 2^(e - BIAS - 2)
    ///
    /// For subnormals (e = 0, m != 0):
    ///   value = sign * (m/4) * 2^(1 - BIAS)
    ///         = sign * m * 2^(1 - BIAS - 2)
    #[inline]
    fn unpack(self) -> UnpackedE5 {
        let bits = self.to_bits();
        let s = (bits >> 7) & 0x1;
        let e = (bits >> 2) & 0x1F;
        let m = bits & 0x03;

        // NaN / Inf
        if e == 0x1F {
            if m == 0 {
                // Inf
                let sign = if s == 0 { 1 } else { -1 };
                return UnpackedE5 {
                    sign,
                    exp: 0,
                    mant: 0,
                    is_zero: false,
                    is_subnormal: false,
                    is_inf: true,
                    is_nan: false,
                };
            } else {
                // NaN
                return UnpackedE5 {
                    sign: 1,
                    exp: 0,
                    mant: 0,
                    is_zero: false,
                    is_subnormal: false,
                    is_inf: false,
                    is_nan: true,
                };
            }
        }

        // Zero / subnormal / normal
        if e == 0 {
            if m == 0 {
                // Zero
                return UnpackedE5 {
                    sign: 1,
                    exp: 0,
                    mant: 0,
                    is_zero: true,
                    is_subnormal: false,
                    is_inf: false,
                    is_nan: false,
                };
            } else {
                // Subnormal
                let sign = if s == 0 { 1 } else { -1 };
                let exp = (1 - FP8_E5M2_BIAS - FP8_E5M2_MAN_BITS as i32) as i16; // -16
                let mant = m as u16;
                return UnpackedE5 {
                    sign,
                    exp,
                    mant,
                    is_zero: false,
                    is_subnormal: true,
                    is_inf: false,
                    is_nan: false,
                };
            }
        }

        // Normal
        let sign = if s == 0 { 1 } else { -1 };
        let exp = (e as i32 - FP8_E5M2_BIAS - FP8_E5M2_MAN_BITS as i32) as i16;
        let mant = (4 + m) as u16; // 4..7

        UnpackedE5 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
            is_inf: false,
            is_nan: false,
        }
    }

    /// Pack integer representation back into Fp8E5M2.
    ///
    /// - Propagates NaN.
    /// - Handles Inf.
    /// - Normalizes mant and exp into encodable range.
    fn pack_unpacked(u: UnpackedE5) -> Self {
        if u.is_nan {
            return Fp8E5M2::CANONICAL_NAN;
        }

        if u.is_inf {
            return if u.sign < 0 {
                Fp8E5M2::NEG_INF
            } else {
                Fp8E5M2::POS_INF
            };
        }

        if u.is_zero || u.mant == 0 {
            return Fp8E5M2::from_bits(0);
        }

        let mut sign = u.sign;
        let mut exp = u.exp;
        let mut mant = u.mant as i32;

        // Normalize mant into ~[4,7].
        while mant > 7 {
            mant = if mant & 1 != 0 {
                (mant + 1) >> 1
            } else {
                mant >> 1
            };
            exp += 1;
        }

        while mant < 4 && mant > 0 && exp > -16 {
            mant <<= 1;
            exp -= 1;
        }

        // finite exponent range for normals in exp-int space is roughly [-16, 13].
        if exp > 13 {
            // overflow -> Inf
            return if sign < 0 {
                Fp8E5M2::NEG_INF
            } else {
                Fp8E5M2::POS_INF
            };
        }
        if exp < -16 {
            // underflow -> zero
            return Fp8E5M2::from_bits(0);
        }

        if mant <= 0 {
            return Fp8E5M2::from_bits(0);
        }
        if mant > 7 {
            mant = 7;
        }

        // convert back: e = exp + bias + MAN_BITS
        let e = (exp + FP8_E5M2_BIAS as i16 + FP8_E5M2_MAN_BITS as i16) as u8; // 1..30
        if e == 0 || e >= 0x1F {
            // overflowed into Inf/NaN space -> Inf
            return if sign < 0 {
                Fp8E5M2::NEG_INF
            } else {
                Fp8E5M2::POS_INF
            };
        }

        let s_bit = if sign < 0 { 1u8 } else { 0u8 };
        let m_bits = (mant as u8).saturating_sub(4) & 0x03;

        let bits = (s_bit << 7) | ((e & 0x1F) << 2) | m_bits;
        Fp8E5M2::from_bits(bits)
    }
}

// internal helper constants
const FP8_E5M2_EXP_BITS: u32 = 5;
const FP8_E5M2_MAN_BITS: u32 = 2;
const FP8_E5M2_BIAS: i32 = 15;

#[inline]
fn fp8e5m2_decode(bits: u8) -> f32 {
    let s = (bits >> 7) & 0x1;
    let e = (bits >> 2) & 0x1F;
    let m = bits & 0x03;

    let sign = if s == 0 { 1.0f32 } else { -1.0f32 };

    if e == 0 {
        // Zero or subnormal
        if m == 0 {
            return sign * 0.0;
        }

        // Subnormal: exponent = 1 - bias
        let exp = 1 - FP8_E5M2_BIAS;
        let frac = (m as f32) / (1u32 << FP8_E5M2_MAN_BITS) as f32;
        return sign * exp2i(exp) * frac;
    }

    if e == 0x1F {
        // Specials: Inf / NaN
        if m == 0 {
            return sign * f32::INFINITY;
        } else {
            return f32::NAN;
        }
    }

    // Normal
    let exp = (e as i32 - FP8_E5M2_BIAS) as i32;
    let frac = 1.0 + (m as f32) / (1u32 << FP8_E5M2_MAN_BITS) as f32;
    sign * exp2i(exp) * frac
}

impl From<Fp8E5M2> for f32 {
    #[inline]
    fn from(x: Fp8E5M2) -> Self {
        fp8e5m2_decode(x.0)
    }
}

impl From<f32> for Fp8E5M2 {
    #[inline]
    fn from(v: f32) -> Self {
        Fp8E5M2::quantize(v)
    }
}

// Display for E5M2
impl fmt::Display for Fp8E5M2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let bits = self.to_bits();
        let val: f32 = (*self).into();

        if bits == Self::CANONICAL_NAN.to_bits() || val.is_nan() {
            write!(f, "Fp8E5M2(NaN, bits=0x{:02X})", bits)
        } else if bits == Self::POS_INF.to_bits() {
            write!(f, "Fp8E5M2(+inf, bits=0x{:02X})", bits)
        } else if bits == Self::NEG_INF.to_bits() {
            write!(f, "Fp8E5M2(-inf, bits=0x{:02X})", bits)
        } else {
            write!(f, "Fp8E5M2({:.8e}, bits=0x{:02X})", val, bits)
        }
    }
}

// ---------- arithmetic for Fp8E5M2 via UnpackedE5 ----------

impl Neg for Fp8E5M2 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let u = self.unpack();

        if u.is_nan {
            Fp8E5M2::CANONICAL_NAN
        } else if u.is_inf {
            Fp8E5M2::pack_unpacked(UnpackedE5 { sign: -u.sign, ..u })
        } else if u.is_zero {
            self
        } else {
            Fp8E5M2::pack_unpacked(UnpackedE5 { sign: -u.sign, ..u })
        }
    }
}

impl Add for Fp8E5M2 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let a = self.unpack();
        let b = rhs.unpack();

        if a.is_nan || b.is_nan {
            return Fp8E5M2::CANONICAL_NAN;
        }

        // Inf rules
        if a.is_inf || b.is_inf {
            if a.is_inf && b.is_inf && a.sign != b.sign {
                return Fp8E5M2::CANONICAL_NAN;
            }
            return if a.is_inf { self } else { rhs };
        }

        if a.is_zero {
            return rhs;
        }
        if b.is_zero {
            return self;
        }

        // Align exponents
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
            return Fp8E5M2::from_bits(0);
        }

        let sign = if sum < 0 { -1 } else { 1 };
        let mant = sum.abs() as u16;

        Fp8E5M2::pack_unpacked(UnpackedE5 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
            is_inf: false,
            is_nan: false,
        })
    }
}

impl Sub for Fp8E5M2 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl Mul for Fp8E5M2 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let a = self.unpack();
        let b = rhs.unpack();

        if a.is_nan || b.is_nan {
            return Fp8E5M2::CANONICAL_NAN;
        }

        // Inf rules
        if a.is_inf || b.is_inf {
            if (a.is_inf && b.is_zero) || (b.is_inf && a.is_zero) {
                return Fp8E5M2::CANONICAL_NAN;
            }

            let sign = a.sign * b.sign;
            return Fp8E5M2::pack_unpacked(UnpackedE5 {
                sign,
                exp: 0,
                mant: 0,
                is_zero: false,
                is_subnormal: false,
                is_inf: true,
                is_nan: false,
            });
        }

        if a.is_zero || b.is_zero {
            return Fp8E5M2::from_bits(0);
        }

        let sign = a.sign * b.sign;
        let exp = a.exp + b.exp;
        let mant = (a.mant as i32 * b.mant as i32) as u16;

        Fp8E5M2::pack_unpacked(UnpackedE5 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
            is_inf: false,
            is_nan: false,
        })
    }
}

impl Div for Fp8E5M2 {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let a = self.unpack();
        let b = rhs.unpack();

        if a.is_nan || b.is_nan {
            return Fp8E5M2::CANONICAL_NAN;
        }

        // Inf / Inf => NaN
        if a.is_inf && b.is_inf {
            return Fp8E5M2::CANONICAL_NAN;
        }

        // Inf / finite -> Inf
        if a.is_inf && !b.is_zero {
            let sign = a.sign * b.sign;
            return Fp8E5M2::pack_unpacked(UnpackedE5 {
                sign,
                exp: 0,
                mant: 0,
                is_zero: false,
                is_subnormal: false,
                is_inf: true,
                is_nan: false,
            });
        }

        // finite / Inf -> 0
        if b.is_inf && !a.is_zero {
            return Fp8E5M2::from_bits(0);
        }

        // finite / 0 -> NaN
        if b.is_zero {
            return Fp8E5M2::CANONICAL_NAN;
        }

        // 0 / finite -> 0
        if a.is_zero {
            return Fp8E5M2::from_bits(0);
        }

        let sign = a.sign * b.sign;
        // `mant_i` is scaled by 2^2 to retain two fraction bits.
        let exp = a.exp - b.exp - FP8_E5M2_MAN_BITS as i16;

        // mant_res ~= (a.mant << 2) / b.mant
        let num = (a.mant as i32) << 2;
        let den = b.mant as i32;

        let mant_i = if num >= 0 {
            (num + den / 2) / den
        } else {
            (num - den / 2) / den
        };

        if mant_i <= 0 {
            return Fp8E5M2::from_bits(0);
        }

        let mant = mant_i as u16;

        Fp8E5M2::pack_unpacked(UnpackedE5 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
            is_inf: false,
            is_nan: false,
        })
    }
}

impl PartialOrd for Fp8E5M2 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let a = self.unpack();
        let b = other.unpack();

        if a.is_nan || b.is_nan {
            return None;
        }

        // Inf
        if a.is_inf || b.is_inf {
            if a.is_inf && b.is_inf {
                if a.sign == b.sign {
                    return Some(Ordering::Equal);
                } else if a.sign < b.sign {
                    return Some(Ordering::Less);
                } else {
                    return Some(Ordering::Greater);
                }
            }
            if a.is_inf {
                return Some(if a.sign > 0 {
                    Ordering::Greater
                } else {
                    Ordering::Less
                });
            } else {
                return Some(if b.sign > 0 {
                    Ordering::Less
                } else {
                    Ordering::Greater
                });
            }
        }

        // zeros
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

        // sign
        if a.sign != b.sign {
            return Some(if a.sign < b.sign {
                Ordering::Less
            } else {
                Ordering::Greater
            });
        }

        // same sign: compare exponent and mantissa
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

//
// ==========================
//   FP8 E4M3
// ==========================
//

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(transparent)]
pub struct Fp8E4M3(u8);

impl Fp8E4M3 {
    /// Canonical NaN: S 1111 111_2 -> 0b0_1111_111 = 0x7F
    pub const CANONICAL_NAN: Self = Self(0x7F);

    /// Max finite: S 1111 110_2 -> 0b0_1111_110 = 0x7E
    pub const MAX_FINITE: Self = Self(0x7E);

    /// Smallest positive normal: S 0001 000_2 = 0b0_0001_000 = 0x08
    pub const MIN_NORMAL_POS: Self = Self(0x08);

    /// Smallest positive subnormal: S 0000 001_2 = 0b0_0000_001 = 0x01
    pub const MIN_SUBNORMAL_POS: Self = Self(0x01);

    #[inline]
    pub const fn from_bits(bits: u8) -> Self {
        Self(bits)
    }

    #[inline]
    pub const fn to_bits(self) -> u8 {
        self.0
    }

    /// Saturating, nearest-neighbor quantization from f32 to FP8 E4M3.
    #[inline]
    pub fn quantize(v: f32) -> Self {
        if v.is_nan() {
            return Self::CANONICAL_NAN;
        }

        // No infinities in E4M3; clamp infinities to huge finite values.
        let v = if v.is_infinite() {
            if v.is_sign_positive() {
                f32::MAX
            } else {
                -f32::MAX
            }
        } else {
            v
        };

        let mut best_bits: u8 = 0;
        let mut best_err: f32 = f32::INFINITY;

        for bits in 0u8..=255 {
            let decoded = fp8e4m3_decode(bits);

            // skip NaN encodings (E = 0xF)
            if decoded.is_nan() {
                continue;
            }

            let err = (decoded - v).abs();
            if err < best_err {
                best_err = err;
                best_bits = bits;
            } else if err == best_err {
                // tie-break: prefer even mantissa LSB
                let cur_lsb = best_bits & 0x1;
                let new_lsb = bits & 0x1;
                if new_lsb == 0 && cur_lsb == 1 {
                    best_bits = bits;
                }
            }
        }

        Self(best_bits)
    }

    /// Unpack into integer representation:
    ///
    /// Normals (1 <= e <= 14):
    ///   value = sign * (1 + m/8) * 2^(e - BIAS)
    ///         = sign * (8 + m) * 2^(e - BIAS - 3)
    ///
    /// Subnormals (e = 0, m != 0):
    ///   value = sign * (m/8) * 2^(1 - BIAS)
    ///         = sign * m * 2^(1 - BIAS - 3)
    #[inline]
    fn unpack(self) -> UnpackedE4 {
        let bits = self.to_bits();
        let s = (bits >> 7) & 0x1;
        let e = (bits >> FP8_E4M3_MAN_BITS) & 0x0F;
        let m = bits & 0x07;

        // Only S 1111 111_2 is NaN: e=0xF, m=0x7
        if e == 0x0F && m == 0x07 {
            return UnpackedE4 {
                sign: 1,
                exp: 0,
                mant: 0,
                is_zero: false,
                is_subnormal: false,
                is_nan: true,
            };
        }

        if e == 0 {
            if m == 0 {
                // Zero
                return UnpackedE4 {
                    sign: 1,
                    exp: 0,
                    mant: 0,
                    is_zero: true,
                    is_subnormal: false,
                    is_nan: false,
                };
            } else {
                // Subnormal:
                // value = sign * (m/8) * 2^(1-bias)
                //       = sign * m * 2^(1-bias-3) = sign * m * 2^-9
                let sign = if s == 0 { 1 } else { -1 };
                let exp = (1 - FP8_E4M3_BIAS - FP8_E4M3_MAN_BITS as i32) as i16; // -9
                let mant = m as u16;
                return UnpackedE4 {
                    sign,
                    exp,
                    mant,
                    is_zero: false,
                    is_subnormal: true,
                    is_nan: false,
                };
            }
        }

        // Normal (including e == 0x0F and m <= 0x06)
        //
        // value = sign * (1 + m/8) * 2^(e - bias)
        //       = sign * (8 + m) * 2^(e - bias - 3)
        let sign = if s == 0 { 1 } else { -1 };
        let exp = (e as i32 - FP8_E4M3_BIAS - FP8_E4M3_MAN_BITS as i32) as i16;
        let mant = (8 + m) as u16; // 8..15

        UnpackedE4 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
            is_nan: false,
        }
    }

    fn pack_unpacked(u: UnpackedE4) -> Self {
        if u.is_nan {
            return Fp8E4M3::CANONICAL_NAN;
        }

        if u.is_zero || u.mant == 0 {
            return Fp8E4M3::from_bits(0);
        }

        let mut sign = u.sign;
        let mut exp = u.exp;
        let mut mant = u.mant as i32;

        // normalize mant into ~[8, 15]
        while mant > 15 {
            mant = if mant & 1 != 0 {
                (mant + 1) >> 1
            } else {
                mant >> 1
            };
            exp += 1;
        }

        while mant < 8 && mant > 0 && exp > -9 {
            mant <<= 1;
            exp -= 1;
        }

        // Representable integer exponent range for normals is roughly [-9, 5]
        if exp > 5 {
            // overflow -> clamp to max finite (no infinities)
            return Fp8E4M3::MAX_FINITE;
        }
        if exp < -9 {
            // underflow -> flush to zero (tiny subnormals)
            return Fp8E4M3::from_bits(0);
        }

        if mant <= 0 {
            return Fp8E4M3::from_bits(0);
        }
        if mant > 15 {
            mant = 15;
        }

        // e = exp + bias + MAN_BITS
        let mut e = (exp + FP8_E4M3_BIAS as i16 + FP8_E4M3_MAN_BITS as i16) as u8; // 1..15
        if e == 0 {
            // shouldn't happen, but guard anyway
            return Fp8E4M3::from_bits(0);
        }
        if e > 0x0F {
            // beyond max exponent -> clamp to max finite
            return Fp8E4M3::MAX_FINITE;
        }

        let s_bit = if sign < 0 { 1u8 } else { 0u8 };
        let mut m_bits = (mant as u8).saturating_sub(8) & 0x07;

        // Avoid generating the NaN pattern S 1111 111_2:
        // if (e=0x0F, m_bits=0x07), clamp to max normal (m_bits=0x06).
        if e == 0x0F && m_bits == 0x07 {
            m_bits = 0x06;
        }

        let bits = (s_bit << 7) | ((e & 0x0F) << FP8_E4M3_MAN_BITS) | m_bits;
        Fp8E4M3::from_bits(bits)
    }
}

// Format constants
const FP8_E4M3_EXP_BITS: u32 = 4;
const FP8_E4M3_MAN_BITS: u32 = 3;
const FP8_E4M3_BIAS: i32 = 7;

#[inline]
fn fp8e4m3_decode(bits: u8) -> f32 {
    let s = (bits >> 7) & 0x1;
    let e = (bits >> FP8_E4M3_MAN_BITS) & 0x0F;
    let m = bits & 0x07;

    let sign = if s == 0 { 1.0f32 } else { -1.0f32 };

    if e == 0 {
        if m == 0 {
            return sign * 0.0;
        }
        // subnormal: exponent = 1 - bias = -6, no implicit 1
        let exp = 1 - FP8_E4M3_BIAS; // -6
        let frac = (m as f32) / (1u32 << FP8_E4M3_MAN_BITS) as f32;
        return sign * exp2i(exp) * frac; // 2^-9 .. 2^-6 * 0.875
    }

    // Only e=0xF, m=0x7 is NaN
    if e == 0x0F && m == 0x07 {
        return f32::NAN;
    }

    // Normal
    let exp = (e as i32 - FP8_E4M3_BIAS) as i32;
    let frac = 1.0 + (m as f32) / (1u32 << FP8_E4M3_MAN_BITS) as f32;
    sign * exp2i(exp) * frac
}

impl From<Fp8E4M3> for f32 {
    #[inline]
    fn from(x: Fp8E4M3) -> Self {
        fp8e4m3_decode(x.0)
    }
}

impl From<f32> for Fp8E4M3 {
    #[inline]
    fn from(v: f32) -> Self {
        Fp8E4M3::quantize(v)
    }
}

// Display for E4M3
impl fmt::Display for Fp8E4M3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let bits = self.to_bits();
        let val: f32 = (*self).into();

        if bits == Self::CANONICAL_NAN.to_bits() || val.is_nan() {
            write!(f, "Fp8E4M3(NaN, bits=0x{:02X})", bits)
        } else {
            write!(f, "Fp8E4M3({:.8e}, bits=0x{:02X})", val, bits)
        }
    }
}

// ---------- arithmetic for Fp8E4M3 via UnpackedE4 ----------

impl Neg for Fp8E4M3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let u = self.unpack();
        if u.is_nan {
            Fp8E4M3::CANONICAL_NAN
        } else if u.is_zero {
            self
        } else {
            Fp8E4M3::pack_unpacked(UnpackedE4 { sign: -u.sign, ..u })
        }
    }
}

impl Add for Fp8E4M3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let a = self.unpack();
        let b = rhs.unpack();

        if a.is_nan || b.is_nan {
            return Fp8E4M3::CANONICAL_NAN;
        }

        if a.is_zero {
            return rhs;
        }
        if b.is_zero {
            return self;
        }

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
            return Fp8E4M3::from_bits(0);
        }

        let sign = if sum < 0 { -1 } else { 1 };
        let mant = sum.abs() as u16;

        Fp8E4M3::pack_unpacked(UnpackedE4 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
            is_nan: false,
        })
    }
}

impl Sub for Fp8E4M3 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl Mul for Fp8E4M3 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let a = self.unpack();
        let b = rhs.unpack();

        if a.is_nan || b.is_nan {
            return Fp8E4M3::CANONICAL_NAN;
        }

        if a.is_zero || b.is_zero {
            return Fp8E4M3::from_bits(0);
        }

        let sign = a.sign * b.sign;
        let exp = a.exp + b.exp;
        let mant = (a.mant as i32 * b.mant as i32) as u16;

        Fp8E4M3::pack_unpacked(UnpackedE4 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
            is_nan: false,
        })
    }
}

impl Div for Fp8E4M3 {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let a = self.unpack();
        let b = rhs.unpack();

        if a.is_nan || b.is_nan {
            return Fp8E4M3::CANONICAL_NAN;
        }

        if b.is_zero {
            return Fp8E4M3::CANONICAL_NAN;
        }

        if a.is_zero {
            return Fp8E4M3::from_bits(0);
        }

        let sign = a.sign * b.sign;
        // `mant_i` is scaled by 2^3 to retain three fraction bits.
        let exp = a.exp - b.exp - FP8_E4M3_MAN_BITS as i16;

        // mant_res ~= (a.mant << 3) / b.mant
        let num = (a.mant as i32) << 3;
        let den = b.mant as i32;
        let mant_i = if num >= 0 {
            (num + den / 2) / den
        } else {
            (num - den / 2) / den
        };

        if mant_i <= 0 {
            return Fp8E4M3::from_bits(0);
        }

        let mant = mant_i as u16;

        Fp8E4M3::pack_unpacked(UnpackedE4 {
            sign,
            exp,
            mant,
            is_zero: false,
            is_subnormal: false,
            is_nan: false,
        })
    }
}

impl PartialOrd for Fp8E4M3 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let a = self.unpack();
        let b = other.unpack();

        if a.is_nan || b.is_nan {
            return None;
        }

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
            return Some(if a.sign < b.sign {
                Ordering::Less
            } else {
                Ordering::Greater
            });
        }

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
