// aether-fpx-gen/src/core/ulp.rs

use core::num::FpCategory;
use core::{f128, f16};

// =============================================================================
// f64
// =============================================================================

#[inline]
pub fn ordered_u64_from_f64(x: f64) -> u64 {
    let u = x.to_bits();

    if (u >> 63) != 0 {
        !u
    } else {
        u | (1u64 << 63)
    }
}

#[inline]
pub fn f64_from_ordered_u64(o: u64) -> f64 {
    let u = if (o >> 63) != 0 {
        o & !(1u64 << 63)
    } else {
        !o
    };

    f64::from_bits(u)
}

#[inline]
fn canonical_zero_f64(x: f64) -> f64 {
    if x == 0.0 {
        0.0
    } else {
        x
    }
}

/// ULP distance for `f64`.
///
/// Behavior:
/// - `NaN` returns `u64::MAX`.
/// - `+0.0` and `-0.0` are canonicalized to `+0.0`.
/// - Same-sign values use ordered-bit distance.
/// - Opposite-sign finite values are measured through zero:
///
/// ```text
/// ulp_dist(-a, b) = ulp_dist(-a, 0.0) + ulp_dist(0.0, b)
/// ```
///
/// This avoids the huge artificial spike caused by comparing ordered integer
/// encodings directly across the sign boundary.
#[inline]
pub fn ulp_dist_f64(a: f64, b: f64) -> u64 {
    if a.is_nan() || b.is_nan() {
        return u64::MAX;
    }

    let a = canonical_zero_f64(a);
    let b = canonical_zero_f64(b);

    if a.is_sign_negative() != b.is_sign_negative() {
        let z = ordered_u64_from_f64(0.0);

        let da = ordered_u64_from_f64(a.abs()).abs_diff(z);
        let db = ordered_u64_from_f64(b.abs()).abs_diff(z);

        da.saturating_add(db)
    } else {
        ordered_u64_from_f64(a).abs_diff(ordered_u64_from_f64(b))
    }
}

#[inline]
pub fn valid_ulp_dist_f64(a: f64, b: f64) -> Option<u64> {
    if a.is_nan() || b.is_nan() {
        None
    } else {
        Some(ulp_dist_f64(a, b))
    }
}

// =============================================================================
// f32
// =============================================================================

#[inline]
pub fn ordered_u32_from_f32(x: f32) -> u32 {
    let u = x.to_bits();

    if (u >> 31) != 0 {
        !u
    } else {
        u | (1u32 << 31)
    }
}

#[inline]
pub fn f32_from_ordered_u32(o: u32) -> f32 {
    let u = if (o >> 31) != 0 {
        o & !(1u32 << 31)
    } else {
        !o
    };

    f32::from_bits(u)
}

#[inline]
fn canonical_zero_f32(x: f32) -> f32 {
    if x == 0.0 {
        0.0
    } else {
        x
    }
}

/// ULP distance for `f32`.
///
/// See `ulp_dist_f64` for sign-crossing behavior.
#[inline]
pub fn ulp_dist_f32(a: f32, b: f32) -> u32 {
    if a.is_nan() || b.is_nan() {
        return u32::MAX;
    }

    let a = canonical_zero_f32(a);
    let b = canonical_zero_f32(b);

    if a.is_sign_negative() != b.is_sign_negative() {
        let z = ordered_u32_from_f32(0.0);

        let da = ordered_u32_from_f32(a.abs()).abs_diff(z);
        let db = ordered_u32_from_f32(b.abs()).abs_diff(z);

        da.saturating_add(db)
    } else {
        ordered_u32_from_f32(a).abs_diff(ordered_u32_from_f32(b))
    }
}

#[inline]
pub fn valid_ulp_dist_f32(a: f32, b: f32) -> Option<u32> {
    if a.is_nan() || b.is_nan() {
        None
    } else {
        Some(ulp_dist_f32(a, b))
    }
}

// =============================================================================
// f16
// =============================================================================

#[inline]
pub fn ordered_u16_from_f16(x: f16) -> u16 {
    let u = x.to_bits();

    if (u >> 15) != 0 {
        !u
    } else {
        u | (1u16 << 15)
    }
}

#[inline]
pub fn f16_from_ordered_u16(o: u16) -> f16 {
    let u = if (o >> 15) != 0 {
        o & !(1u16 << 15)
    } else {
        !o
    };

    f16::from_bits(u)
}

#[inline]
fn canonical_zero_f16(x: f16) -> f16 {
    if x == f16::from_bits(0) {
        f16::from_bits(0)
    } else {
        x
    }
}

/// ULP distance for `f16`.
///
/// See `ulp_dist_f64` for sign-crossing behavior.
#[inline]
pub fn ulp_dist_f16(a: f16, b: f16) -> u16 {
    if a.is_nan() || b.is_nan() {
        return u16::MAX;
    }

    let a = canonical_zero_f16(a);
    let b = canonical_zero_f16(b);

    if a.is_sign_negative() != b.is_sign_negative() {
        let z = ordered_u16_from_f16(f16::from_bits(0));

        let da = ordered_u16_from_f16(a.abs()).abs_diff(z);
        let db = ordered_u16_from_f16(b.abs()).abs_diff(z);

        da.saturating_add(db)
    } else {
        ordered_u16_from_f16(a).abs_diff(ordered_u16_from_f16(b))
    }
}

#[inline]
pub fn valid_ulp_dist_f16(a: f16, b: f16) -> Option<u16> {
    if a.is_nan() || b.is_nan() {
        None
    } else {
        Some(ulp_dist_f16(a, b))
    }
}

// =============================================================================
// f128
// =============================================================================

#[inline]
pub fn ordered_u128_from_f128(x: f128) -> u128 {
    let u = x.to_bits();

    if (u >> 127) != 0 {
        !u
    } else {
        u | (1u128 << 127)
    }
}

#[inline]
pub fn f128_from_ordered_u128(o: u128) -> f128 {
    let u = if (o >> 127) != 0 {
        o & !(1u128 << 127)
    } else {
        !o
    };

    f128::from_bits(u)
}

#[inline]
fn canonical_zero_f128(x: f128) -> f128 {
    if x == 0.0 {
        0.0
    } else {
        x
    }
}

/// ULP distance for `f128`.
///
/// See `ulp_dist_f64` for sign-crossing behavior.
#[inline]
pub fn ulp_dist_f128(a: f128, b: f128) -> u128 {
    if a.is_nan() || b.is_nan() {
        return u128::MAX;
    }

    let a = canonical_zero_f128(a);
    let b = canonical_zero_f128(b);

    if a.is_sign_negative() != b.is_sign_negative() {
        let z = ordered_u128_from_f128(0.0);

        let da = ordered_u128_from_f128(a.abs()).abs_diff(z);
        let db = ordered_u128_from_f128(b.abs()).abs_diff(z);

        da.saturating_add(db)
    } else {
        ordered_u128_from_f128(a).abs_diff(ordered_u128_from_f128(b))
    }
}

#[inline]
pub fn valid_ulp_dist_f128(a: f128, b: f128) -> Option<u128> {
    if a.is_nan() || b.is_nan() {
        None
    } else {
        Some(ulp_dist_f128(a, b))
    }
}

// =============================================================================
// Debug structs
// =============================================================================

#[derive(Clone, Copy, Debug)]
pub struct UlpDebugF64 {
    pub a: f64,
    pub b: f64,
    pub a_bits: u64,
    pub b_bits: u64,
    pub a_class: FpCategory,
    pub b_class: FpCategory,
    pub ulp: u64,
}

#[inline]
pub fn ulp_debug_f64(a: f64, b: f64) -> UlpDebugF64 {
    let a_bits = a.to_bits();
    let b_bits = b.to_bits();
    let a_class = a.classify();
    let b_class = b.classify();
    let ulp = ulp_dist_f64(a, b);

    UlpDebugF64 {
        a,
        b,
        a_bits,
        b_bits,
        a_class,
        b_class,
        ulp,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct UlpDebugF32 {
    pub a: f32,
    pub b: f32,
    pub a_bits: u32,
    pub b_bits: u32,
    pub a_class: FpCategory,
    pub b_class: FpCategory,
    pub ulp: u32,
}

#[inline]
pub fn ulp_debug_f32(a: f32, b: f32) -> UlpDebugF32 {
    let a_bits = a.to_bits();
    let b_bits = b.to_bits();
    let a_class = a.classify();
    let b_class = b.classify();
    let ulp = ulp_dist_f32(a, b);

    UlpDebugF32 {
        a,
        b,
        a_bits,
        b_bits,
        a_class,
        b_class,
        ulp,
    }
}

// =============================================================================
// Stepping
// =============================================================================

#[inline]
pub fn next_f64(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }

    let o = if x == 0.0 {
        ordered_u64_from_f64(0.0)
    } else {
        ordered_u64_from_f64(x)
    };

    f64_from_ordered_u64(o + 1)
}

#[inline]
pub fn prev_f64(x: f64) -> f64 {
    if x.is_nan() || x == f64::NEG_INFINITY {
        return x;
    }
    if x == 0.0 {
        return f64::from_bits((1_u64 << 63) | 1);
    }

    let o = ordered_u64_from_f64(x);

    f64_from_ordered_u64(o - 1)
}

#[inline]
pub fn next_f32(x: f32) -> f32 {
    if x.is_nan() || x == f32::INFINITY {
        return x;
    }

    let o = if x == 0.0 {
        ordered_u32_from_f32(0.0)
    } else {
        ordered_u32_from_f32(x)
    };

    f32_from_ordered_u32(o + 1)
}

#[inline]
pub fn prev_f32(x: f32) -> f32 {
    if x.is_nan() || x == f32::NEG_INFINITY {
        return x;
    }
    if x == 0.0 {
        return f32::from_bits((1_u32 << 31) | 1);
    }

    let o = ordered_u32_from_f32(x);

    f32_from_ordered_u32(o - 1)
}

#[inline]
pub fn next_f16(x: f16) -> f16 {
    if x.is_nan() || x == f16::INFINITY {
        return x;
    }

    let o = if x == f16::from_bits(0) {
        ordered_u16_from_f16(f16::from_bits(0))
    } else {
        ordered_u16_from_f16(x)
    };

    f16_from_ordered_u16(o + 1)
}

#[inline]
pub fn prev_f16(x: f16) -> f16 {
    if x.is_nan() || x == f16::NEG_INFINITY {
        return x;
    }
    if x == f16::from_bits(0) {
        return f16::from_bits((1_u16 << 15) | 1);
    }

    let o = ordered_u16_from_f16(x);

    f16_from_ordered_u16(o - 1)
}

#[inline]
pub fn next_f128(x: f128) -> f128 {
    if x.is_nan() || x == f128::INFINITY {
        return x;
    }

    let o = if x == 0.0 {
        ordered_u128_from_f128(0.0)
    } else {
        ordered_u128_from_f128(x)
    };

    f128_from_ordered_u128(o + 1)
}

#[inline]
pub fn prev_f128(x: f128) -> f128 {
    if x.is_nan() || x == f128::NEG_INFINITY {
        return x;
    }
    if x == 0.0 {
        return f128::from_bits((1_u128 << 127) | 1);
    }

    let o = ordered_u128_from_f128(x);

    f128_from_ordered_u128(o - 1)
}

// =============================================================================
// Spans
// =============================================================================

/// Raw ordered-representation span.
///
/// This is not sign-boundary-corrected. For error metrics, prefer
/// `ulp_dist_f64`.
#[inline]
pub fn ulp_span_f64(a: f64, b: f64) -> u64 {
    ordered_u64_from_f64(a).abs_diff(ordered_u64_from_f64(b))
}

/// Raw ordered-representation span.
///
/// This is not sign-boundary-corrected. For error metrics, prefer
/// `ulp_dist_f32`.
#[inline]
pub fn ulp_span_f32(a: f32, b: f32) -> u32 {
    ordered_u32_from_f32(a).abs_diff(ordered_u32_from_f32(b))
}

/// Raw ordered-representation span.
///
/// This is not sign-boundary-corrected. For error metrics, prefer
/// `ulp_dist_f16`.
#[inline]
pub fn ulp_span_f16(a: f16, b: f16) -> u16 {
    ordered_u16_from_f16(a).abs_diff(ordered_u16_from_f16(b))
}

/// Raw ordered-representation span.
///
/// This is not sign-boundary-corrected. For error metrics, prefer
/// `ulp_dist_f128`.
#[inline]
pub fn ulp_span_f128(a: f128, b: f128) -> u128 {
    ordered_u128_from_f128(a).abs_diff(ordered_u128_from_f128(b))
}

// =============================================================================
// Generators
// =============================================================================

pub struct OrderedF64Iter {
    cur: u64,
    end: u64,
    done: bool,
}

impl OrderedF64Iter {
    #[inline]
    pub fn new(lo: f64, hi: f64) -> Self {
        let a = ordered_u64_from_f64(lo);
        let b = ordered_u64_from_f64(hi);

        Self {
            cur: a,
            end: b,
            done: a > b,
        }
    }
}

impl Iterator for OrderedF64Iter {
    type Item = f64;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let x = f64_from_ordered_u64(self.cur);

        if self.cur == self.end {
            self.done = true;
        } else {
            self.cur += 1;
        }

        Some(x)
    }
}

pub struct OrderedF32Iter {
    cur: u32,
    end: u32,
    done: bool,
}

impl OrderedF32Iter {
    #[inline]
    pub fn new(lo: f32, hi: f32) -> Self {
        let a = ordered_u32_from_f32(lo);
        let b = ordered_u32_from_f32(hi);

        Self {
            cur: a,
            end: b,
            done: a > b,
        }
    }
}

impl Iterator for OrderedF32Iter {
    type Item = f32;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let x = f32_from_ordered_u32(self.cur);

        if self.cur == self.end {
            self.done = true;
        } else {
            self.cur += 1;
        }

        Some(x)
    }
}

pub struct OrderedF16Iter {
    cur: u16,
    end: u16,
    done: bool,
}

impl OrderedF16Iter {
    #[inline]
    pub fn new(lo: f16, hi: f16) -> Self {
        let a = ordered_u16_from_f16(lo);
        let b = ordered_u16_from_f16(hi);

        Self {
            cur: a,
            end: b,
            done: a > b,
        }
    }
}

impl Iterator for OrderedF16Iter {
    type Item = f16;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let x = f16_from_ordered_u16(self.cur);

        if self.cur == self.end {
            self.done = true;
        } else {
            self.cur += 1;
        }

        Some(x)
    }
}

// Backward-compatible aliases from your original file.

#[inline]
pub fn ordered_f64_bits(x: f64) -> u64 {
    ordered_u64_from_f64(x)
}

#[inline]
pub fn f64_from_ordered_bits(u: u64) -> f64 {
    f64_from_ordered_u64(u)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f64_order_roundtrip() {
        let values = [
            f64::NEG_INFINITY,
            -1.0,
            -f64::MIN_POSITIVE,
            -0.0,
            0.0,
            f64::MIN_POSITIVE,
            1.0,
            f64::INFINITY,
        ];

        for x in values {
            let o = ordered_u64_from_f64(x);
            let y = f64_from_ordered_u64(o);

            if x == 0.0 {
                assert_eq!(y, 0.0);
            } else {
                assert_eq!(x.to_bits(), y.to_bits());
            }
        }
    }

    #[test]
    fn f32_order_roundtrip() {
        let values = [
            f32::NEG_INFINITY,
            -1.0,
            -f32::MIN_POSITIVE,
            -0.0,
            0.0,
            f32::MIN_POSITIVE,
            1.0,
            f32::INFINITY,
        ];

        for x in values {
            let o = ordered_u32_from_f32(x);
            let y = f32_from_ordered_u32(o);

            if x == 0.0 {
                assert_eq!(y, 0.0);
            } else {
                assert_eq!(x.to_bits(), y.to_bits());
            }
        }
    }

    #[test]
    fn valid_ulp_is_none_for_nan() {
        assert_eq!(valid_ulp_dist_f32(f32::NAN, 1.0), None);
        assert_eq!(valid_ulp_dist_f32(1.0, f32::NAN), None);

        assert_eq!(valid_ulp_dist_f64(f64::NAN, 1.0), None);
        assert_eq!(valid_ulp_dist_f64(1.0, f64::NAN), None);
    }

    #[test]
    fn valid_ulp_is_defined_for_nonzero_sign_crossings() {
        assert!(valid_ulp_dist_f32(1.0e-7_f32, -1.0e-7_f32).is_some());
        assert!(valid_ulp_dist_f64(1.0e-15_f64, -1.0e-15_f64).is_some());
    }

    #[test]
    fn valid_ulp_is_defined_when_one_side_is_zero() {
        assert_eq!(valid_ulp_dist_f32(0.0_f32, -0.0_f32), Some(0));
        assert!(valid_ulp_dist_f32(0.0_f32, 1.0e-7_f32).is_some());

        assert_eq!(valid_ulp_dist_f64(0.0_f64, -0.0_f64), Some(0));
        assert!(valid_ulp_dist_f64(0.0_f64, 1.0e-15_f64).is_some());
    }

    #[test]
    fn sign_crossing_distance_is_through_zero_f32() {
        let x = 1.0e-7_f32;
        let y = 2.0e-7_f32;

        let direct = ulp_dist_f32(-x, y);
        let through_zero = ulp_dist_f32(-x, 0.0) + ulp_dist_f32(0.0, y);

        assert_eq!(direct, through_zero);
    }

    #[test]
    fn sign_crossing_distance_is_through_zero_f64() {
        let x = 1.0e-15_f64;
        let y = 2.0e-15_f64;

        let direct = ulp_dist_f64(-x, y);
        let through_zero = ulp_dist_f64(-x, 0.0) + ulp_dist_f64(0.0, y);

        assert_eq!(direct, through_zero);
    }

    #[test]
    fn next_prev_f64_around_zero() {
        let p = next_f64(0.0);
        let n = prev_f64(0.0);

        assert!(p > 0.0);
        assert!(n < 0.0);

        assert_eq!(ulp_dist_f64(0.0, p), 1);
        assert_eq!(ulp_dist_f64(0.0, n), 1);
    }

    #[test]
    fn next_prev_f32_around_zero() {
        let p = next_f32(0.0);
        let n = prev_f32(0.0);

        assert!(p > 0.0);
        assert!(n < 0.0);

        assert_eq!(ulp_dist_f32(0.0, p), 1);
        assert_eq!(ulp_dist_f32(0.0, n), 1);
    }

    #[test]
    fn next_prev_f16_and_f128_around_zero() {
        let p16 = next_f16(f16::from_bits(0));
        let n16 = prev_f16(f16::from_bits(0));
        assert!(p16 > f16::from_bits(0));
        assert!(n16 < f16::from_bits(0));
        assert_eq!(ulp_dist_f16(f16::from_bits(0), p16), 1);
        assert_eq!(ulp_dist_f16(f16::from_bits(0), n16), 1);

        let p128 = next_f128(0.0);
        let n128 = prev_f128(0.0);
        assert!(p128 > 0.0);
        assert!(n128 < 0.0);
        assert_eq!(ulp_dist_f128(0.0, p128), 1);
        assert_eq!(ulp_dist_f128(0.0, n128), 1);
    }

    #[test]
    fn debug_matches_metric() {
        let d32 = ulp_debug_f32(-1.0e-7, 1.0e-7);
        assert_eq!(d32.ulp, ulp_dist_f32(-1.0e-7, 1.0e-7));

        let d64 = ulp_debug_f64(-1.0e-15, 1.0e-15);
        assert_eq!(d64.ulp, ulp_dist_f64(-1.0e-15, 1.0e-15));
    }
}
