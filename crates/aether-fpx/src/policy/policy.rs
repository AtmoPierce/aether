use core::f64::consts::FRAC_PI_4;

/* ============================================================
 * Reduction result
 * ============================================================ */

#[derive(Copy, Clone)]
pub struct ReductF64 {
    pub x: f64,     // reduced argument in [0, PI/4]
    pub octant: u8, // 0..7
}

#[derive(Copy, Clone)]
pub struct ReductF32 {
    pub x: f32,
    pub octant: u8,
}

/* ============================================================
 * Helpers
 * ============================================================ */

#[inline(always)]
fn prev_f64(x: f64) -> f64 {
    f64::from_bits(x.to_bits() - 1)
}

#[inline(always)]
fn prev_f32(x: f32) -> f32 {
    f32::from_bits(x.to_bits() - 1)
}

/* ============================================================
 * Classification
 * ============================================================ */

#[inline(always)]
fn classify_f64(x: f64) -> Option<f64> {
    if x.is_finite() {
        Some(x)
    } else {
        None
    }
}

#[inline(always)]
fn classify_f32(x: f32) -> Option<f32> {
    if x.is_finite() {
        Some(x)
    } else {
        None
    }
}

/* ============================================================
 * SIN / COS reduction to [0, PI/4]
 * ============================================================ */
// ---------- f64 ----------

// 1 / (pi/4)
const INV_FRAC_PI_4_F64: f64 = 1.0_f64 / core::f64::consts::FRAC_PI_4;

// pi/4 split (hi + lo)
const FRAC_PI_4_HI_F64: f64 =
    f64::from_bits(core::f64::consts::FRAC_PI_4.to_bits() & 0xffff_ffff_f800_0000);
const FRAC_PI_4_LO_F64: f64 = core::f64::consts::FRAC_PI_4 - FRAC_PI_4_HI_F64;

#[inline(always)]
pub fn reduce_sin_cos_f64(x: f64) -> Option<ReductF64> {
    let x = classify_f64(x)?;

    // k = round(x / (pi/4))
    let k = (x * INV_FRAC_PI_4_F64).round() as i64;
    let oct = (k & 7) as u8;

    // r = x - k*(pi/4) using hi/lo split
    let kd = k as f64;
    let r = (x - kd * FRAC_PI_4_HI_F64) - kd * FRAC_PI_4_LO_F64;

    // fold into [0, pi/4]
    let xr = if (oct & 1) == 0 {
        r
    } else {
        FRAC_PI_4_HI_F64 - r
    };

    Some(ReductF64 { x: xr, octant: oct })
}

// ---------- f32 ----------

// 1 / (pi/4)
const INV_FRAC_PI_4_F32: f32 = 1.0_f32 / core::f32::consts::FRAC_PI_4;

const FRAC_PI_4_F32: f32 = core::f32::consts::FRAC_PI_4;

// pi/4 split (hi + lo)
const FRAC_PI_4_HI_F32: f32 = f32::from_bits(core::f32::consts::FRAC_PI_4.to_bits() & 0xffff_0000);
const FRAC_PI_4_LO_F32: f32 = core::f32::consts::FRAC_PI_4 - FRAC_PI_4_HI_F32;

#[inline(always)]
pub fn reduce_sin_cos_f32(x: f32) -> Option<ReductF32> {
    let x = classify_f32(x)?;

    // k = round(x / (pi/4))
    let k = (x * INV_FRAC_PI_4_F32).round() as i32;
    let oct = (k & 7) as u8;

    // r = x - k*(pi/4) using hi/lo split
    let kd = k as f32;
    let r = (x - kd * FRAC_PI_4_HI_F32) - kd * FRAC_PI_4_LO_F32;

    // fold into [0, pi/4]
    let xr = if (oct & 1) == 0 {
        r
    } else {
        (FRAC_PI_4_HI_F32 - r) + FRAC_PI_4_LO_F32
    };

    Some(ReductF32 { x: xr, octant: oct })
}

/* ============================================================
 * TAN policy
 * ============================================================ */

#[inline(always)]
pub fn apply_policy_tan_f64(x: f64) -> f64 {
    let x = match classify_f64(x) {
        Some(v) => v,
        None => return f64::NAN,
    };

    if x == FRAC_PI_4 {
        prev_f64(x)
    } else if x == -FRAC_PI_4 {
        -prev_f64(-x)
    } else {
        x
    }
}

#[inline(always)]
pub fn apply_policy_tan_f32(x: f32) -> f32 {
    let x = match classify_f32(x) {
        Some(v) => v,
        None => return f32::NAN,
    };

    if x == FRAC_PI_4_F32 {
        prev_f32(x)
    } else if x == -FRAC_PI_4_F32 {
        -prev_f32(-x)
    } else {
        x
    }
}

/* ============================================================
 * SIN / COS policy
 * ============================================================ */

#[inline(always)]
pub fn apply_policy_sin_f64(x: f64) -> f64 {
    match classify_f64(x) {
        Some(v) => v,
        None => f64::NAN,
    }
}

#[inline(always)]
pub fn apply_policy_cos_f64(x: f64) -> f64 {
    match classify_f64(x) {
        Some(v) => v,
        None => f64::NAN,
    }
}

#[inline(always)]
pub fn apply_policy_sin_f32(x: f32) -> f32 {
    match classify_f32(x) {
        Some(v) => v,
        None => f32::NAN,
    }
}

#[inline(always)]
pub fn apply_policy_cos_f32(x: f32) -> f32 {
    match classify_f32(x) {
        Some(v) => v,
        None => f32::NAN,
    }
}

/* ============================================================
 * Reconstruction: SIN / COS (f64)
 * ============================================================ */

#[inline(always)]
pub fn reconstruct_sin_f64(oct: u8, sin_x: f64, cos_x: f64) -> f64 {
    match oct & 7 {
        0 => sin_x,
        1 => cos_x,
        2 => cos_x,
        3 => sin_x,
        4 => -sin_x,
        5 => -cos_x,
        6 => -cos_x,
        7 => -sin_x,
        _ => unreachable!(),
    }
}

#[inline(always)]
pub fn reconstruct_cos_f64(oct: u8, sin_x: f64, cos_x: f64) -> f64 {
    match oct & 7 {
        0 => cos_x,
        1 => sin_x,
        2 => -sin_x,
        3 => -cos_x,
        4 => -cos_x,
        5 => -sin_x,
        6 => sin_x,
        7 => cos_x,
        _ => unreachable!(),
    }
}

/* ============================================================
 * Reconstruction: SIN / COS (f32)
 * ============================================================ */

#[inline(always)]
pub fn reconstruct_sin_f32(oct: u8, sin_x: f32, cos_x: f32) -> f32 {
    match oct & 7 {
        0 => sin_x,
        1 => cos_x,
        2 => cos_x,
        3 => sin_x,
        4 => -sin_x,
        5 => -cos_x,
        6 => -cos_x,
        7 => -sin_x,
        _ => unreachable!(),
    }
}

#[inline(always)]
pub fn reconstruct_cos_f32(oct: u8, sin_x: f32, cos_x: f32) -> f32 {
    match oct & 7 {
        0 => cos_x,
        1 => sin_x,
        2 => -sin_x,
        3 => -cos_x,
        4 => -cos_x,
        5 => -sin_x,
        6 => sin_x,
        7 => cos_x,
        _ => unreachable!(),
    }
}
