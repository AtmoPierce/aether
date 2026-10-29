/// Wrap x into [a, b) for f32
/// Assumes b > a and p = b-a is finite and nonzero.
#[inline(always)]
pub fn wrap_to_interval_f32(x: f32, a: f32, b: f32) -> f32 {
    let p = b - a;
    let y = x - a;

    // k = floor(y/p) computed via trunc cast after biasing for sign.
    // This avoids calling libm floor(). For negative y, trunc rounds toward 0,
    // so subtract 1 when there is a fractional part.
    let q = y / p;
    let mut k = q as i32; // trunc toward 0
    if q < 0.0 && (k as f32) != q {
        k -= 1;
    }

    // r = y - k*p in [0, p)
    let r = y - (k as f32) * p;
    r + a
}

/// Wrap x into [a, b) for f64
/// Assumes b > a and p = b-a is finite and nonzero.
#[inline(always)]
pub fn wrap_to_interval_f64(x: f64, a: f64, b: f64) -> f64 {
    let p = b - a;
    let y = x - a;

    let q = y / p;
    let mut k = q as i64; // trunc toward 0
    if q < 0.0 && (k as f64) != q {
        k -= 1;
    }

    let r = y - (k as f64) * p;
    r + a
}

/// Fast wrap x into [0, TAU) for f32.
/// No division; uses multiply by INV_TAU + round + FMA/sub + 1–2 fixups.
#[inline(always)]
pub fn wrap_0_2pi_f32_fast(x: f32) -> f32 {
    const TAU: f32 = core::f32::consts::TAU;
    const INV_TAU: f32 = 0.15915494309189535_f32; // 1/(2PI)

    // k ~= round(x / TAU)
    let k = (x * INV_TAU).round();
    // r = x - k*TAU
    let mut r = x - k * TAU;

    // Map into [0, TAU)
    if r < 0.0 {
        r += TAU;
    }
    if r >= TAU {
        r -= TAU;
    }
    r
}

/// Fast wrap x into [0, TAU) for f64 (unsplit version).
#[inline(always)]
pub fn wrap_0_2pi_f64_fast(x: f64) -> f64 {
    const TAU: f64 = core::f64::consts::TAU;
    const INV_TAU: f64 = 0.15915494309189535_f64;

    let k = (x * INV_TAU).round();
    let mut r = x - k * TAU;

    if r < 0.0 {
        r += TAU;
    }
    if r >= TAU {
        r -= TAU;
    }
    r
}

/// Fast wrap x into [0, PI) for tan reducers (f32).
#[inline(always)]
pub fn wrap_0_pi_f32_fast(x: f32) -> f32 {
    const PI: f32 = core::f32::consts::PI;
    const INV_PI: f32 = 0.3183098861837907_f32; // 1/PI

    let k = (x * INV_PI).round();
    let mut r = x - k * PI;

    if r < 0.0 {
        r += PI;
    }
    if r >= PI {
        r -= PI;
    }
    r
}

/// Fast wrap x into [0, PI) for tan reducers (f64).
#[inline(always)]
pub fn wrap_0_pi_f64_fast(x: f64) -> f64 {
    const PI: f64 = core::f64::consts::PI;
    const INV_PI: f64 = 0.3183098861837907_f64;

    let k = (x * INV_PI).round();
    let mut r = x - k * PI;

    if r < 0.0 {
        r += PI;
    }
    if r >= PI {
        r -= PI;
    }
    r
}

#[inline(always)]
pub fn wrap_0_2pi_f64(x: f64) -> f64 {
    wrap_0_2pi_f64_fast(x)
}

#[inline(always)]
pub fn wrap_0_2pi_f32(x: f32) -> f32 {
    wrap_0_2pi_f32_fast(x)
}
