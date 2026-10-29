use core::{f128, f16};

#[inline(always)]
pub fn rem_pio2_f64(x: f64) -> (f64, u8) {
    crate::trig::payne_hanek::rem_pio2_f64(x)
}

#[inline(always)]
pub fn rem_pio2_f32(x: f32) -> (f32, u8) {
    crate::trig::payne_hanek::rem_pio2_f32(x)
}

#[inline(always)]
pub fn rem_pio2_f128(x: f128) -> (f128, u8) {
    crate::trig::payne_hanek::rem_pio2_f128(x)
}

#[inline(always)]
pub fn rem_pio2_f16(x: f16) -> (f16, u8) {
    crate::trig::payne_hanek::rem_pio2_f16(x)
}

#[inline(always)]
pub fn rem_pio2_tan_f64(x: f64) -> (f64, u8) {
    crate::trig::payne_hanek::rem_pio2_tan_f64(x)
}

#[inline(always)]
pub fn rem_pio2_tan_f32(x: f32) -> (f32, u8) {
    crate::trig::payne_hanek::rem_pio2_tan_f32(x)
}

#[inline(always)]
pub fn rem_pio2_tan_f128(x: f128) -> (f128, u8) {
    crate::trig::payne_hanek::rem_pio2_tan_f128(x)
}

#[inline(always)]
pub fn rem_pio2_tan_f16(x: f16) -> (f16, u8) {
    crate::trig::payne_hanek::rem_pio2_tan_f16(x)
}

#[inline(always)]
pub fn eval_horner_f64_slice<const N: usize>(x: f64, cs: &[f64; N]) -> f64 {
    let mut p = 0.0_f64;
    for &c in cs.iter().rev() {
        p = p.mul_add(x, c);
    }
    p
}

#[inline(always)]
pub fn eval_estrin_f64_slice<const N: usize>(x: f64, cs: &[f64; N]) -> f64 {
    if N == 0 {
        return 0.0_f64;
    }

    let mut scratch = [0.0_f64; N];
    scratch.copy_from_slice(cs);

    let mut len = N;
    let mut x_pow = x;
    while len > 1 {
        let mut write = 0;
        let mut read = 0;

        while read + 1 < len {
            scratch[write] = scratch[read] + scratch[read + 1] * x_pow;
            write += 1;
            read += 2;
        }

        if read < len {
            scratch[write] = scratch[read];
            write += 1;
        }

        len = write;
        x_pow *= x_pow;
    }

    scratch[0]
}

#[inline(always)]
pub fn eval_horner_f32_slice<const N: usize>(x: f32, cs: &[f32; N]) -> f32 {
    let mut p = 0.0_f32;
    for &c in cs.iter().rev() {
        p = p.mul_add(x, c);
    }
    p
}

#[inline(always)]
pub fn eval_estrin_f32_slice<const N: usize>(x: f32, cs: &[f32; N]) -> f32 {
    if N == 0 {
        return 0.0_f32;
    }

    let mut scratch = [0.0_f32; N];
    scratch.copy_from_slice(cs);

    let mut len = N;
    let mut x_pow = x;
    while len > 1 {
        let mut write = 0;
        let mut read = 0;

        while read + 1 < len {
            scratch[write] = scratch[read] + scratch[read + 1] * x_pow;
            write += 1;
            read += 2;
        }

        if read < len {
            scratch[write] = scratch[read];
            write += 1;
        }

        len = write;
        x_pow *= x_pow;
    }

    scratch[0]
}

#[inline(always)]
pub fn eval_horner_f16_slice<const N: usize>(x: f16, cs: &[f16; N]) -> f16 {
    let mut p = 0.0_f16;
    for &c in cs.iter().rev() {
        p = p * x + c;
    }
    p
}

#[inline(always)]
pub fn eval_horner_f128_slice<const N: usize>(x: f128, cs: &[f128; N]) -> f128 {
    let mut p = 0.0_f128;
    for &c in cs.iter().rev() {
        p = p * x + c;
    }
    p
}
