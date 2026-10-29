use core::f128;

// Translated from glibc
// x = (-1)^s * 2^(exp - 16383) * (1 + mant/2^112)
// x = m * 2^e, with m in [0.5, 1)
fn frexp_f128(x: f128) -> (f128, i32) {
    let bits = x.to_bits();

    const EXP_MASK: u128 = 0x7FFF_0000_0000_0000_0000_0000_0000_0000;
    const MANT_MASK: u128 = 0x0000_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF;
    const SIGN_MASK: u128 = 0x8000_0000_0000_0000_0000_0000_0000_0000;

    let exp = ((bits & EXP_MASK) >> 112) as i32;
    let mant = bits & MANT_MASK;
    let sign = bits & SIGN_MASK;

    // zero / subnormal: for our sqrt range (>0) this is basically irrelevant,
    // but let's be reasonable.
    if exp == 0 {
        return (x, 0);
    }

    // unbiased exponent
    let e_unbiased = exp - 16383;

    // m1 in [1,2): set exponent to 16383, keep mantissa
    let m1_bits = sign | (16383u128 << 112) | mant;
    let m1 = f128::from_bits(m1_bits);

    // Shift by 1 bit: m = m1 / 2 in [0.5,1), e_out = e_unbiased + 1
    let m = m1 * (0.5_f128);
    let e_out = e_unbiased + 1;

    (m, e_out)
}

pub fn dsqrt_f128(x: f128) -> f128 {
    let zero = f128::from_bits(0);

    if x == zero {
        return zero;
    }
    if x < zero {
        return f128::NAN;
    }

    // x = m * 2^e, with m in [0.5,1)
    let (m, mut e) = frexp_f128(x);

    // Initial guess from f64 sqrt of mantissa
    let mut y = (m as f64).sqrt() as f128;
    let half: f128 = 0.5f64 as f128;

    // Newton iterations on m:
    // y_{n+1} = 0.5 * (y + m / y)
    for _ in 0..10 {
        y = half * (y + m / y);
    }

    let mut result = y;

    // If e is odd, pull out a sqrt(2) factor
    if (e & 1) != 0 {
        const SQRT2_BITS: u128 = 0x3FFF_6A09_E667_F3BC_C908_B2FB_1366_EA95;
        let sqrt2 = f128::from_bits(SQRT2_BITS);
        result *= sqrt2;
        e -= 1;
    }

    // Now e is even: sqrt(x) = result * 2^(e/2)
    let two: f128 = 2.0f64 as f128;
    result * two.powi(e / 2)
}

pub fn sqrt_f(x: f64) -> f64 {
    {
        dsqrt_f128(x as f128) as f64
    }
}

pub fn sqrt_df(x: f64) -> f64 {
    {
        let y = dsqrt_f128(x as f128);
        let one: f128 = 1.0_f128;
        let two: f128 = 2.0_f128;
        (one / (two * y)) as f64
    }
}
