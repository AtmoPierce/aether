use core::arch::aarch64::*;

const F32_LANES: usize = 4;
const F64_LANES: usize = 2;

#[inline]
#[target_feature(enable = "neon")]
pub unsafe fn eval_blocked_horner_f32<const N: usize>(x: f32, coefficients: &[f32; N]) -> f32 {
    if N == 0 {
        return 0.0;
    }
    let mut t = 1.0_f32;
    for _ in 0..F32_LANES {
        t *= x;
    }
    let t_vector = vdupq_n_f32(t);
    let mut partials = vdupq_n_f32(0.0);
    let block_count = N.div_ceil(F32_LANES);
    for block in (0..block_count).rev() {
        let start = block * F32_LANES;
        let packed_coefficients = if start + F32_LANES <= N {
            unsafe { vld1q_f32(coefficients.as_ptr().add(start)) }
        } else {
            let mut tail = [0.0_f32; F32_LANES];
            tail[..N - start].copy_from_slice(&coefficients[start..]);
            unsafe { vld1q_f32(tail.as_ptr()) }
        };
        partials = vfmaq_f32(packed_coefficients, partials, t_vector);
    }

    let mut powers = [1.0_f32; F32_LANES];
    for lane in 1..F32_LANES {
        powers[lane] = powers[lane - 1] * x;
    }
    let weights = unsafe { vld1q_f32(powers.as_ptr()) };
    let weighted = vmulq_f32(partials, weights);
    let mut lanes = [0.0_f32; F32_LANES];
    unsafe { vst1q_f32(lanes.as_mut_ptr(), weighted) };
    lanes.into_iter().sum()
}

#[inline]
#[target_feature(enable = "neon")]
pub unsafe fn eval_horizontal_estrin_f32<const N: usize>(x: f32, coefficients: &[f32; N]) -> f32 {
    if N == 0 {
        return 0.0;
    }
    let mut scratch = *coefficients;
    let mut len = N;
    let mut x_power = x;
    while len > 1 {
        let next_len = len.div_ceil(2);
        let mut output = 0;
        while output < next_len {
            let lane_count = (next_len - output).min(F32_LANES);
            let mut even = [0.0_f32; F32_LANES];
            let mut odd = [0.0_f32; F32_LANES];
            for lane in 0..lane_count {
                let input = 2 * (output + lane);
                even[lane] = scratch[input];
                if input + 1 < len {
                    odd[lane] = scratch[input + 1];
                }
            }
            let even = unsafe { vld1q_f32(even.as_ptr()) };
            let odd = unsafe { vld1q_f32(odd.as_ptr()) };
            let combined = vfmaq_f32(even, odd, vdupq_n_f32(x_power));
            let mut result = [0.0_f32; F32_LANES];
            unsafe { vst1q_f32(result.as_mut_ptr(), combined) };
            scratch[output..output + lane_count].copy_from_slice(&result[..lane_count]);
            output += lane_count;
        }
        len = next_len;
        x_power *= x_power;
    }
    scratch[0]
}

#[inline]
#[target_feature(enable = "neon")]
pub unsafe fn eval_blocked_horner_f64<const N: usize>(x: f64, coefficients: &[f64; N]) -> f64 {
    if N == 0 {
        return 0.0;
    }
    let t = x * x;
    let t_vector = vdupq_n_f64(t);
    let mut partials = vdupq_n_f64(0.0);
    let block_count = N.div_ceil(F64_LANES);
    for block in (0..block_count).rev() {
        let start = block * F64_LANES;
        let packed_coefficients = if start + F64_LANES <= N {
            unsafe { vld1q_f64(coefficients.as_ptr().add(start)) }
        } else {
            let mut tail = [0.0_f64; F64_LANES];
            tail[..N - start].copy_from_slice(&coefficients[start..]);
            unsafe { vld1q_f64(tail.as_ptr()) }
        };
        partials = vfmaq_f64(packed_coefficients, partials, t_vector);
    }
    let weights = [1.0_f64, x];
    let weighted = vmulq_f64(partials, unsafe { vld1q_f64(weights.as_ptr()) });
    let mut lanes = [0.0_f64; F64_LANES];
    unsafe { vst1q_f64(lanes.as_mut_ptr(), weighted) };
    lanes.into_iter().sum()
}

#[inline]
#[target_feature(enable = "neon")]
pub unsafe fn eval_horizontal_estrin_f64<const N: usize>(x: f64, coefficients: &[f64; N]) -> f64 {
    if N == 0 {
        return 0.0;
    }
    let mut scratch = *coefficients;
    let mut len = N;
    let mut x_power = x;
    while len > 1 {
        let next_len = len.div_ceil(2);
        let mut output = 0;
        while output < next_len {
            let lane_count = (next_len - output).min(F64_LANES);
            let mut even = [0.0_f64; F64_LANES];
            let mut odd = [0.0_f64; F64_LANES];
            for lane in 0..lane_count {
                let input = 2 * (output + lane);
                even[lane] = scratch[input];
                if input + 1 < len {
                    odd[lane] = scratch[input + 1];
                }
            }
            let even = unsafe { vld1q_f64(even.as_ptr()) };
            let odd = unsafe { vld1q_f64(odd.as_ptr()) };
            let combined = vfmaq_f64(even, odd, vdupq_n_f64(x_power));
            let mut result = [0.0_f64; F64_LANES];
            unsafe { vst1q_f64(result.as_mut_ptr(), combined) };
            scratch[output..output + lane_count].copy_from_slice(&result[..lane_count]);
            output += lane_count;
        }
        len = next_len;
        x_power *= x_power;
    }
    scratch[0]
}

#[inline]
#[target_feature(enable = "neon")]
pub unsafe fn eval_horner_f32<const N: usize>(values: &mut [f32], coefficients: &[f32; N]) {
    let packed_len = values.len() / F32_LANES * F32_LANES;
    for chunk in values[..packed_len].chunks_exact_mut(F32_LANES) {
        let x = unsafe { vld1q_f32(chunk.as_ptr()) };
        let mut p = vdupq_n_f32(0.0);
        for &coefficient in coefficients.iter().rev() {
            p = vfmaq_f32(vdupq_n_f32(coefficient), p, x);
        }
        unsafe { vst1q_f32(chunk.as_mut_ptr(), p) };
    }
    super::scalar_horner_f32(&mut values[packed_len..], coefficients);
}

#[inline]
#[target_feature(enable = "neon")]
pub unsafe fn eval_estrin_f32<const N: usize>(values: &mut [f32], coefficients: &[f32; N]) {
    let packed_len = values.len() / F32_LANES * F32_LANES;
    for chunk in values[..packed_len].chunks_exact_mut(F32_LANES) {
        let x = unsafe { vld1q_f32(chunk.as_ptr()) };
        let mut scratch = [vdupq_n_f32(0.0); N];
        for (slot, &coefficient) in scratch.iter_mut().zip(coefficients) {
            *slot = vdupq_n_f32(coefficient);
        }

        let mut len = N;
        let mut x_pow = x;
        while len > 1 {
            let mut write = 0;
            let mut read = 0;
            while read + 1 < len {
                scratch[write] = vfmaq_f32(scratch[read], scratch[read + 1], x_pow);
                write += 1;
                read += 2;
            }
            if read < len {
                scratch[write] = scratch[read];
                write += 1;
            }
            len = write;
            x_pow = vmulq_f32(x_pow, x_pow);
        }

        let result = if N == 0 { vdupq_n_f32(0.0) } else { scratch[0] };
        unsafe { vst1q_f32(chunk.as_mut_ptr(), result) };
    }
    super::scalar_estrin_f32(&mut values[packed_len..], coefficients);
}

#[inline]
#[target_feature(enable = "neon")]
pub unsafe fn eval_horner_f64<const N: usize>(values: &mut [f64], coefficients: &[f64; N]) {
    let packed_len = values.len() / F64_LANES * F64_LANES;
    for chunk in values[..packed_len].chunks_exact_mut(F64_LANES) {
        let x = unsafe { vld1q_f64(chunk.as_ptr()) };
        let mut p = vdupq_n_f64(0.0);
        for &coefficient in coefficients.iter().rev() {
            p = vfmaq_f64(vdupq_n_f64(coefficient), p, x);
        }
        unsafe { vst1q_f64(chunk.as_mut_ptr(), p) };
    }
    super::scalar_horner_f64(&mut values[packed_len..], coefficients);
}

#[inline]
#[target_feature(enable = "neon")]
pub unsafe fn eval_estrin_f64<const N: usize>(values: &mut [f64], coefficients: &[f64; N]) {
    let packed_len = values.len() / F64_LANES * F64_LANES;
    for chunk in values[..packed_len].chunks_exact_mut(F64_LANES) {
        let x = unsafe { vld1q_f64(chunk.as_ptr()) };
        let mut scratch = [vdupq_n_f64(0.0); N];
        for (slot, &coefficient) in scratch.iter_mut().zip(coefficients) {
            *slot = vdupq_n_f64(coefficient);
        }

        let mut len = N;
        let mut x_pow = x;
        while len > 1 {
            let mut write = 0;
            let mut read = 0;
            while read + 1 < len {
                scratch[write] = vfmaq_f64(scratch[read], scratch[read + 1], x_pow);
                write += 1;
                read += 2;
            }
            if read < len {
                scratch[write] = scratch[read];
                write += 1;
            }
            len = write;
            x_pow = vmulq_f64(x_pow, x_pow);
        }

        let result = if N == 0 { vdupq_n_f64(0.0) } else { scratch[0] };
        unsafe { vst1q_f64(chunk.as_mut_ptr(), result) };
    }
    super::scalar_estrin_f64(&mut values[packed_len..], coefficients);
}
