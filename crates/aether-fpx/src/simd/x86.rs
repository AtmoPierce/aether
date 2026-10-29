#[cfg(target_arch = "x86")]
use core::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

const F32_LANES: usize = 8;
const F64_LANES: usize = 4;

#[inline]
#[target_feature(enable = "avx2,fma")]
pub unsafe fn eval_blocked_horner_f32<const N: usize>(x: f32, coefficients: &[f32; N]) -> f32 {
    if N == 0 {
        return 0.0;
    }
    let mut t = 1.0_f32;
    for _ in 0..F32_LANES {
        t *= x;
    }
    let t_vector = _mm256_set1_ps(t);
    let mut partials = _mm256_setzero_ps();
    let block_count = N.div_ceil(F32_LANES);
    for block in (0..block_count).rev() {
        let start = block * F32_LANES;
        let packed_coefficients = if start + F32_LANES <= N {
            unsafe { _mm256_loadu_ps(coefficients.as_ptr().add(start)) }
        } else {
            let mut tail = [0.0_f32; F32_LANES];
            tail[..N - start].copy_from_slice(&coefficients[start..]);
            unsafe { _mm256_loadu_ps(tail.as_ptr()) }
        };
        partials = _mm256_fmadd_ps(partials, t_vector, packed_coefficients);
    }

    let mut powers = [1.0_f32; F32_LANES];
    for lane in 1..F32_LANES {
        powers[lane] = powers[lane - 1] * x;
    }
    let weights = unsafe { _mm256_loadu_ps(powers.as_ptr()) };
    let weighted = _mm256_mul_ps(partials, weights);
    let mut lanes = [0.0_f32; F32_LANES];
    unsafe { _mm256_storeu_ps(lanes.as_mut_ptr(), weighted) };
    lanes.into_iter().sum()
}

#[inline]
#[target_feature(enable = "avx2,fma")]
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
            let even = unsafe { _mm256_loadu_ps(even.as_ptr()) };
            let odd = unsafe { _mm256_loadu_ps(odd.as_ptr()) };
            let combined = _mm256_fmadd_ps(odd, _mm256_set1_ps(x_power), even);
            let mut result = [0.0_f32; F32_LANES];
            unsafe { _mm256_storeu_ps(result.as_mut_ptr(), combined) };
            scratch[output..output + lane_count].copy_from_slice(&result[..lane_count]);
            output += lane_count;
        }
        len = next_len;
        x_power *= x_power;
    }
    scratch[0]
}

#[inline]
#[target_feature(enable = "avx2,fma")]
pub unsafe fn eval_blocked_horner_f64<const N: usize>(x: f64, coefficients: &[f64; N]) -> f64 {
    if N == 0 {
        return 0.0;
    }
    let mut t = 1.0_f64;
    for _ in 0..F64_LANES {
        t *= x;
    }
    let t_vector = _mm256_set1_pd(t);
    let mut partials = _mm256_setzero_pd();
    let block_count = N.div_ceil(F64_LANES);
    for block in (0..block_count).rev() {
        let start = block * F64_LANES;
        let packed_coefficients = if start + F64_LANES <= N {
            unsafe { _mm256_loadu_pd(coefficients.as_ptr().add(start)) }
        } else {
            let mut tail = [0.0_f64; F64_LANES];
            tail[..N - start].copy_from_slice(&coefficients[start..]);
            unsafe { _mm256_loadu_pd(tail.as_ptr()) }
        };
        partials = _mm256_fmadd_pd(partials, t_vector, packed_coefficients);
    }

    let mut powers = [1.0_f64; F64_LANES];
    for lane in 1..F64_LANES {
        powers[lane] = powers[lane - 1] * x;
    }
    let weights = unsafe { _mm256_loadu_pd(powers.as_ptr()) };
    let weighted = _mm256_mul_pd(partials, weights);
    let mut lanes = [0.0_f64; F64_LANES];
    unsafe { _mm256_storeu_pd(lanes.as_mut_ptr(), weighted) };
    lanes.into_iter().sum()
}

#[inline]
#[target_feature(enable = "avx2,fma")]
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
            let even = unsafe { _mm256_loadu_pd(even.as_ptr()) };
            let odd = unsafe { _mm256_loadu_pd(odd.as_ptr()) };
            let combined = _mm256_fmadd_pd(odd, _mm256_set1_pd(x_power), even);
            let mut result = [0.0_f64; F64_LANES];
            unsafe { _mm256_storeu_pd(result.as_mut_ptr(), combined) };
            scratch[output..output + lane_count].copy_from_slice(&result[..lane_count]);
            output += lane_count;
        }
        len = next_len;
        x_power *= x_power;
    }
    scratch[0]
}

#[inline]
#[target_feature(enable = "avx2,fma")]
pub unsafe fn eval_horner_f32<const N: usize>(values: &mut [f32], coefficients: &[f32; N]) {
    let packed_len = values.len() / F32_LANES * F32_LANES;
    for chunk in values[..packed_len].chunks_exact_mut(F32_LANES) {
        let x = unsafe { _mm256_loadu_ps(chunk.as_ptr()) };
        let mut p = _mm256_setzero_ps();
        for &coefficient in coefficients.iter().rev() {
            p = _mm256_fmadd_ps(p, x, _mm256_set1_ps(coefficient));
        }
        unsafe { _mm256_storeu_ps(chunk.as_mut_ptr(), p) };
    }
    super::scalar_horner_f32(&mut values[packed_len..], coefficients);
}

#[inline]
#[target_feature(enable = "avx2,fma")]
pub unsafe fn eval_estrin_f32<const N: usize>(values: &mut [f32], coefficients: &[f32; N]) {
    let packed_len = values.len() / F32_LANES * F32_LANES;
    for chunk in values[..packed_len].chunks_exact_mut(F32_LANES) {
        let x = unsafe { _mm256_loadu_ps(chunk.as_ptr()) };
        let mut scratch = [_mm256_setzero_ps(); N];
        for (slot, &coefficient) in scratch.iter_mut().zip(coefficients) {
            *slot = _mm256_set1_ps(coefficient);
        }

        let mut len = N;
        let mut x_pow = x;
        while len > 1 {
            let mut write = 0;
            let mut read = 0;
            while read + 1 < len {
                scratch[write] = _mm256_fmadd_ps(scratch[read + 1], x_pow, scratch[read]);
                write += 1;
                read += 2;
            }
            if read < len {
                scratch[write] = scratch[read];
                write += 1;
            }
            len = write;
            x_pow = _mm256_mul_ps(x_pow, x_pow);
        }

        let result = if N == 0 {
            _mm256_setzero_ps()
        } else {
            scratch[0]
        };
        unsafe { _mm256_storeu_ps(chunk.as_mut_ptr(), result) };
    }
    super::scalar_estrin_f32(&mut values[packed_len..], coefficients);
}

#[inline]
#[target_feature(enable = "avx2,fma")]
pub unsafe fn eval_horner_f64<const N: usize>(values: &mut [f64], coefficients: &[f64; N]) {
    let packed_len = values.len() / F64_LANES * F64_LANES;
    for chunk in values[..packed_len].chunks_exact_mut(F64_LANES) {
        let x = unsafe { _mm256_loadu_pd(chunk.as_ptr()) };
        let mut p = _mm256_setzero_pd();
        for &coefficient in coefficients.iter().rev() {
            p = _mm256_fmadd_pd(p, x, _mm256_set1_pd(coefficient));
        }
        unsafe { _mm256_storeu_pd(chunk.as_mut_ptr(), p) };
    }
    super::scalar_horner_f64(&mut values[packed_len..], coefficients);
}

#[inline]
#[target_feature(enable = "avx2,fma")]
pub unsafe fn eval_estrin_f64<const N: usize>(values: &mut [f64], coefficients: &[f64; N]) {
    let packed_len = values.len() / F64_LANES * F64_LANES;
    for chunk in values[..packed_len].chunks_exact_mut(F64_LANES) {
        let x = unsafe { _mm256_loadu_pd(chunk.as_ptr()) };
        let mut scratch = [_mm256_setzero_pd(); N];
        for (slot, &coefficient) in scratch.iter_mut().zip(coefficients) {
            *slot = _mm256_set1_pd(coefficient);
        }

        let mut len = N;
        let mut x_pow = x;
        while len > 1 {
            let mut write = 0;
            let mut read = 0;
            while read + 1 < len {
                scratch[write] = _mm256_fmadd_pd(scratch[read + 1], x_pow, scratch[read]);
                write += 1;
                read += 2;
            }
            if read < len {
                scratch[write] = scratch[read];
                write += 1;
            }
            len = write;
            x_pow = _mm256_mul_pd(x_pow, x_pow);
        }

        let result = if N == 0 {
            _mm256_setzero_pd()
        } else {
            scratch[0]
        };
        unsafe { _mm256_storeu_pd(chunk.as_mut_ptr(), result) };
    }
    super::scalar_estrin_f64(&mut values[packed_len..], coefficients);
}
