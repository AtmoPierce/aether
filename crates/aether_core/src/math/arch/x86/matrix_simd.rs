use crate::math::{Matrix, Vector};

#[cfg(target_arch = "x86")]
use core::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

#[target_feature(enable = "avx,fma")]
pub unsafe fn mul_vec6_avx_fma_f64<const M: usize>(
    matrix: &Matrix<f64, M, 6>,
    rhs: &Vector<f64, 6>,
) -> Vector<f64, M> {
    let mut result = Vector { data: [0.0; M] };
    let v_ptr = rhs.data.as_ptr();

    let v0123 = unsafe { _mm256_loadu_pd(v_ptr.add(0)) };
    let v45 = unsafe { _mm_loadu_pd(v_ptr.add(4)) };

    for i in 0..M {
        let row_ptr = matrix.data[i].as_ptr();
        let r0123 = unsafe { _mm256_loadu_pd(row_ptr.add(0)) };
        let r45 = unsafe { _mm_loadu_pd(row_ptr.add(4)) };

        let p0123 = _mm256_fmadd_pd(r0123, v0123, _mm256_setzero_pd());
        let p45 = _mm_mul_pd(r45, v45);

        let mut lanes4 = [0.0_f64; 4];
        unsafe { _mm256_storeu_pd(lanes4.as_mut_ptr(), p0123) };
        let mut lanes2 = [0.0_f64; 2];
        unsafe { _mm_storeu_pd(lanes2.as_mut_ptr(), p45) };

        result.data[i] = lanes4[0] + lanes4[1] + lanes4[2] + lanes4[3] + lanes2[0] + lanes2[1];
    }

    result
}

#[target_feature(enable = "avx")]
pub unsafe fn mul_vec6_avx_f64<const M: usize>(
    matrix: &Matrix<f64, M, 6>,
    rhs: &Vector<f64, 6>,
) -> Vector<f64, M> {
    let mut result = Vector { data: [0.0; M] };
    let v_ptr = rhs.data.as_ptr();

    let v0123 = unsafe { _mm256_loadu_pd(v_ptr.add(0)) };
    let v45 = unsafe { _mm_loadu_pd(v_ptr.add(4)) };

    for i in 0..M {
        let row_ptr = matrix.data[i].as_ptr();
        let r0123 = unsafe { _mm256_loadu_pd(row_ptr.add(0)) };
        let r45 = unsafe { _mm_loadu_pd(row_ptr.add(4)) };

        let p0123 = _mm256_mul_pd(r0123, v0123);
        let p45 = _mm_mul_pd(r45, v45);

        let mut lanes4 = [0.0_f64; 4];
        unsafe { _mm256_storeu_pd(lanes4.as_mut_ptr(), p0123) };
        let mut lanes2 = [0.0_f64; 2];
        unsafe { _mm_storeu_pd(lanes2.as_mut_ptr(), p45) };

        result.data[i] = lanes4[0] + lanes4[1] + lanes4[2] + lanes4[3] + lanes2[0] + lanes2[1];
    }

    result
}

#[target_feature(enable = "sse2,fma")]
pub unsafe fn mul_vec6_fma128_f64<const M: usize>(
    matrix: &Matrix<f64, M, 6>,
    rhs: &Vector<f64, 6>,
) -> Vector<f64, M> {
    let mut result = Vector { data: [0.0; M] };
    let v_ptr = rhs.data.as_ptr();

    let v01 = unsafe { _mm_loadu_pd(v_ptr.add(0)) };
    let v23 = unsafe { _mm_loadu_pd(v_ptr.add(2)) };
    let v45 = unsafe { _mm_loadu_pd(v_ptr.add(4)) };

    for i in 0..M {
        let row_ptr = matrix.data[i].as_ptr();
        let r01 = unsafe { _mm_loadu_pd(row_ptr.add(0)) };
        let r23 = unsafe { _mm_loadu_pd(row_ptr.add(2)) };
        let r45 = unsafe { _mm_loadu_pd(row_ptr.add(4)) };

        let mut acc = _mm_setzero_pd();
        acc = _mm_fmadd_pd(r01, v01, acc);
        acc = _mm_fmadd_pd(r23, v23, acc);
        acc = _mm_fmadd_pd(r45, v45, acc);

        let mut lanes = [0.0_f64; 2];
        unsafe { _mm_storeu_pd(lanes.as_mut_ptr(), acc) };
        result.data[i] = lanes[0] + lanes[1];
    }

    result
}

#[target_feature(enable = "sse2")]
pub unsafe fn mul_vec6_sse2_f64<const M: usize>(
    matrix: &Matrix<f64, M, 6>,
    rhs: &Vector<f64, 6>,
) -> Vector<f64, M> {
    let mut result = Vector { data: [0.0; M] };
    let v_ptr = rhs.data.as_ptr();

    let v01 = unsafe { _mm_loadu_pd(v_ptr.add(0)) };
    let v23 = unsafe { _mm_loadu_pd(v_ptr.add(2)) };
    let v45 = unsafe { _mm_loadu_pd(v_ptr.add(4)) };

    for i in 0..M {
        let row_ptr = matrix.data[i].as_ptr();
        let r01 = unsafe { _mm_loadu_pd(row_ptr.add(0)) };
        let r23 = unsafe { _mm_loadu_pd(row_ptr.add(2)) };
        let r45 = unsafe { _mm_loadu_pd(row_ptr.add(4)) };

        let p01 = _mm_mul_pd(r01, v01);
        let p23 = _mm_mul_pd(r23, v23);
        let p45 = _mm_mul_pd(r45, v45);
        let sum = _mm_add_pd(_mm_add_pd(p01, p23), p45);

        let mut lanes = [0.0_f64; 2];
        unsafe { _mm_storeu_pd(lanes.as_mut_ptr(), sum) };
        result.data[i] = lanes[0] + lanes[1];
    }

    result
}

#[target_feature(enable = "avx,fma")]
pub unsafe fn mul_vec4_avx_fma_f64<const M: usize>(
    matrix: &Matrix<f64, M, 4>,
    rhs: &Vector<f64, 4>,
) -> Vector<f64, M> {
    let mut result = Vector { data: [0.0; M] };
    let v = unsafe { _mm256_loadu_pd(rhs.data.as_ptr()) };

    for i in 0..M {
        let r = unsafe { _mm256_loadu_pd(matrix.data[i].as_ptr()) };
        let p = _mm256_fmadd_pd(r, v, _mm256_setzero_pd());
        let mut lanes = [0.0_f64; 4];
        unsafe { _mm256_storeu_pd(lanes.as_mut_ptr(), p) };
        result.data[i] = lanes[0] + lanes[1] + lanes[2] + lanes[3];
    }
    result
}

#[target_feature(enable = "avx")]
pub unsafe fn mul_vec4_avx_f64<const M: usize>(
    matrix: &Matrix<f64, M, 4>,
    rhs: &Vector<f64, 4>,
) -> Vector<f64, M> {
    let mut result = Vector { data: [0.0; M] };
    let v = unsafe { _mm256_loadu_pd(rhs.data.as_ptr()) };

    for i in 0..M {
        let r = unsafe { _mm256_loadu_pd(matrix.data[i].as_ptr()) };
        let p = _mm256_mul_pd(r, v);
        let mut lanes = [0.0_f64; 4];
        unsafe { _mm256_storeu_pd(lanes.as_mut_ptr(), p) };
        result.data[i] = lanes[0] + lanes[1] + lanes[2] + lanes[3];
    }
    result
}

#[target_feature(enable = "avx,fma")]
pub unsafe fn mul_vec6_avx_fma_f32<const M: usize>(
    matrix: &Matrix<f32, M, 6>,
    rhs: &Vector<f32, 6>,
) -> Vector<f32, M> {
    let mut result = Vector { data: [0.0; M] };
    let v = _mm256_set_ps(
        0.0,
        0.0,
        rhs.data[5],
        rhs.data[4],
        rhs.data[3],
        rhs.data[2],
        rhs.data[1],
        rhs.data[0],
    );

    for i in 0..M {
        let row = &matrix.data[i];
        let r = _mm256_set_ps(0.0, 0.0, row[5], row[4], row[3], row[2], row[1], row[0]);
        let p = _mm256_fmadd_ps(r, v, _mm256_setzero_ps());

        let mut lanes = [0.0_f32; 8];
        unsafe { _mm256_storeu_ps(lanes.as_mut_ptr(), p) };
        result.data[i] = lanes[0] + lanes[1] + lanes[2] + lanes[3] + lanes[4] + lanes[5];
    }

    result
}

#[target_feature(enable = "avx")]
pub unsafe fn mul_vec6_avx_f32<const M: usize>(
    matrix: &Matrix<f32, M, 6>,
    rhs: &Vector<f32, 6>,
) -> Vector<f32, M> {
    let mut result = Vector { data: [0.0; M] };
    let v = _mm256_set_ps(
        0.0,
        0.0,
        rhs.data[5],
        rhs.data[4],
        rhs.data[3],
        rhs.data[2],
        rhs.data[1],
        rhs.data[0],
    );

    for i in 0..M {
        let row = &matrix.data[i];
        let r = _mm256_set_ps(0.0, 0.0, row[5], row[4], row[3], row[2], row[1], row[0]);
        let p = _mm256_mul_ps(r, v);

        let mut lanes = [0.0_f32; 8];
        unsafe { _mm256_storeu_ps(lanes.as_mut_ptr(), p) };
        result.data[i] = lanes[0] + lanes[1] + lanes[2] + lanes[3] + lanes[4] + lanes[5];
    }

    result
}

#[target_feature(enable = "sse,fma")]
pub unsafe fn mul_vec6_fma128_f32<const M: usize>(
    matrix: &Matrix<f32, M, 6>,
    rhs: &Vector<f32, 6>,
) -> Vector<f32, M> {
    let mut result = Vector { data: [0.0; M] };

    let v0123 = unsafe { _mm_loadu_ps(rhs.data.as_ptr()) };

    for i in 0..M {
        let row = &matrix.data[i];
        let r0123 = unsafe { _mm_loadu_ps(row.as_ptr()) };
        let mut acc = _mm_setzero_ps();
        acc = _mm_fmadd_ps(r0123, v0123, acc);

        let mut lanes = [0.0_f32; 4];
        unsafe { _mm_storeu_ps(lanes.as_mut_ptr(), acc) };
        result.data[i] = lanes[0] + lanes[1] + lanes[2] + lanes[3] + row[4] * rhs.data[4] + row[5] * rhs.data[5];
    }

    result
}

#[target_feature(enable = "sse")]
pub unsafe fn mul_vec6_sse_f32<const M: usize>(
    matrix: &Matrix<f32, M, 6>,
    rhs: &Vector<f32, 6>,
) -> Vector<f32, M> {
    let mut result = Vector { data: [0.0; M] };

    let v0123 = unsafe { _mm_loadu_ps(rhs.data.as_ptr()) };

    for i in 0..M {
        let row = &matrix.data[i];
        let r0123 = unsafe { _mm_loadu_ps(row.as_ptr()) };
        let p = _mm_mul_ps(r0123, v0123);

        let mut lanes = [0.0_f32; 4];
        unsafe { _mm_storeu_ps(lanes.as_mut_ptr(), p) };
        result.data[i] = lanes[0] + lanes[1] + lanes[2] + lanes[3] + row[4] * rhs.data[4] + row[5] * rhs.data[5];
    }

    result
}

#[target_feature(enable = "sse,fma")]
pub unsafe fn mul_vec4_sse_fma_f32<const M: usize>(
    matrix: &Matrix<f32, M, 4>,
    rhs: &Vector<f32, 4>,
) -> Vector<f32, M> {
    let mut result = Vector { data: [0.0; M] };
    let v = unsafe { _mm_loadu_ps(rhs.data.as_ptr()) };

    for i in 0..M {
        let r = unsafe { _mm_loadu_ps(matrix.data[i].as_ptr()) };
        let p = _mm_fmadd_ps(r, v, _mm_setzero_ps());
        let mut lanes = [0.0_f32; 4];
        unsafe { _mm_storeu_ps(lanes.as_mut_ptr(), p) };
        result.data[i] = lanes[0] + lanes[1] + lanes[2] + lanes[3];
    }
    result
}

#[target_feature(enable = "sse")]
pub unsafe fn mul_vec4_sse_f32<const M: usize>(
    matrix: &Matrix<f32, M, 4>,
    rhs: &Vector<f32, 4>,
) -> Vector<f32, M> {
    let mut result = Vector { data: [0.0; M] };
    let v = unsafe { _mm_loadu_ps(rhs.data.as_ptr()) };

    for i in 0..M {
        let r = unsafe { _mm_loadu_ps(matrix.data[i].as_ptr()) };
        let p = _mm_mul_ps(r, v);
        let mut lanes = [0.0_f32; 4];
        unsafe { _mm_storeu_ps(lanes.as_mut_ptr(), p) };
        result.data[i] = lanes[0] + lanes[1] + lanes[2] + lanes[3];
    }
    result
}
