use crate::math::{Matrix, Vector};
use core::arch::aarch64::*;

pub const NEON_LANES_U8: usize = 16;
pub const NEON_LANES_U16: usize = 8;
pub const NEON_LANES_U32: usize = 4;
pub const NEON_LANES_U64: usize = 2;
pub const NEON_LANES_F32: usize = NEON_LANES_U32;
pub const NEON_LANES_F64: usize = NEON_LANES_U64;

#[target_feature(enable = "neon")]
pub unsafe fn dot4_neon_f32(a: &Vector<f32, 4>, b: &Vector<f32, 4>) -> f32 {
    let va = unsafe { vld1q_f32(a.data.as_ptr()) };
    let vb = unsafe { vld1q_f32(b.data.as_ptr()) };
    let prod = vmulq_f32(va, vb);
    vaddvq_f32(prod)
}

#[target_feature(enable = "neon")]
pub unsafe fn dot4_neon_f64(a: &Vector<f64, 4>, b: &Vector<f64, 4>) -> f64 {
    let va0 = unsafe { vld1q_f64(a.data.as_ptr()) };
    let vb0 = unsafe { vld1q_f64(b.data.as_ptr()) };
    let va1 = unsafe { vld1q_f64(a.data.as_ptr().add(NEON_LANES_F64)) };
    let vb1 = unsafe { vld1q_f64(b.data.as_ptr().add(NEON_LANES_F64)) };
    vaddvq_f64(vmulq_f64(va0, vb0)) + vaddvq_f64(vmulq_f64(va1, vb1))
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_vec4_neon_f32<const M: usize>(
    matrix: &Matrix<f32, M, 4>,
    rhs: &Vector<f32, 4>,
) -> Vector<f32, M> {
    let mut out = Vector { data: [0.0; M] };
    let v = unsafe { vld1q_f32(rhs.data.as_ptr()) };

    for i in 0..M {
        let r = unsafe { vld1q_f32(matrix.data[i].as_ptr()) };
        out.data[i] = vaddvq_f32(vmulq_f32(r, v));
    }

    out
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_vec4_neon_f64<const M: usize>(
    matrix: &Matrix<f64, M, 4>,
    rhs: &Vector<f64, 4>,
) -> Vector<f64, M> {
    let mut out = Vector { data: [0.0; M] };

    let v0 = unsafe { vld1q_f64(rhs.data.as_ptr()) };
    let v1 = unsafe { vld1q_f64(rhs.data.as_ptr().add(NEON_LANES_F64)) };

    for i in 0..M {
        let r0 = unsafe { vld1q_f64(matrix.data[i].as_ptr()) };
        let r1 = unsafe { vld1q_f64(matrix.data[i].as_ptr().add(NEON_LANES_F64)) };
        out.data[i] = vaddvq_f64(vmulq_f64(r0, v0)) + vaddvq_f64(vmulq_f64(r1, v1));
    }

    out
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_matrix_neon_f32<const M: usize, const N: usize, const P: usize>(
    lhs: &Matrix<f32, M, N>,
    rhs: &Matrix<f32, N, P>,
) -> Matrix<f32, M, P> {
    let mut out = Matrix { data: [[0.0; P]; M] };

    for i in 0..M {
        for j in 0..P {
            let mut acc = 0.0_f32;
            let mut k = 0;

            while k + NEON_LANES_F32 <= N {
                let l = unsafe { vld1q_f32(lhs.data[i].as_ptr().add(k)) };
                let r_col = [
                    rhs.data[k][j],
                    rhs.data[k + 1][j],
                    rhs.data[k + 2][j],
                    rhs.data[k + 3][j],
                ];
                let r = unsafe { vld1q_f32(r_col.as_ptr()) };
                acc += vaddvq_f32(vmulq_f32(l, r));
                k += NEON_LANES_F32;
            }

            while k < N {
                acc += lhs.data[i][k] * rhs.data[k][j];
                k += 1;
            }

            out.data[i][j] = acc;
        }
    }

    out
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_matrix_neon_f64<const M: usize, const N: usize, const P: usize>(
    lhs: &Matrix<f64, M, N>,
    rhs: &Matrix<f64, N, P>,
) -> Matrix<f64, M, P> {
    let mut out = Matrix { data: [[0.0; P]; M] };

    for i in 0..M {
        for j in 0..P {
            let mut acc = 0.0_f64;
            let mut k = 0;

            while k + NEON_LANES_F64 <= N {
                let l = unsafe { vld1q_f64(lhs.data[i].as_ptr().add(k)) };
                let r_col = [rhs.data[k][j], rhs.data[k + 1][j]];
                let r = unsafe { vld1q_f64(r_col.as_ptr()) };
                acc += vaddvq_f64(vmulq_f64(l, r));
                k += NEON_LANES_F64;
            }

            while k < N {
                acc += lhs.data[i][k] * rhs.data[k][j];
                k += 1;
            }

            out.data[i][j] = acc;
        }
    }

    out
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_mat3_neon_f32(
    lhs: &Matrix<f32, 3, 3>,
    rhs: &Matrix<f32, 3, 3>,
) -> Matrix<f32, 3, 3> {
    unsafe { mul_matrix_neon_f32(lhs, rhs) }
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_mat3_neon_f64(
    lhs: &Matrix<f64, 3, 3>,
    rhs: &Matrix<f64, 3, 3>,
) -> Matrix<f64, 3, 3> {
    unsafe { mul_matrix_neon_f64(lhs, rhs) }
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_mat4_neon_f32(
    lhs: &Matrix<f32, 4, 4>,
    rhs: &Matrix<f32, 4, 4>,
) -> Matrix<f32, 4, 4> {
    unsafe { mul_matrix_neon_f32(lhs, rhs) }
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_mat4_neon_f64(
    lhs: &Matrix<f64, 4, 4>,
    rhs: &Matrix<f64, 4, 4>,
) -> Matrix<f64, 4, 4> {
    unsafe { mul_matrix_neon_f64(lhs, rhs) }
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_mat6_neon_f32(
    lhs: &Matrix<f32, 6, 6>,
    rhs: &Matrix<f32, 6, 6>,
) -> Matrix<f32, 6, 6> {
    unsafe { mul_matrix_neon_f32(lhs, rhs) }
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_mat6_neon_f64(
    lhs: &Matrix<f64, 6, 6>,
    rhs: &Matrix<f64, 6, 6>,
) -> Matrix<f64, 6, 6> {
    unsafe { mul_matrix_neon_f64(lhs, rhs) }
}
