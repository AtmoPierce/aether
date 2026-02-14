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
pub unsafe fn mul_vec6_neon_f32<const M: usize>(
    matrix: &Matrix<f32, M, 6>,
    rhs: &Vector<f32, 6>,
) -> Vector<f32, M> {
    let mut out = Vector { data: [0.0; M] };
    let v0123 = unsafe { vld1q_f32(rhs.data.as_ptr()) };
    let v45xx = [rhs.data[4], rhs.data[5], 0.0_f32, 0.0_f32];
    let v45xx_vec = unsafe { vld1q_f32(v45xx.as_ptr()) };

    for i in 0..M {
        let r0123 = unsafe { vld1q_f32(matrix.data[i].as_ptr()) };
        let r45xx = [matrix.data[i][4], matrix.data[i][5], 0.0_f32, 0.0_f32];
        let r45xx_vec = unsafe { vld1q_f32(r45xx.as_ptr()) };
        out.data[i] = vaddvq_f32(vmulq_f32(r0123, v0123)) + vaddvq_f32(vmulq_f32(r45xx_vec, v45xx_vec));
    }

    out
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_vec6_neon_f64<const M: usize>(
    matrix: &Matrix<f64, M, 6>,
    rhs: &Vector<f64, 6>,
) -> Vector<f64, M> {
    let mut out = Vector { data: [0.0; M] };

    let v01 = unsafe { vld1q_f64(rhs.data.as_ptr()) };
    let v23 = unsafe { vld1q_f64(rhs.data.as_ptr().add(2)) };
    let v45 = unsafe { vld1q_f64(rhs.data.as_ptr().add(4)) };

    for i in 0..M {
        let r01 = unsafe { vld1q_f64(matrix.data[i].as_ptr()) };
        let r23 = unsafe { vld1q_f64(matrix.data[i].as_ptr().add(2)) };
        let r45 = unsafe { vld1q_f64(matrix.data[i].as_ptr().add(4)) };

        out.data[i] = vaddvq_f64(vmulq_f64(r01, v01))
            + vaddvq_f64(vmulq_f64(r23, v23))
            + vaddvq_f64(vmulq_f64(r45, v45));
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
    let mut out = Matrix { data: [[0.0_f32; 3]; 3] };

    for i in 0..3 {
        let mut acc = vdupq_n_f32(0.0_f32);
        for k in 0..3 {
            let b = [rhs.data[k][0], rhs.data[k][1], rhs.data[k][2], 0.0_f32];
            let b_vec = unsafe { vld1q_f32(b.as_ptr()) };
            acc = vmlaq_n_f32(acc, b_vec, lhs.data[i][k]);
        }
        let mut lanes = [0.0_f32; 4];
        unsafe { vst1q_f32(lanes.as_mut_ptr(), acc) };
        out.data[i][0] = lanes[0];
        out.data[i][1] = lanes[1];
        out.data[i][2] = lanes[2];
    }

    out
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_mat3_neon_f64(
    lhs: &Matrix<f64, 3, 3>,
    rhs: &Matrix<f64, 3, 3>,
) -> Matrix<f64, 3, 3> {
    let mut out = Matrix { data: [[0.0_f64; 3]; 3] };

    for i in 0..3 {
        let mut acc01 = vdupq_n_f64(0.0_f64);
        let mut acc2x = vdupq_n_f64(0.0_f64);
        for k in 0..3 {
            let b01 = [rhs.data[k][0], rhs.data[k][1]];
            let b2x = [rhs.data[k][2], 0.0_f64];
            let b01_vec = unsafe { vld1q_f64(b01.as_ptr()) };
            let b2x_vec = unsafe { vld1q_f64(b2x.as_ptr()) };
            acc01 = vfmaq_n_f64(acc01, b01_vec, lhs.data[i][k]);
            acc2x = vfmaq_n_f64(acc2x, b2x_vec, lhs.data[i][k]);
        }
        out.data[i][0] = vgetq_lane_f64(acc01, 0);
        out.data[i][1] = vgetq_lane_f64(acc01, 1);
        out.data[i][2] = vgetq_lane_f64(acc2x, 0);
    }

    out
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_mat4_neon_f32(
    lhs: &Matrix<f32, 4, 4>,
    rhs: &Matrix<f32, 4, 4>,
) -> Matrix<f32, 4, 4> {
    let mut out = Matrix { data: [[0.0_f32; 4]; 4] };

    for i in 0..4 {
        let mut acc = vdupq_n_f32(0.0_f32);
        for k in 0..4 {
            let b_vec = unsafe { vld1q_f32(rhs.data[k].as_ptr()) };
            acc = vmlaq_n_f32(acc, b_vec, lhs.data[i][k]);
        }
        unsafe { vst1q_f32(out.data[i].as_mut_ptr(), acc) };
    }

    out
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_mat4_neon_f64(
    lhs: &Matrix<f64, 4, 4>,
    rhs: &Matrix<f64, 4, 4>,
) -> Matrix<f64, 4, 4> {
    let mut out = Matrix { data: [[0.0_f64; 4]; 4] };

    for i in 0..4 {
        let mut acc01 = vdupq_n_f64(0.0_f64);
        let mut acc23 = vdupq_n_f64(0.0_f64);
        for k in 0..4 {
            let b01_vec = unsafe { vld1q_f64(rhs.data[k].as_ptr()) };
            let b23_vec = unsafe { vld1q_f64(rhs.data[k].as_ptr().add(2)) };
            acc01 = vfmaq_n_f64(acc01, b01_vec, lhs.data[i][k]);
            acc23 = vfmaq_n_f64(acc23, b23_vec, lhs.data[i][k]);
        }
        unsafe {
            vst1q_f64(out.data[i].as_mut_ptr(), acc01);
            vst1q_f64(out.data[i].as_mut_ptr().add(2), acc23);
        }
    }

    out
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_mat6_neon_f32(
    lhs: &Matrix<f32, 6, 6>,
    rhs: &Matrix<f32, 6, 6>,
) -> Matrix<f32, 6, 6> {
    let mut out = Matrix { data: [[0.0_f32; 6]; 6] };

    for i in 0..6 {
        let mut acc0123 = vdupq_n_f32(0.0_f32);
        let mut acc45xx = vdupq_n_f32(0.0_f32);
        for k in 0..6 {
            let b0123 = unsafe { vld1q_f32(rhs.data[k].as_ptr()) };
            let b45xx = [rhs.data[k][4], rhs.data[k][5], 0.0_f32, 0.0_f32];
            let b45xx_vec = unsafe { vld1q_f32(b45xx.as_ptr()) };
            acc0123 = vmlaq_n_f32(acc0123, b0123, lhs.data[i][k]);
            acc45xx = vmlaq_n_f32(acc45xx, b45xx_vec, lhs.data[i][k]);
        }

        unsafe {
            vst1q_f32(out.data[i].as_mut_ptr(), acc0123);
        }
        let mut tail = [0.0_f32; 4];
        unsafe { vst1q_f32(tail.as_mut_ptr(), acc45xx) };
        out.data[i][4] = tail[0];
        out.data[i][5] = tail[1];
    }

    out
}

#[target_feature(enable = "neon")]
pub unsafe fn mul_mat6_neon_f64(
    lhs: &Matrix<f64, 6, 6>,
    rhs: &Matrix<f64, 6, 6>,
) -> Matrix<f64, 6, 6> {
    let mut out = Matrix { data: [[0.0_f64; 6]; 6] };

    for i in 0..6 {
        let mut acc01 = vdupq_n_f64(0.0_f64);
        let mut acc23 = vdupq_n_f64(0.0_f64);
        let mut acc45 = vdupq_n_f64(0.0_f64);
        for k in 0..6 {
            let b01_vec = unsafe { vld1q_f64(rhs.data[k].as_ptr()) };
            let b23_vec = unsafe { vld1q_f64(rhs.data[k].as_ptr().add(2)) };
            let b45_vec = unsafe { vld1q_f64(rhs.data[k].as_ptr().add(4)) };
            acc01 = vfmaq_n_f64(acc01, b01_vec, lhs.data[i][k]);
            acc23 = vfmaq_n_f64(acc23, b23_vec, lhs.data[i][k]);
            acc45 = vfmaq_n_f64(acc45, b45_vec, lhs.data[i][k]);
        }
        unsafe {
            vst1q_f64(out.data[i].as_mut_ptr(), acc01);
            vst1q_f64(out.data[i].as_mut_ptr().add(2), acc23);
            vst1q_f64(out.data[i].as_mut_ptr().add(4), acc45);
        }
    }

    out
}
