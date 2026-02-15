use crate::math::{Matrix, Vector};

#[inline(always)]
pub fn mul_vec4_m33_f64<const M: usize>(
    matrix: &Matrix<f64, M, 4>,
    rhs: &Vector<f64, 4>,
) -> Vector<f64, M> {
    let mut out = Vector { data: [0.0; M] };
    for i in 0..M {
        let row = &matrix.data[i];
        out.data[i] = row[0] * rhs.data[0]
            + row[1] * rhs.data[1]
            + row[2] * rhs.data[2]
            + row[3] * rhs.data[3];
    }
    out
}

#[inline(always)]
pub fn mul_vec4_m33_f32<const M: usize>(
    matrix: &Matrix<f32, M, 4>,
    rhs: &Vector<f32, 4>,
) -> Vector<f32, M> {
    let mut out = Vector { data: [0.0; M] };
    for i in 0..M {
        let row = &matrix.data[i];
        out.data[i] = row[0] * rhs.data[0]
            + row[1] * rhs.data[1]
            + row[2] * rhs.data[2]
            + row[3] * rhs.data[3];
    }
    out
}

#[inline(always)]
pub fn dot4_m33_f64(a: &Vector<f64, 4>, b: &Vector<f64, 4>) -> f64 {
    a.data[0] * b.data[0] + a.data[1] * b.data[1] + a.data[2] * b.data[2] + a.data[3] * b.data[3]
}

#[inline(always)]
pub fn dot4_m33_f32(a: &Vector<f32, 4>, b: &Vector<f32, 4>) -> f32 {
    a.data[0] * b.data[0] + a.data[1] * b.data[1] + a.data[2] * b.data[2] + a.data[3] * b.data[3]
}
