use crate::math::{
    arch::x86::{matrix_simd, vector_simd},
    Matrix, Vector,
};

fn dot4_f32_scalar(a: &Vector<f32, 4>, b: &Vector<f32, 4>) -> f32 {
    a.data[0] * b.data[0] + a.data[1] * b.data[1] + a.data[2] * b.data[2] + a.data[3] * b.data[3]
}

fn dot4_f64_scalar(a: &Vector<f64, 4>, b: &Vector<f64, 4>) -> f64 {
    a.data[0] * b.data[0] + a.data[1] * b.data[1] + a.data[2] * b.data[2] + a.data[3] * b.data[3]
}

fn mul_vec4_scalar_f32<const M: usize>(
    matrix: &Matrix<f32, M, 4>,
    rhs: &Vector<f32, 4>,
) -> Vector<f32, M> {
    let mut out = Vector { data: [0.0; M] };
    for i in 0..M {
        out.data[i] = matrix.data[i][0] * rhs.data[0]
            + matrix.data[i][1] * rhs.data[1]
            + matrix.data[i][2] * rhs.data[2]
            + matrix.data[i][3] * rhs.data[3];
    }
    out
}

fn mul_vec4_scalar_f64<const M: usize>(
    matrix: &Matrix<f64, M, 4>,
    rhs: &Vector<f64, 4>,
) -> Vector<f64, M> {
    let mut out = Vector { data: [0.0; M] };
    for i in 0..M {
        out.data[i] = matrix.data[i][0] * rhs.data[0]
            + matrix.data[i][1] * rhs.data[1]
            + matrix.data[i][2] * rhs.data[2]
            + matrix.data[i][3] * rhs.data[3];
    }
    out
}

#[test]
fn dot4_sse_f32_matches_scalar() {
    if !std::is_x86_feature_detected!("sse") {
        return;
    }

    let a = Vector {
        data: [1.0_f32, -2.0, 3.5, 0.25],
    };
    let b = Vector {
        data: [2.0_f32, 1.5, -4.0, 10.0],
    };
    let expected = dot4_f32_scalar(&a, &b);
    let got = unsafe { vector_simd::dot4_sse_f32(&a, &b) };

    assert!((got - expected).abs() < 1.0e-6);
}

#[test]
fn dot4_avx_f64_matches_scalar() {
    if !std::is_x86_feature_detected!("avx") {
        return;
    }

    let a = Vector {
        data: [1.0_f64, -2.0, 3.5, 0.25],
    };
    let b = Vector {
        data: [2.0_f64, 1.5, -4.0, 10.0],
    };
    let expected = dot4_f64_scalar(&a, &b);
    let got = unsafe { vector_simd::dot4_avx_f64(&a, &b) };

    assert!((got - expected).abs() < 1.0e-12);
}

#[test]
fn dot4_avx_fma_f64_matches_scalar() {
    if !std::is_x86_feature_detected!("avx") || !std::is_x86_feature_detected!("fma") {
        return;
    }

    let a = Vector {
        data: [1.0_f64, -2.0, 3.5, 0.25],
    };
    let b = Vector {
        data: [2.0_f64, 1.5, -4.0, 10.0],
    };
    let expected = dot4_f64_scalar(&a, &b);
    let got = unsafe { vector_simd::dot4_avx_fma_f64(&a, &b) };

    assert!((got - expected).abs() < 1.0e-12);
}

#[test]
fn mul_vec4_sse_f32_matches_scalar() {
    if !std::is_x86_feature_detected!("sse") {
        return;
    }

    let matrix = Matrix {
        data: [
            [1.0_f32, 2.0, 3.0, 4.0],
            [-1.5, 0.25, 10.0, -2.0],
            [0.0, 0.5, -0.5, 2.0],
        ],
    };
    let rhs = Vector {
        data: [2.0_f32, -1.0, 0.5, 4.0],
    };

    let expected = mul_vec4_scalar_f32(&matrix, &rhs);
    let got = unsafe { matrix_simd::mul_vec4_sse_f32(&matrix, &rhs) };

    for i in 0..3 {
        assert!((got.data[i] - expected.data[i]).abs() < 1.0e-5);
    }
}

#[test]
fn mul_vec4_avx_f64_matches_scalar() {
    if !std::is_x86_feature_detected!("avx") {
        return;
    }

    let matrix = Matrix {
        data: [
            [1.0_f64, 2.0, 3.0, 4.0],
            [-1.5, 0.25, 10.0, -2.0],
            [0.0, 0.5, -0.5, 2.0],
        ],
    };
    let rhs = Vector {
        data: [2.0_f64, -1.0, 0.5, 4.0],
    };

    let expected = mul_vec4_scalar_f64(&matrix, &rhs);
    let got = unsafe { matrix_simd::mul_vec4_avx_f64(&matrix, &rhs) };

    for i in 0..3 {
        assert!((got.data[i] - expected.data[i]).abs() < 1.0e-12);
    }
}

#[test]
fn mul_vec4_avx_fma_f64_matches_scalar() {
    if !std::is_x86_feature_detected!("avx") || !std::is_x86_feature_detected!("fma") {
        return;
    }

    let matrix = Matrix {
        data: [
            [1.0_f64, 2.0, 3.0, 4.0],
            [-1.5, 0.25, 10.0, -2.0],
            [0.0, 0.5, -0.5, 2.0],
        ],
    };
    let rhs = Vector {
        data: [2.0_f64, -1.0, 0.5, 4.0],
    };

    let expected = mul_vec4_scalar_f64(&matrix, &rhs);
    let got = unsafe { matrix_simd::mul_vec4_avx_fma_f64(&matrix, &rhs) };

    for i in 0..3 {
        assert!((got.data[i] - expected.data[i]).abs() < 1.0e-12);
    }
}
