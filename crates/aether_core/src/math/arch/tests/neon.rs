use crate::math::{
    arch::arm::neon,
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

fn mul_vec6_scalar_f32<const M: usize>(
    matrix: &Matrix<f32, M, 6>,
    rhs: &Vector<f32, 6>,
) -> Vector<f32, M> {
    let mut out = Vector { data: [0.0; M] };
    for i in 0..M {
        out.data[i] = matrix.data[i][0] * rhs.data[0]
            + matrix.data[i][1] * rhs.data[1]
            + matrix.data[i][2] * rhs.data[2]
            + matrix.data[i][3] * rhs.data[3]
            + matrix.data[i][4] * rhs.data[4]
            + matrix.data[i][5] * rhs.data[5];
    }
    out
}

fn mul_vec6_scalar_f64<const M: usize>(
    matrix: &Matrix<f64, M, 6>,
    rhs: &Vector<f64, 6>,
) -> Vector<f64, M> {
    let mut out = Vector { data: [0.0; M] };
    for i in 0..M {
        out.data[i] = matrix.data[i][0] * rhs.data[0]
            + matrix.data[i][1] * rhs.data[1]
            + matrix.data[i][2] * rhs.data[2]
            + matrix.data[i][3] * rhs.data[3]
            + matrix.data[i][4] * rhs.data[4]
            + matrix.data[i][5] * rhs.data[5];
    }
    out
}

fn mul_mat_scalar_f32<const M: usize, const N: usize, const P: usize>(
    matrix_a: &Matrix<f32, M, N>,
    matrix_b: &Matrix<f32, N, P>,
) -> Matrix<f32, M, P> {
    let mut out = Matrix { data: [[0.0; P]; M] };
    for i in 0..M {
        for j in 0..P {
            let mut acc = 0.0;
            for k in 0..N {
                acc += matrix_a.data[i][k] * matrix_b.data[k][j];
            }
            out.data[i][j] = acc;
        }
    }
    out
}

fn mul_mat_scalar_f64<const M: usize, const N: usize, const P: usize>(
    matrix_a: &Matrix<f64, M, N>,
    matrix_b: &Matrix<f64, N, P>,
) -> Matrix<f64, M, P> {
    let mut out = Matrix { data: [[0.0; P]; M] };
    for i in 0..M {
        for j in 0..P {
            let mut acc = 0.0;
            for k in 0..N {
                acc += matrix_a.data[i][k] * matrix_b.data[k][j];
            }
            out.data[i][j] = acc;
        }
    }
    out
}

#[test]
fn dot4_neon_f32_matches_scalar() {
    let a = Vector {
        data: [1.0_f32, -2.0, 3.5, 0.25],
    };
    let b = Vector {
        data: [2.0_f32, 1.5, -4.0, 10.0],
    };

    let expected = dot4_f32_scalar(&a, &b);
    let got = unsafe { neon::dot4_neon_f32(&a, &b) };

    assert!((got - expected).abs() < 1.0e-6);
}

#[test]
fn dot4_neon_f64_matches_scalar() {
    let a = Vector {
        data: [1.0_f64, -2.0, 3.5, 0.25],
    };
    let b = Vector {
        data: [2.0_f64, 1.5, -4.0, 10.0],
    };

    let expected = dot4_f64_scalar(&a, &b);
    let got = unsafe { neon::dot4_neon_f64(&a, &b) };

    assert!((got - expected).abs() < 1.0e-12);
}

#[test]
fn mul_vec4_neon_f32_matches_scalar() {
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
    let got = unsafe { neon::mul_vec4_neon_f32(&matrix, &rhs) };

    for i in 0..3 {
        assert!((got.data[i] - expected.data[i]).abs() < 1.0e-5);
    }
}

#[test]
fn mul_vec4_neon_f64_matches_scalar() {
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
    let got = unsafe { neon::mul_vec4_neon_f64(&matrix, &rhs) };

    for i in 0..3 {
        assert!((got.data[i] - expected.data[i]).abs() < 1.0e-12);
    }
}

#[test]
fn mul_vec6_neon_f32_matches_scalar() {
    let matrix = Matrix {
        data: [
            [1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0],
            [-1.5, 0.25, 10.0, -2.0, 3.0, -1.0],
            [0.0, 0.5, -0.5, 2.0, -3.0, 1.0],
        ],
    };
    let rhs = Vector {
        data: [2.0_f32, -1.0, 0.5, 4.0, -0.5, 1.5],
    };

    let expected = mul_vec6_scalar_f32(&matrix, &rhs);
    let got = unsafe { neon::mul_vec6_neon_f32(&matrix, &rhs) };

    for i in 0..3 {
        assert!((got.data[i] - expected.data[i]).abs() < 1.0e-5);
    }
}

#[test]
fn mul_vec6_neon_f64_matches_scalar() {
    let matrix = Matrix {
        data: [
            [1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0],
            [-1.5, 0.25, 10.0, -2.0, 3.0, -1.0],
            [0.0, 0.5, -0.5, 2.0, -3.0, 1.0],
        ],
    };
    let rhs = Vector {
        data: [2.0_f64, -1.0, 0.5, 4.0, -0.5, 1.5],
    };

    let expected = mul_vec6_scalar_f64(&matrix, &rhs);
    let got = unsafe { neon::mul_vec6_neon_f64(&matrix, &rhs) };

    for i in 0..3 {
        assert!((got.data[i] - expected.data[i]).abs() < 1.0e-12);
    }
}

#[test]
fn mul_mat3_neon_f32_matches_scalar() {
    let a = Matrix {
        data: [[1.0_f32, 2.0, 3.0], [0.5, -1.0, 2.0], [3.0, 0.0, -2.0]],
    };
    let b = Matrix {
        data: [[2.0_f32, 0.0, 1.0], [1.0, 3.0, -2.0], [0.0, 1.0, 4.0]],
    };

    let expected = mul_mat_scalar_f32(&a, &b);
    let got = unsafe { neon::mul_mat3_neon_f32(&a, &b) };

    for i in 0..3 {
        for j in 0..3 {
            assert!((got.data[i][j] - expected.data[i][j]).abs() < 1.0e-5);
        }
    }
}

#[test]
fn mul_mat4_neon_f64_matches_scalar() {
    let a = Matrix {
        data: [
            [1.0_f64, 2.0, 3.0, 4.0],
            [0.5, -1.0, 2.0, 1.5],
            [3.0, 0.0, -2.0, 1.0],
            [2.0, 1.0, 0.0, -1.0],
        ],
    };
    let b = Matrix {
        data: [
            [2.0_f64, 0.0, 1.0, -1.0],
            [1.0, 3.0, -2.0, 0.5],
            [0.0, 1.0, 4.0, 2.0],
            [1.5, -1.0, 0.5, 2.5],
        ],
    };

    let expected = mul_mat_scalar_f64(&a, &b);
    let got = unsafe { neon::mul_mat4_neon_f64(&a, &b) };

    for i in 0..4 {
        for j in 0..4 {
            assert!((got.data[i][j] - expected.data[i][j]).abs() < 1.0e-12);
        }
    }
}

#[test]
fn mul_mat6_neon_f64_matches_scalar() {
    let a = Matrix {
        data: [
            [1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0],
            [2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            [3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
            [4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
            [5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            [6.0, 7.0, 8.0, 9.0, 10.0, 11.0],
        ],
    };
    let b = Matrix {
        data: [
            [11.0_f64, 10.0, 9.0, 8.0, 7.0, 6.0],
            [10.0, 9.0, 8.0, 7.0, 6.0, 5.0],
            [9.0, 8.0, 7.0, 6.0, 5.0, 4.0],
            [8.0, 7.0, 6.0, 5.0, 4.0, 3.0],
            [7.0, 6.0, 5.0, 4.0, 3.0, 2.0],
            [6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
        ],
    };

    let expected = mul_mat_scalar_f64(&a, &b);
    let got = unsafe { neon::mul_mat6_neon_f64(&a, &b) };

    for i in 0..6 {
        for j in 0..6 {
            assert!((got.data[i][j] - expected.data[i][j]).abs() < 1.0e-12);
        }
    }
}

#[test]
fn mul_matrix_neon_f32_larger_matches_scalar() {
    let a = Matrix {
        data: [
            [1.0_f32, 2.0, 3.0, 4.0, 5.0],
            [2.0, 3.0, 4.0, 5.0, 6.0],
            [3.0, 4.0, 5.0, 6.0, 7.0],
            [4.0, 5.0, 6.0, 7.0, 8.0],
            [5.0, 6.0, 7.0, 8.0, 9.0],
            [6.0, 7.0, 8.0, 9.0, 10.0],
            [7.0, 8.0, 9.0, 10.0, 11.0],
        ],
    };
    let b = Matrix {
        data: [
            [1.0_f32, 0.0, 2.0, 1.0, -1.0, 3.0],
            [0.0, 1.0, 1.0, 2.0, 3.0, -2.0],
            [2.0, 1.0, 0.0, -1.0, 2.0, 1.0],
            [1.0, 2.0, 3.0, 0.0, -2.0, 2.0],
            [3.0, -1.0, 2.0, 1.0, 0.0, 1.0],
        ],
    };

    let expected = mul_mat_scalar_f32(&a, &b);
    let got = unsafe { neon::mul_matrix_neon_f32::<7, 5, 6>(&a, &b) };

    for i in 0..7 {
        for j in 0..6 {
            assert!((got.data[i][j] - expected.data[i][j]).abs() < 1.0e-5);
        }
    }
}
