use crate::math::Matrix;

fn assert_matrix_close_f32<const M: usize, const N: usize>(
    lhs: &Matrix<f32, M, N>,
    rhs: &Matrix<f32, M, N>,
) {
    for i in 0..M {
        for j in 0..N {
            assert!((lhs.data[i][j] - rhs.data[i][j]).abs() < 1.0e-4);
        }
    }
}

fn assert_matrix_close_f64<const M: usize, const N: usize>(
    lhs: &Matrix<f64, M, N>,
    rhs: &Matrix<f64, M, N>,
) {
    for i in 0..M {
        for j in 0..N {
            assert!((lhs.data[i][j] - rhs.data[i][j]).abs() < 1.0e-10);
        }
    }
}

#[test]
fn test_matrix_mul_4x4_simd_wrapper_f64() {
    let a = Matrix {
        data: [
            [1.0_f64, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ],
    };
    let b = Matrix {
        data: [
            [16.0_f64, 15.0, 14.0, 13.0],
            [12.0, 11.0, 10.0, 9.0],
            [8.0, 7.0, 6.0, 5.0],
            [4.0, 3.0, 2.0, 1.0],
        ],
    };

    let simd = a.mul_mat4_simd(&b);
    let scalar = a * b;
    assert_matrix_close_f64(&simd, &scalar);
}

#[test]
fn test_matrix_mul_4x4_simd_wrapper_f32() {
    let a = Matrix {
        data: [
            [1.0_f32, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ],
    };
    let b = Matrix {
        data: [
            [16.0_f32, 15.0, 14.0, 13.0],
            [12.0, 11.0, 10.0, 9.0],
            [8.0, 7.0, 6.0, 5.0],
            [4.0, 3.0, 2.0, 1.0],
        ],
    };

    let simd = a.mul_mat4_simd(&b);
    let scalar = a * b;
    assert_matrix_close_f32(&simd, &scalar);
}

#[test]
fn test_matrix_mul_3x3_simd_wrapper_f64() {
    let a = Matrix {
        data: [[1.0_f64, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
    };
    let b = Matrix {
        data: [[9.0_f64, 8.0, 7.0], [6.0, 5.0, 4.0], [3.0, 2.0, 1.0]],
    };

    let simd = a.mul_mat3_simd(&b);
    let scalar = a * b;
    assert_matrix_close_f64(&simd, &scalar);
}

#[test]
fn test_matrix_mul_3x3_simd_wrapper_f32() {
    let a = Matrix {
        data: [[1.0_f32, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
    };
    let b = Matrix {
        data: [[9.0_f32, 8.0, 7.0], [6.0, 5.0, 4.0], [3.0, 2.0, 1.0]],
    };

    let simd = a.mul_mat3_simd(&b);
    let scalar = a * b;
    assert_matrix_close_f32(&simd, &scalar);
}

#[test]
fn test_matrix_mul_6x6_simd_wrapper_f64() {
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

    let simd = a.mul_mat6_simd(&b);
    let scalar = a * b;
    assert_matrix_close_f64(&simd, &scalar);
}

#[test]
fn test_matrix_mul_6x6_simd_wrapper_f32() {
    let a = Matrix {
        data: [
            [1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0],
            [2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            [3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
            [4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
            [5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            [6.0, 7.0, 8.0, 9.0, 10.0, 11.0],
        ],
    };
    let b = Matrix {
        data: [
            [11.0_f32, 10.0, 9.0, 8.0, 7.0, 6.0],
            [10.0, 9.0, 8.0, 7.0, 6.0, 5.0],
            [9.0, 8.0, 7.0, 6.0, 5.0, 4.0],
            [8.0, 7.0, 6.0, 5.0, 4.0, 3.0],
            [7.0, 6.0, 5.0, 4.0, 3.0, 2.0],
            [6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
        ],
    };

    let simd = a.mul_mat6_simd(&b);
    let scalar = a * b;
    assert_matrix_close_f32(&simd, &scalar);
}

#[test]
fn test_matrix_mul_simd_generic_larger_f64() {
    let a = Matrix {
        data: [
            [1.0_f64, 2.0, 3.0, 4.0, 5.0],
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
            [1.0_f64, 0.0, 2.0, 1.0, -1.0, 3.0],
            [0.0, 1.0, 1.0, 2.0, 3.0, -2.0],
            [2.0, 1.0, 0.0, -1.0, 2.0, 1.0],
            [1.0, 2.0, 3.0, 0.0, -2.0, 2.0],
            [3.0, -1.0, 2.0, 1.0, 0.0, 1.0],
        ],
    };

    let simd = a.mul_matrix_simd(&b);
    let scalar = a * b;
    assert_matrix_close_f64(&simd, &scalar);
}

#[test]
fn test_matrix_mul_simd_generic_larger_f32() {
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

    let simd = a.mul_matrix_simd(&b);
    let scalar = a * b;
    assert_matrix_close_f32(&simd, &scalar);
}
