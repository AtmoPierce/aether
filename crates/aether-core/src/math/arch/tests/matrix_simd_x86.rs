use crate::math::{Matrix, Vector};

#[test]
fn test_matrix_vector_mul_vec4_simd_f64() {
    let m = Matrix {
        data: [[1.0_f64, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0]],
    };
    let v = Vector {
        data: [1.0_f64, -1.0, 2.0, -2.0],
    };

    let r = m.mul_vec4_simd(&v);
    assert_eq!(r.data, [-3.0, -3.0]);
}

#[test]
fn test_matrix_vector_mul_vec4_simd_f32() {
    let m = Matrix {
        data: [[1.0_f32, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0]],
    };
    let v = Vector {
        data: [1.0_f32, -1.0, 2.0, -2.0],
    };

    let r = m.mul_vec4_simd(&v);
    assert_eq!(r.data, [-3.0, -3.0]);
}

#[test]
fn test_matrix_vector_mul_vec6_simd_f64() {
    let m = Matrix {
        data: [
            [1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0],
            [2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        ],
    };
    let v = Vector {
        data: [1.0_f64, -1.0, 2.0, -2.0, 3.0, -3.0],
    };

    let r = m.mul_vec6_simd(&v);
    assert_eq!(r.data, [-6.0, -6.0]);
}

#[test]
fn test_matrix_vector_mul_vec6_simd_f32() {
    let m = Matrix {
        data: [
            [1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0],
            [2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        ],
    };
    let v = Vector {
        data: [1.0_f32, -1.0, 2.0, -2.0, 3.0, -3.0],
    };

    let r = m.mul_vec6_simd(&v);
    assert_eq!(r.data, [-6.0, -6.0]);
}
