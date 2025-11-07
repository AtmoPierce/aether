use aether_core::math::Matrix;
use crate::assert::matrix::*;

#[test]
fn test_matrix_exact_equality() {
    let a = Matrix::new([[1.0, 2.0], [3.0, 4.0]]);
    let b = Matrix::new([[1.0, 2.0], [3.0, 4.0]]);
    assert_matrix_approx_eq(&a, &b, 1e-12);
}

#[test]
fn test_matrix_tolerance_pass() {
    let a = Matrix::new([[1.0, 2.0], [3.0, 4.0]]);
    let b = Matrix::new([[1.0009, 2.001], [2.999, 3.999]]);
    assert_matrix_approx_eq(&a, &b, 1e-2);
}

#[test]
#[should_panic]
fn test_matrix_tolerance_fail() {
    let a = Matrix::new([[1.0, 2.0], [3.0, 4.0]]);
    let b = Matrix::new([[1.0, 2.0], [3.0, 5.0]]);
    assert_matrix_approx_eq(&a, &b, 1e-4);
}

#[test]
fn test_matrix_finite() {
    let a = Matrix::new([[1.0, 2.0], [3.0, 4.0]]);
    assert_matrix_finite(&a);
}

#[test]
#[should_panic]
fn test_matrix_not_finite() {
    let a = Matrix::new([[1.0, f64::INFINITY], [3.0, 4.0]]);
    assert_matrix_finite(&a);
}

#[test]
fn test_matrix_norms() {
    let a = Matrix::new([[1.0, 2.0], [3.0, 4.0]]);
    let b = Matrix::new([[1.1, 1.9], [2.9, 4.2]]);
    let l_inf = matrix_max_abs_diff(&a, &b);
    let frob = matrix_diff_frobenius(&a, &b);
    assert!(l_inf > 0.0 && frob > 0.0);
}

#[test]
fn test_matrix_symmetric() {
    let a = Matrix::new([[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]]);
    assert_matrix_symmetric(&a, 1e-12);
}

#[test]
#[should_panic]
fn test_matrix_not_symmetric() {
    let a = Matrix::new([[1.0, 2.0], [3.01, 4.0]]);
    assert_matrix_symmetric(&a, 1e-3);
}
