use aether_core::math::Vector;
use crate::assert::vector::*;

#[test]
fn test_vector_exact_equality() {
    let a = Vector::from([1.0, 2.0, 3.0]);
    let b = Vector::from([1.0, 2.0, 3.0]);
    assert_vector_approx_eq(&a, &b, 1e-12);
}

#[test]
fn test_vector_tolerance_pass() {
    let a = Vector::from([1.0, 2.0, 3.0]);
    let b = Vector::from([1.001, 1.999, 3.002]);
    assert_vector_approx_eq(&a, &b, 1e-2);
}

#[test]
#[should_panic]
fn test_vector_tolerance_fail() {
    let a = Vector::from([1.0, 2.0, 3.0]);
    let b = Vector::from([1.0, 2.0, 4.0]);
    assert_vector_approx_eq(&a, &b, 1e-4);
}

#[test]
fn test_vector_finite() {
    let v = Vector::from([1.0, 2.0, 3.0]);
    assert_vector_finite(&v);
}

#[test]
#[should_panic]
fn test_vector_not_finite() {
    let v = Vector::from([1.0, f64::NAN, 3.0]);
    assert_vector_finite(&v);
}
