use crate::assert::coding::*;

#[cfg(feature = "bincode")]
#[test]
fn test_vector_equivalence() {
    use aether_core::math::Vector;

    let v = Vector::from([1.0, 2.0, 3.0]);
    assert_vector_equivalence(&v);
}

#[cfg(feature = "bincode")]
#[test]
fn test_quaternion_equivalence() {
    use aether_core::attitude::Quaternion;
    use aether_core::reference_frame::{ICRF,Body};

    let q: Quaternion<f64, ICRF<f64>,Body<f64>> = Quaternion::new(0.707, 0.707, 0.0, 0.0);
    assert_quaternion_equivalence(&q);

    let q: Quaternion<f64, ICRF<f64>,Body<f64>> = Quaternion::identity();
    assert_quaternion_equivalence(&q);
}