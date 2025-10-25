use super::*;
use crate::attitude::{DirectionCosineMatrix, Euler, Quaternion};
use crate::reference_frame::Body;
use crate::{math::Matrix, matrix};
use approx::assert_relative_eq;
use num_traits::Float;

pub fn matrices_approx_eq<T: Float + std::fmt::Debug>(
    a: &Matrix<T, 3, 3>,
    b: &Matrix<T, 3, 3>,
    epsilon: T,
) -> bool {
    let mut ok = true;
    for r in 0..3 {
        for c in 0..3 {
            let diff = (a.data[r][c] - b.data[r][c]).abs();
            if diff > epsilon {
                println!(
                    "Mismatch at ({}, {}): a = {:?}, b = {:?}, diff = {:?}",
                    r, c, a.data[r][c], b.data[r][c], diff
                );
                ok = false;
            }
        }
    }
    ok
}

pub fn test_round_trip<T>(dcm: DirectionCosineMatrix<T, Body<f64>, Body<f64>>, epsilon: T)
where
    T: Float + std::fmt::Debug,
{
    let q: Quaternion<T> = (&dcm).try_into().unwrap();
    let q_n = q.normalized();
    let dcm_rt: DirectionCosineMatrix<T, Body<f64>, Body<f64>> = DirectionCosineMatrix::from(q_n);
    assert!(matrices_approx_eq(
        dcm.as_matrix(),
        dcm_rt.as_matrix(),
        epsilon,
    ));
}
