#[cfg(test)]
mod tests {
    use super::*;
    use crate::attitude::tests::test_utils::*;
    use crate::attitude::{DirectionCosineMatrix, Euler, Quaternion};
    use crate::math::Matrix;
    use crate::matrix;
    use crate::reference_frame::Body;
    use approx::assert_relative_eq;
    use crate::real::Real;

    const EPSILON: f64 = 1e-6;

    #[test]
    fn test_dcm_c1_quaternion_roundtrip() {
        let dcm1 = matrix![
            0.0, -1.0, 0.0;
            1.0,  0.0, 0.0;
            0.0,  0.0, 1.0
        ];
        // trigger c1
        test_round_trip(DirectionCosineMatrix::from(&dcm1), EPSILON);
    }

    #[test]
    fn test_dcm_c2_quaternion_roundtrip() {
        let dcm2: Matrix<f64, 3, 3> = matrix![
            0.0,  0.0,  1.0;
            0.0, -1.0,  0.0;
            1.0,  0.0,  0.0
        ];
        // trigger c2
        test_round_trip(DirectionCosineMatrix::from(&dcm2), EPSILON);
    }

    #[test]
    fn test_dcm_c3_quaternion_roundtrip() {
        let dcm3 = matrix![
            -1.0,  0.0,  0.0;
             0.0,  1.0,  0.0;
             0.0,  0.0, -1.0
        ];
        // trigger c3
        test_round_trip(DirectionCosineMatrix::from(&dcm3), EPSILON);
    }

    #[test]
    fn test_dcm_c4_quaternion_roundtrip() {
        let dcm4 = matrix![
            -1.0,  0.0,  0.0;
             0.0, -1.0,  0.0;
             0.0,  0.0,  1.0
        ];
        // trigger c4
        test_round_trip(DirectionCosineMatrix::from(&dcm4), EPSILON);
    }

    #[test]
    fn test_dcm_try_from_quaternion_roundtrip() {
        // Generic, non-degenerate DCM
        let dcm = matrix![
             0.36, -0.48,  0.80;
             0.80,  0.60,  0.00;
            -0.48,  0.64,  0.60
        ];

        let dcm1:DirectionCosineMatrix<f64, Body<f64>, Body<f64>>  = DirectionCosineMatrix::from(&dcm);

        // Exercise your exact TryFrom implementation
        let q = Quaternion::try_from(&dcm1)
            .expect("TryFrom<DCM> failed");

        let dcm2 = DirectionCosineMatrix::from(&q);

        let m1 = dcm1.as_matrix();
        let m2 = dcm2.as_matrix();

        for r in 0..3 {
            for c in 0..3 {
                assert!(
                    (m1[r][c] - m2[r][c]).abs() <= EPSILON,
                    "Mismatch at ({}, {}): {} vs {}",
                    r, c, m1[r][c], m2[r][c]
                );
            }
        }
    }

    #[test]
    fn test_dcm_euler_roundtrip() {
        let euler: Euler<f64, Body<f64>, Body<f64>> = Euler::new(-0.7, 0.3, 1.2);

        let dcm1 = DirectionCosineMatrix::from(euler);
        let euler_back = Euler::from(&dcm1);
        let dcm2 = DirectionCosineMatrix::from(euler_back);

        let m1 = dcm1.as_matrix();
        let m2 = dcm2.as_matrix();

        for r in 0..3 {
            for c in 0..3 {
                assert!(
                    (m1[r][c] - m2[r][c]).abs() <= 1e-1,
                    "Mismatch at ({}, {}): {} vs {}",
                    r, c, m1[r][c], m2[r][c]
                );
            }
        }
    }
}
