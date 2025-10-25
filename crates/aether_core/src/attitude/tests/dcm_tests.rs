#[cfg(test)]
mod tests {
    use super::*;
    use crate::attitude::tests::test_utils::*;
    use crate::attitude::{DirectionCosineMatrix, Euler, Quaternion};
    use crate::math::Matrix;
    use crate::matrix;
    use approx::assert_relative_eq;
    use num_traits::Float;

    const EPSILON: f64 = 1e-6;

    #[test]
    fn test_dcm_c1_quaternion_roundtrip() {
        let dcm1 = matrix![
            0.0, -1.0, 0.0;
            1.0, 0.0, 0.0;
            0.0, 0.0, 1.0
        ];
        // trigger c1
        test_round_trip(DirectionCosineMatrix::from(&dcm1), EPSILON);
    }

    #[test]
    fn test_dcm_c2_quaternion_roundtrip() {
        let dcm2: Matrix<f64, 3, 3> = matrix![
            0.0, 0.0,  1.0;
            0.0, -1.0,  0.0;
            1.0, 0.0, 0.0
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
        ]; // trigger c4
        test_round_trip(DirectionCosineMatrix::from(&dcm4), EPSILON);
    }

    #[test]
    fn test_dcm_euler_roundtrip() {}
}
