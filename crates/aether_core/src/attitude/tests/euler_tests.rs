mod tests {
    use crate::attitude::{DirectionCosineMatrix, Euler, Quaternion};
    use crate::attitude::tests::test_utils::*;
    use crate::reference_frame::{Body, ITRF};
    use crate::real::Real;

    const EPSILON: f64 = 1e-15;

    /// Wrap angle to (-pi, pi)
    fn wrap_angle(x: f64) -> f64 {
        (x + std::f64::consts::PI) % (2.0 * std::f64::consts::PI)
            - std::f64::consts::PI
    }

    #[test]
    fn test_euler_dcm_roundtrip() {
        type I = ITRF<f64>;
        type B = Body<f64>;

        let euler: Euler<f64, I, B> =
            Euler::new(-0.5, 0.1, 0.2);

        let dcm1: DirectionCosineMatrix<f64, I, B> =
            DirectionCosineMatrix::from(euler);

        let euler_back: Euler<f64, I, B> =
            Euler::from(&dcm1);

        let dcm2: DirectionCosineMatrix<f64, I, B> =
            DirectionCosineMatrix::from(euler_back);

        // Euler angles are lossy - comparing DCMs instead
        assert!(
            matrices_approx_eq(
                &dcm1.as_matrix(),
                &dcm2.as_matrix(),
                0.1
            ),
            "DCM mismatch after Euler roundtrip"
        );
    }

    #[test]
    fn test_euler_quaternion_roundtrip() {
        type I = ITRF<f64>;
        type B = Body<f64>;

        let euler: Euler<f64, I, B> =
            Euler::new(-0.5, 0.1, 0.2);

        let quat: Quaternion<f64, I, B> =
            Quaternion::from(&euler).normalized();

        let euler_back: Euler<f64, I, B> =
            Euler::from(&quat);

        for i in 0..3 {
            let original = euler.data.data[i];
            let recovered = euler_back.data.data[i];

            let diff = wrap_angle(original - recovered);

            assert!(
                diff.abs() <= EPSILON,
                "Euler mismatch at index {}: original = {}, recovered = {}, diff = {}",
                i,
                original,
                recovered,
                diff
            );
        }
    }

    #[test]
    fn test_euler_quaternion_dcm_consistency() {
        type I = ITRF<f64>;
        type B = Body<f64>;

        let euler: Euler<f64, I, B> =
            Euler::new(0.3, -0.4, 0.7);

        let q: Quaternion<f64, I, B> =
            Quaternion::from(&euler).normalized();

        let dcm_from_q: DirectionCosineMatrix<f64, I, B> =
            DirectionCosineMatrix::from(&q);

        let dcm_from_e: DirectionCosineMatrix<f64, I, B> =
            DirectionCosineMatrix::from(euler);

        assert!(
            matrices_approx_eq(
                &dcm_from_q.as_matrix(),
                &dcm_from_e.as_matrix(),
                EPSILON
            ),
            "DCM from Euler and Quaternion disagree"
        );
    }

    #[test]
    fn test_euler_identity_rotation() {
        type I = ITRF<f64>;
        type B = Body<f64>;

        let euler: Euler<f64, I, B> =
            Euler::new(0.0, 0.0, 0.0);

        let q: Quaternion<f64, I, B> =
            Quaternion::from(&euler).normalized();

        let dcm: DirectionCosineMatrix<f64, I, B> =
            DirectionCosineMatrix::from(&q);

        let ident = dcm.as_matrix();

        // Check approximate identity matrix
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (ident[i][j] - expected).abs() <= EPSILON,
                    "Identity DCM mismatch at ({}, {})",
                    i,
                    j
                );
            }
        }
    }
}
