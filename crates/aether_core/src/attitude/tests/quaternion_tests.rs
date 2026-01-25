mod tests {
    use crate::attitude::{DirectionCosineMatrix, Euler, Quaternion};
    use crate::math::Vector;
    use crate::reference_frame::{Body, ITRF, ReferenceFrame};
    use crate::real::Real;

    const EPSILON: f64 = 1e-15;

    /// q and -q represent the same rotation
    fn quat_equivalent<T, A, B>(
        a: &Quaternion<T, A, B>,
        b: &Quaternion<T, A, B>,
    )
    where
        T: Real,
        A: ReferenceFrame,
        B: ReferenceFrame,
    {
        let dot = a.dot(b).abs();
        assert!(
            (dot - T::ONE).abs() <= T::from_f64(EPSILON),
            "Quaternions not equivalent: |dot| = {:?}",
            dot
        );
    }

    #[test]
    fn quaternion_euler_roundtrip_itrf_to_body() {
        type I = ITRF<f64>;
        type B = Body<f64>;

        let q_ib: Quaternion<f64, I, B> =
            Quaternion::new(0.9, 0.1, 0.2, 0.3).normalized();

        let euler: Euler<f64, I, B> = Euler::from(&q_ib);

        let q_back: Quaternion<f64, I, B> =
            Quaternion::from(&euler).normalized();

        quat_equivalent(&q_ib, &q_back);
    }

    #[test]
    fn quaternion_dcm_roundtrip_itrf_to_body() {
        type I = ITRF<f64>;
        type B = Body<f64>;

        let q_ib: Quaternion<f64, I, B> =
            Quaternion::new(0.9, 0.1, 0.2, 0.3).normalized();

        let dcm_ib: DirectionCosineMatrix<f64, I, B> =
            DirectionCosineMatrix::from(&q_ib);

        let q_back: Quaternion<f64, I, B> =
            Quaternion::try_from(&dcm_ib)
                .unwrap()
                .normalized();

        quat_equivalent(&q_ib, &q_back);
    }

    #[test]
    fn quaternion_sign_invariance_rotation() {
        type I = ITRF<f64>;
        type B = Body<f64>;

        let q_ib: Quaternion<f64, I, B> =
            Quaternion::new(0.9, 0.1, 0.2, 0.3).normalized();

        let v_i = Vector::new([1.0, 2.0, 3.0]);

        let v_b1 = q_ib.rotate_vector(v_i);
        let v_b2 = (-q_ib).rotate_vector(v_i);

        assert!(
            (v_b1 - v_b2).norm() <= EPSILON,
            "Rotation not invariant under sign flip"
        );
    }

    #[test]
    fn quaternion_vs_dcm_vector_rotation_itrf_to_body() {
        use crate::coordinate::Cartesian;
        
        type I = ITRF<f64>;
        type B = Body<f64>;

        let q_ib: Quaternion<f64, I, B> = Quaternion::new(0.7, -0.2, 0.3, 0.6).normalized();
        let dcm_ib: DirectionCosineMatrix<f64, I, B> = DirectionCosineMatrix::from(&q_ib);

        let v_i = Cartesian::<f64, I>::new(4.0, -1.0, 2.0);

        let v_b_q = q_ib * v_i;
        let v_b_dcm = dcm_ib * v_i;

        assert!(
            (v_b_q - v_b_dcm).norm() <= EPSILON,
            "Quaternion and DCM rotations disagree"
        );
    }

    #[test]
    fn dcm_euler_quaternion_full_roundtrip_itrf_to_body() {
        type I = ITRF<f64>;
        type B = Body<f64>;

        let q0: Quaternion<f64, I, B> =
            Quaternion::new(0.4, 0.5, -0.1, 0.75).normalized();

        let dcm: DirectionCosineMatrix<f64, I, B> =
            DirectionCosineMatrix::from(&q0);

        let euler: Euler<f64, I, B> = Euler::from(&dcm);

        let q1: Quaternion<f64, I, B> =
            Quaternion::from(&euler).normalized();

        quat_equivalent(&q0, &q1);
    }
}
