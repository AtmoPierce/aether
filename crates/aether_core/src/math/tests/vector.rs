mod vector {}

#[cfg(test)]
mod arithmetic_tests {
    use crate::math::vector::Vector;

    #[test]
    fn borrowed_vector_dot_mul_works() {
        let a = Vector::new([1.0_f64, 2.0, 3.0]);
        let b = Vector::new([4.0_f64, 5.0, 6.0]);

        let dot_ref = &a * &b;
        let dot_owned = a * b;

        assert_eq!(dot_ref, 32.0);
        assert_eq!(dot_owned, 32.0);
    }

    #[test]
    fn vector_dot_and_norm_n4_fast_path() {
        let a = Vector::new([1.0_f64, 2.0, 3.0, 4.0]);
        let b = Vector::new([5.0_f64, 6.0, 7.0, 8.0]);

        let dot = a.dot(&b);
        assert_eq!(dot, 70.0);

        let norm = a.norm();
        assert!((norm - (30.0_f64).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn vector_dot_n4_fast_path_f32() {
        let a = Vector::new([1.0_f32, 2.0, 3.0, 4.0]);
        let b = Vector::new([5.0_f32, 6.0, 7.0, 8.0]);

        let dot = a.dot(&b);
        assert_eq!(dot, 70.0_f32);
    }

    #[test]
    fn vector_dot4_simd_f64() {
        let a = Vector::new([1.0_f64, 2.0, 3.0, 4.0]);
        let b = Vector::new([5.0_f64, 6.0, 7.0, 8.0]);

        let dot = a.dot4_simd(&b);
        assert_eq!(dot, 70.0_f64);
    }

    #[test]
    fn vector_dot4_simd_f32() {
        let a = Vector::new([1.0_f32, 2.0, 3.0, 4.0]);
        let b = Vector::new([5.0_f32, 6.0, 7.0, 8.0]);

        let dot = a.dot4_simd(&b);
        assert_eq!(dot, 70.0_f32);
    }

    #[test]
    fn borrowed_add_assign_sub_assign_work() {
        let mut a = Vector::new([1.0_f64, 2.0, 3.0]);
        let b = Vector::new([4.0_f64, 5.0, 6.0]);

        a += &b;
        assert_eq!(a.data, [5.0, 7.0, 9.0]);

        a -= &b;
        assert_eq!(a.data, [1.0, 2.0, 3.0]);
    }
}

// In some test or module compiled with std:
#[cfg(feature = "std")]
mod thread_safety_checks {
    use crate::math::vector::Vector;

    // Compile-time checks
    fn assert_send_sync<T: Send + Sync>() {}
    #[test]
    fn vector_f64_is_send_sync() {
        assert_send_sync::<Vector<f64, 3>>();
    }

    // Runtime smoke test: move a Vector across threads and back
    #[test]
    fn can_move_vector_across_threads() {
        let v = Vector::new([1.0_f64, 2.0, 3.0]);
        let handle = std::thread::spawn(move || {
            // Owned 'v' moved into thread
            let mut w = v;
            w[0] += 1.0;
            w
        });
        let out = handle.join().unwrap();
        assert_eq!(out.data, [2.0, 2.0, 3.0]);
    }
}
