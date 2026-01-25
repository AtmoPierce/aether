mod vector {}

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
