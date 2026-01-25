#[cfg(all(test, feature = "std"))]
mod stochastic_parallel_tests {
    use crate::gradient_descent::GDStochasticParallel;
    use aether::math::Vector;

    // Converges on f(x) = 0.5||x||^2 with noisy grads.
    #[test]
    fn parallel_sgd_noisy_quadratic_converges() {
        let mut sgd: GDStochasticParallel<f64, 2> = GDStochasticParallel::new(0xFEED_FACEu64);
        sgd.learning_rate(0.05)
            .decay(5e-4)
            .momentum(0.6)
            .batch_size(256)
            .ema(0.95)
            .tolerances(1e-3, 5e-9)
            .tol_f(1e-8)
            .min_iters(200)
            .max_iters(50_000)
            .threads(2); // 0/None => auto

        let x0 = Vector::new([4.0, -2.0]);
        let (x, fval, _iters, conv) = sgd.minimize(
            x0,
            |x| 0.5 * (x[0] * x[0] + x[1] * x[1]),
            |x, rng| {
                // unbiased noisy grad: x + noise in [-0.1, 0.1]
                let n0 = (rng.next_f64() - 0.5) * 0.2;
                let n1 = (rng.next_f64() - 0.5) * 0.2;
                Vector::new([x[0] + n0, x[1] + n1])
            },
        );

        assert!(conv, "did not converge: f={fval}, x={x:?}");
        assert!(fval < 1e-6, "objective not small enough: {fval}");
        assert!(
            x[0].abs() < 1e-3 && x[1].abs() < 1e-3,
            "solution too far: {x:?}"
        );
    }
}
