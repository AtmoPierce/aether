#![cfg(feature = "std")]
use std::string::ToString;

// ------------------ tests (need std) ------------------
#[cfg(all(test, feature = "std"))]
mod gd_stochastic_tests {
    use crate::gradient_descent::GDStochastic;
    use aether_core::math::Vector;
    use aether_rand::randomizers::xor_shift::XorShift64Star;

    // Noisy quadratic: f(x) = 0.5 ||x||^2, grad = x + noise``
    #[test]
    fn sgd_noisy_quadratic_converges() {
        let mut sgd: GDStochastic<2> = GDStochastic::new(0xDEADBEEF);
        sgd.learning_rate(0.05)
            .decay(5e-4)
            .momentum(0.6)
            .batch_size(8)
            .ema(0.95)
            .tolerances(2e-3, 5e-9) // tol_grad_ema, tol_step
            .tol_f(1e-6) // also stop on tiny f
            .min_iters(500) // avoid premature stop
            .max_iters(200_000);

        let x0 = Vector::new([5.0, -3.0]);
        let (x, fval, iters, conv) = sgd.minimize(
            x0,
            |x| 0.5 * (x[0] * x[0] + x[1] * x[1]),
            |x, rng| {
                // unbiased noisy grad sample
                Vector::new([
                    x[0] + (rng.next_f64() - 0.5) * 0.2,
                    x[1] + (rng.next_f64() - 0.5) * 0.2,
                ])
            },
        );

        println!("SGD result: x={x:?}, fval={fval}, conv={conv}");
        assert!(conv, "did not meet tolerances");
        assert!(fval < 1e-3, "fval too high: {fval}");
        assert!(x[0].abs() < 0.05 && x[1].abs() < 0.05, "x too far: {x:?}");
    }

    // Same problem with bounds and pure SGD (batch_size=1)
    #[test]
    fn sgd_with_bounds() {
        let mut sgd: GDStochastic<3> = GDStochastic::new(123);
        sgd.learning_rate(0.05)
            .decay(5e-4)
            .momentum(0.6)
            .batch_size(8)
            .ema(0.95)
            .tolerances(2e-3, 5e-9) // tol_grad_ema, tol_step
            .tol_f(1e-6) // also stop on tiny f
            .min_iters(500) // avoid premature stop
            .max_iters(200_000);

        let x0 = Vector::new([0.9, 0.4, 0.09]);
        let (x, _fval, _it, conv) = sgd.minimize(
            x0,
            |x| 0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]),
            |x, rng| {
                // grad = x + noise (σ=0.05)
                Vector::new([
                    x[0] + rng.normal01() * 0.05,
                    x[1] + rng.normal01() * 0.05,
                    x[2] + rng.normal01() * 0.05,
                ])
            },
        );

        assert!(conv);
        // Should be near 0 and within bounds
        for i in 0..3 {
            assert!(x[i].abs() <= 0.2, "axis {i} too far: {}", x[i]);
        }
    }
}
