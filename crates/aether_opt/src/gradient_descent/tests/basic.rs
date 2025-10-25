#![cfg(feature = "std")]
use std::string::ToString;

#[cfg(test)]
mod gd_tests {
    extern crate std;
    use crate::gradient_descent::GradientDescentGeneric;
    use aether_core::math::Vector;

    // f(x) = 0.5 * ||x||^2, grad = x
    #[test]
    fn quad_fixed_step() {
        let mut gd: GradientDescentGeneric<f64, 2> = GradientDescentGeneric::new();
        gd.backtracking(false)
            .learning_rate(0.2)
            .momentum(0.0)
            .tolerances(1e-10, 1e-10)
            .max_iters(10_000);

        let x0 = Vector::new([10.0, -5.0]);
        let (x, fval, iters, conv) = gd.minimize(
            x0,
            |x| 0.5 * (x[0] * x[0] + x[1] * x[1]),
            |x| Vector::new([x[0], x[1]]),
        );

        assert!(conv, "did not converge in {iters}");
        assert!(fval < 1e-16);
        assert!(x[0].abs() < 1e-8 && x[1].abs() < 1e-8);
    }
}
