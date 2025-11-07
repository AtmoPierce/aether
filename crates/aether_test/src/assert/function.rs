/// Assert that the maximum absolute error between `f1(x)` and `f0(x)`
/// over a uniform grid on [a, b] is ≤ `tol`. -- This function allows you to test two functions against each other.
///
/// Returns the (x*, err*) of the worst point if you want to log it in the caller.
pub fn assert_max_abs_error_on_grid<F, G>(
    f0: F,
    f1: G,
    a: f64,
    b: f64,
    samples: usize,
    tol: f64,
) -> (f64, f64)
where
    F: Fn(f64) -> f64,
    G: Fn(f64) -> f64,
{
    assert!(samples > 0, "samples must be > 0");
    let mut worst_x = a;
    let mut max_e = 0.0_f64;

    for i in 0..=samples {
        let t = i as f64 / samples as f64;
        let x = a + (b - a) * t;
        let e = (f1(x) - f0(x)).abs();
        if e > max_e {
            max_e = e;
            worst_x = x;
        }
    }

    if max_e > tol {
        panic!(
            "max |f1 - f0| = {:.3e} at x = {:.15e} exceeds tol = {:.3e} on [{:.6},{:.6}] (n_samples={})",
            max_e, worst_x, tol, a, b, samples
        );
    }

    (worst_x, max_e)
}

