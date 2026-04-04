use core::ops::{Add, Mul, Sub};

/// Linearly interpolates between `a` and `b` using parameter `t`.
///
/// `t = 0.0` returns `a`, `t = 1.0` returns `b`, and values outside
/// that range extrapolate linearly.
#[inline]
pub fn lerp<T>(a: T, b: T, t: f64) -> T
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Mul<f64, Output = T>,
{
    a + (b - a) * t
}

/// Computes the normalized interpolation parameter of `value` in `[a, b]`.
#[inline]
pub fn inverse_lerp(a: f64, b: f64, value: f64) -> f64 {
    let denom = b - a;
    if denom.abs() <= f64::EPSILON {
        0.0
    } else {
        (value - a) / denom
    }
}

/// Remaps `value` from `[in_a, in_b]` into `[out_a, out_b]`.
#[inline]
pub fn remap(value: f64, in_a: f64, in_b: f64, out_a: f64, out_b: f64) -> f64 {
    lerp(out_a, out_b, inverse_lerp(in_a, in_b, value))
}

#[cfg(test)]
mod tests {
    use super::*;
    use aether_core::math::Vector;

    #[test]
    fn lerp_scalar_midpoint() {
        assert_eq!(lerp(0.0_f64, 10.0_f64, 0.25), 2.5);
    }

    #[test]
    fn lerp_vector_midpoint() {
        let a = Vector::new([0.0_f64, 0.0, 0.0]);
        let b = Vector::new([10.0_f64, -4.0, 2.0]);
        let v = lerp(a, b, 0.5);
        assert_eq!(v, Vector::new([5.0, -2.0, 1.0]));
    }

    #[test]
    fn inverse_lerp_degenerate_range() {
        assert_eq!(inverse_lerp(1.0, 1.0, 5.0), 0.0);
    }
}
