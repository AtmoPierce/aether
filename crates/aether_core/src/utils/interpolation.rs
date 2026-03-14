#[inline]
pub fn lerp_i128(x0: i128, y0: i128, x1: i128, y1: i128, x: i128) -> i128 {
    if x1 == x0 {
        return y1;
    }
    y0 + ((x - x0) * (y1 - y0)) / (x1 - x0)
}

#[inline]
pub fn lerp_f64(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    let dx = x1 - x0;
    if dx == 0.0 {
        return y1;
    }
    y0 + (x - x0) * (y1 - y0) / dx
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lerp_i128_midpoint() {
        assert_eq!(lerp_i128(0, 0, 10, 100, 5), 50);
    }

    #[test]
    fn lerp_f64_midpoint() {
        let y = lerp_f64(0.0, 0.0, 10.0, 100.0, 5.0);
        assert!((y - 50.0).abs() < 1e-12);
    }
}
