#[cfg(test)]
mod tests {
    use crate::design::functions::sqrt_f;
    use crate::functions::tests::ulp_diff_f64;

    #[test]
    fn compare_sqrt_f_to_core_on_design_interval() {
        let a = 0.0_f64;
        let b = 100.0_f64;
        let samples = 100_000usize;

        let mut max_abs = 0.0_f64;
        let mut max_rel = 0.0_f64;
        let mut max_ulp = 0u64;
        let mut x = a;
        let mut worst_x = 0.0;
        let mut worst_ours = 0.0;
        let mut worst_core = 0.0;

        for i in 0..=samples {
            let t = i as f64 / samples as f64;
            let x = a + (b - a) * t;

            let ours = sqrt_f(x);
            let core = x.sqrt();
            let abs = (ours - core).abs();
            let ulp = ulp_diff_f64(ours, core);
            if ulp > max_ulp {
                max_ulp = ulp;
                max_abs = abs;
                worst_x = x;
                worst_ours = ours;
                worst_core = core;
            }
        }
        assert!(max_ulp <= 1, "too many ULPs off: {}", max_ulp);
        assert!(max_rel < 1e-15, "relative error too large: {}", max_rel);
    }
}
