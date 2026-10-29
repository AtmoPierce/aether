mod quantities {
    use crate::ocp::fp6::Fp6E3M2;

    fn approx_eq(a: f32, b: f32, eps: f32) -> bool {
        let approx = (a - b).abs();
        println!("Approximation: {}", approx);
        approx <= eps
    }

    // -------------------------
    // FP6 E3M2 tests
    // -------------------------

    #[test]
    fn fp6e3m2_zero_and_sign() {
        // +0: s=0, e=0, m=0 => 0b0_000_00 = 0x00
        // -0: s=1, e=0, m=0 => 0b1_000_00 = 0x20
        let pz = Fp6E3M2::from_bits(0x00);
        let nz = Fp6E3M2::from_bits(0x20);

        let pzf: f32 = pz.into();
        let nzf: f32 = nz.into();

        println!("Fp6 +0 = {}, f32 = {}", pz, pzf);
        println!("Fp6 -0 = {}, f32 = {}", nz, nzf);

        assert_eq!(pzf, 0.0);
        assert_eq!(nzf, -0.0);
        assert_eq!(pzf.to_bits(), 0.0f32.to_bits());
        assert_eq!(nzf.to_bits(), (-0.0f32).to_bits());
    }

    #[test]
    fn fp6e3m2_subnormal_and_min_normal_decode() {
        let sub = Fp6E3M2::MIN_SUBNORMAL_POS; // 0b0_000_01
        let nor = Fp6E3M2::MIN_NORMAL_POS; // 0b0_001_00

        let sub_f: f32 = sub.into();
        let nor_f: f32 = nor.into();

        println!("Fp6 MIN_SUBNORMAL_POS = {}, f32 = {:.8e}", sub, sub_f);
        println!("Fp6 MIN_NORMAL_POS    = {}, f32 = {:.8e}", nor, nor_f);

        // For E3M2, bias=3:
        // subnormal: e=0, m=1..3:
        //   value = (m/4) * 2^(1-bias) = (m/4) * 2^-2
        // m=1 => 1/4 * 1/4 = 1/16
        assert!(approx_eq(sub_f, 1.0 / 16.0, 1e-7));

        // min normal: e=1, m=0:
        //   value = (1 + 0/4) * 2^(1-3) = 1 * 2^-2 = 1/4
        assert!(approx_eq(nor_f, 0.25, 1e-7));
    }

    #[test]
    fn fp6e3m2_max_finite_decode() {
        let max = Fp6E3M2::MAX_FINITE;
        let v: f32 = max.into();

        println!("Fp6 MAX_FINITE => {} (f32 = {:.8e})", max, v);

        // Expected: (1 + 3/4) * 2^(6-3) = 1.75 * 8 = 14
        assert!(approx_eq(v, 28.0, 1e-6));
    }

    #[test]
    fn fp6e3m2_roundtrip_simple_values() {
        let vals = [0.0f32, 1.0, -1.0, 0.5, -0.5, 2.0, 5.0, 10.0];

        for &v in &vals {
            let q = Fp6E3M2::from(v);
            let r: f32 = q.into();
            let err = (r - v).abs();

            println!("Fp6 roundtrip: v = {v:.8e}, q = {q}, r = {r:.8e}, err = {err:.3e}");

            assert!(r.is_finite() || r.is_nan());
            // coarse format; allow big relative error but make sure it’s not absurd
            let tol = 0.5 * v.abs().max(1.0);
            assert!(
                err <= tol,
                "Fp6 roundtrip too large: v={v}, r={r}, err={err}, tol={tol}"
            );
        }
    }

    #[test]
    fn fp6e3m2_division_preserves_scale() {
        let one = Fp6E3M2::from(1.0);
        let half = Fp6E3M2::from(0.5);

        assert_eq!(f32::from(one / one), 1.0);
        assert_eq!(f32::from(one / half), 2.0);
        assert_eq!(f32::from(half / one), 0.5);
        assert!(Fp6E3M2::from(0.0) < half);
        assert!(-half < Fp6E3M2::from(0.0));
    }
}
