mod quantities {
    use crate::ocp::exp2i::*;
    use crate::ocp::fp8::*;
    use crate::ocp::scale::*;

    fn approx_eq(a: f32, b: f32, eps: f32) -> bool {
        (a - b).abs() <= eps
    }

    // Helper for tests: print error then assert.
    fn assert_approx(label: &str, a: f32, b: f32, eps: f32) {
        let err = (a - b).abs();
        println!("{label}: a = {a:.8e}, b = {b:.8e}, err = {err:.3e}, eps = {eps:.3e}");
        assert!(
            err <= eps,
            "{label} failed: a = {a}, b = {b}, err = {err}, eps = {eps}"
        );
    }

    // -------------------------
    // E8M0 tests
    // -------------------------

    #[test]
    fn e8m0_basic_powers_of_two() {
        let s0 = E8M0::new(0);
        let s1 = E8M0::new(1);
        let sm1 = E8M0::new(-1);
        let s10 = E8M0::new(10);

        println!("s0  = {}", s0);
        println!("s1  = {}", s1);
        println!("sm1 = {}", sm1);
        println!("s10 = {}", s10);

        assert_eq!(s0.scale_f32(), 1.0);
        assert_eq!(s1.scale_f32(), 2.0);
        assert_eq!(sm1.scale_f32(), 0.5);
        assert_eq!(s10.scale_f32(), 1024.0);
    }

    #[test]
    fn e8m0_nan_sentinel_behaves() {
        let n = E8M0::NAN;
        println!("E8M0 NAN = {}", n);

        assert!(n.is_nan());
        let s = n.scale_f32();
        assert!(s.is_nan());
        assert_eq!(n.exponent(), None);
        assert_eq!(n.bits(), i8::MIN);
    }

    // -------------------------
    // exp2i tests
    // -------------------------

    #[test]
    fn exp2i_exact_powers() {
        for &k in &[0, 1, 10, -1] {
            let v = exp2i(k);
            println!("exp2i({k}) = {v}");
        }

        assert_eq!(exp2i(0), 1.0);
        assert_eq!(exp2i(1), 2.0);
        assert_eq!(exp2i(10), 1024.0);
        assert_eq!(exp2i(-1), 0.5);
    }

    // -------------------------
    // FP8 E5M2 tests
    // -------------------------

    #[test]
    fn fp8e5m2_zero_and_sign() {
        let pz = Fp8E5M2::from_bits(0x00);
        let nz = Fp8E5M2::from_bits(0x80);
        let pzf: f32 = pz.into();
        let nzf: f32 = nz.into();

        println!("E5M2 pz = {}, f32 = {}", pz, pzf);
        println!("E5M2 nz = {}, f32 = {}", nz, nzf);

        assert_eq!(pzf, 0.0);
        assert_eq!(nzf, -0.0);
        assert_eq!(pzf.to_bits(), 0.0f32.to_bits());
        assert_eq!(nzf.to_bits(), (-0.0f32).to_bits());
    }

    #[test]
    fn fp8e5m2_min_normal_and_subnormal_decode() {
        let sub = Fp8E5M2::MIN_SUBNORMAL_POS;
        let nor = Fp8E5M2::MIN_NORMAL_POS;

        let sub_f: f32 = sub.into();
        let nor_f: f32 = nor.into();

        println!(
            "E5M2 MIN_SUBNORMAL_POS bits=0x{:02X} => {} (f32 = {:.8e})",
            sub.to_bits(),
            sub,
            sub_f
        );
        println!(
            "E5M2 MIN_NORMAL_POS    bits=0x{:02X} => {} (f32 = {:.8e})",
            nor.to_bits(),
            nor,
            nor_f
        );

        // With bias=15, mant bits=2:
        // sub: e=0, m=1 => (1/4)*2^(1-15) = 2^-16 = 1/65536
        assert_approx("fp8e5m2_min_subnormal", sub_f, 1.0 / 65536.0, 1e-10);

        // min normal: e=1, m=0 => 1 * 2^(1-15) = 2^-14 = 1/16384
        assert_approx("fp8e5m2_min_normal", nor_f, 1.0 / 16384.0, 1e-10);
    }

    #[test]
    fn fp8e5m2_infinity_and_nan_decode() {
        let pinf = Fp8E5M2::POS_INF;
        let ninf = Fp8E5M2::NEG_INF;
        let nan = Fp8E5M2::CANONICAL_NAN;

        let pinf_f: f32 = pinf.into();
        let ninf_f: f32 = ninf.into();
        let nan_f: f32 = nan.into();

        println!("E5M2 POS_INF = {}, f32 = {}", pinf, pinf_f);
        println!("E5M2 NEG_INF = {}, f32 = {}", ninf, ninf_f);
        println!("E5M2 NAN     = {}, f32 = {}", nan, nan_f);

        assert!(pinf_f.is_infinite() && pinf_f.is_sign_positive());
        assert!(ninf_f.is_infinite() && ninf_f.is_sign_negative());
        assert!(nan_f.is_nan());
    }

    #[test]
    fn fp8e5m2_max_finite_decode() {
        let max = Fp8E5M2::MAX_FINITE;
        let v: f32 = max.into();

        println!(
            "E5M2 MAX_FINITE bits=0x{:02X} => {} (f32 = {:.8e})",
            max.to_bits(),
            max,
            v
        );

        // Expected: (1 + 3/4) * 2^15 = 1.75 * 32768 = 57344
        assert_approx("fp8e5m2_max_finite", v, 57_344.0, 1.0);
    }

    #[test]
    fn fp8e5m2_roundtrip_simple_values() {
        let vals = [0.0f32, 1.0, -1.0, 0.5, -0.5, 10.0];

        for &v in &vals {
            let q = Fp8E5M2::from(v);
            let r: f32 = q.into();
            let err = (r - v).abs();

            println!("E5M2 roundtrip: v = {v:.8e}, q = {q}, r = {r:.8e}, err = {err:.3e}");

            assert!(r.is_finite());
            let tol = 0.5 * v.abs().max(1.0);
            assert!(
                err <= tol,
                "E5M2 roundtrip too large: v={v}, r={r}, err={err}, tol={tol}"
            );
        }
    }

    // -------------------------
    // FP8 E4M3 tests
    // -------------------------

    #[test]
    fn fp8e4m3_zero_and_sign() {
        let pz = Fp8E4M3::from_bits(0x00);
        let nz = Fp8E4M3::from_bits(0x80);
        let pzf: f32 = pz.into();
        let nzf: f32 = nz.into();

        println!("E4M3 pz = {}, f32 = {}", pz, pzf);
        println!("E4M3 nz = {}, f32 = {}", nz, nzf);

        assert_eq!(pzf, 0.0);
        assert_eq!(nzf, -0.0);
        assert_eq!(pzf.to_bits(), 0.0f32.to_bits());
        assert_eq!(nzf.to_bits(), (-0.0f32).to_bits());
    }

    #[test]
    fn fp8e4m3_subnormal_and_min_normal_decode() {
        let sub = Fp8E4M3::MIN_SUBNORMAL_POS; // bits=0x01
        let nor = Fp8E4M3::MIN_NORMAL_POS; // bits=0x08

        let sub_f: f32 = sub.into();
        let nor_f: f32 = nor.into();

        println!(
            "E4M3 MIN_SUBNORMAL_POS bits=0x{:02X} => {} (f32 = {:.8e})",
            sub.to_bits(),
            sub,
            sub_f
        );
        println!(
            "E4M3 MIN_NORMAL_POS    bits=0x{:02X} => {} (f32 = {:.8e})",
            nor.to_bits(),
            nor,
            nor_f
        );

        // bias=7, mant bits=3:
        // sub: e=0,m=1 => (1/8)*2^(1-7) = 2^-9 = 1/512
        assert_approx("fp8e4m3_min_subnormal", sub_f, 1.0 / 512.0, 1e-6);

        // min normal: e=1,m=0 => 1 * 2^(1-7) = 2^-6 = 1/64
        assert_approx("fp8e4m3_min_normal", nor_f, 1.0 / 64.0, 1e-6);
    }

    #[test]
    fn fp8e4m3_max_finite_decode() {
        let max = Fp8E4M3::MAX_FINITE;
        let v: f32 = max.into();

        println!(
            "E4M3 MAX_FINITE bits=0x{:02X} => {} (f32 = {:.8e})",
            max.to_bits(),
            max,
            v
        );

        assert_approx("fp8e4m3_max_finite", v, 448.0, 0.5);
    }

    #[test]
    fn fp8e4m3_nan_decode() {
        let nan = Fp8E4M3::CANONICAL_NAN;
        let v: f32 = nan.into();

        println!(
            "E4M3 NAN bits=0x{:02X} => {} (f32 = {v})",
            nan.to_bits(),
            nan
        );

        assert!(v.is_nan());
    }

    #[test]
    fn fp8e4m3_roundtrip_simple_values() {
        let vals = [0.0f32, 1.0, -1.0, 0.5, -0.5, 10.0];

        for &v in &vals {
            let q = Fp8E4M3::from(v);
            let r: f32 = q.into();
            let err = (r - v).abs();

            println!("E4M3 roundtrip: v = {v:.8e}, q = {q}, r = {r:.8e}, err = {err:.3e}");

            assert!(r.is_finite() || r.is_nan());
            let tol = 0.5 * v.abs().max(1.0);
            assert!(
                err <= tol,
                "E4M3 roundtrip too large: v={v}, r={r}, err={err}, tol={tol}"
            );
        }
    }
}

mod operations {
    use crate::ocp::fp8::{Fp8E4M3, Fp8E5M2};

    #[test]
    fn fp8e5m2_division_preserves_scale() {
        let one = Fp8E5M2::from(1.0);
        let half = Fp8E5M2::from(0.5);

        assert_eq!(f32::from(one / one), 1.0);
        assert_eq!(f32::from(one / half), 2.0);
        assert_eq!(f32::from(half / one), 0.5);
        assert!(Fp8E5M2::from(0.0) < half);
        assert!(-half < Fp8E5M2::from(0.0));
    }

    #[test]
    fn fp8e4m3_division_preserves_scale() {
        let one = Fp8E4M3::from(1.0);
        let half = Fp8E4M3::from(0.5);

        assert_eq!(f32::from(one / one), 1.0);
        assert_eq!(f32::from(one / half), 2.0);
        assert_eq!(f32::from(half / one), 0.5);
        assert!(Fp8E4M3::from(0.0) < half);
        assert!(-half < Fp8E4M3::from(0.0));
    }
}
