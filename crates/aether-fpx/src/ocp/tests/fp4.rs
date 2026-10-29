mod quantities {
    use crate::ocp::fp4::Fp4E2M1;

    fn approx_eq(a: f32, b: f32, eps: f32) -> bool {
        (a - b).abs() <= eps
    }

    // -------------------------
    // FP4 E2M1 tests
    // -------------------------

    #[test]
    fn fp4e2m1_zero_and_sign() {
        use crate::ocp::fp4::Fp4E2M1;

        // +0.0: bits = 0x0
        // NaN sentinel: bits = 0x8
        // A clearly negative finite value: bits = 0xF (~= -6.0)
        let pz = Fp4E2M1::from_bits(0x0);
        let neg = Fp4E2M1::from_bits(0xF);

        let pzf: f32 = pz.into();
        let negf: f32 = neg.into();

        println!("Fp4 +0 bits=0x00 => {}, f32={}", pz, pzf);
        println!("Fp4 neg bits=0x0F => {}, f32={}", neg, negf);

        // Only +0 is representable as zero in the encoding
        assert_eq!(pzf.to_bits(), 0.0f32.to_bits());
        // Check if negative finite value
        assert!(negf.is_finite() && negf < 0.0);
    }

    #[test]
    fn fp4e2m1_full_table() {
        // Exact value table for E2M1 (bias=1, s e e m layout)
        let cases: &[(u8, f32)] = &[
            (0x0, 0.0),
            (0x1, 0.5),
            (0x2, 1.0),
            (0x3, 1.5),
            (0x4, 2.0),
            (0x5, 3.0),
            (0x6, 4.0),
            (0x7, 6.0),
            (0x8, -0.0),
            (0x9, -0.5),
            (0xA, -1.0),
            (0xB, -1.5),
            (0xC, -2.0),
            (0xD, -3.0),
            (0xE, -4.0),
            (0xF, -6.0),
        ];

        for &(bits, expected) in cases {
            let x = Fp4E2M1::from_bits(bits);
            let v: f32 = x.into();
            let err = (v - expected).abs();

            println!(
                "Fp4 bits=0x{:X} => {} (f32={:.8e}), expected={:.8e}, err={:.3e}",
                bits, x, v, expected, err
            );

            // These are exactly representable by construction; allow tiny FP noise.
            assert!(
                approx_eq(v, expected, 1e-6),
                "decode mismatch: bits=0x{bits:X}, got {v}, expected {expected}"
            );
        }
    }

    #[test]
    fn fp4e2m1_roundtrip_simple_values() {
        // A mix of exactly representable and in-between values.
        let vals = [
            0.0f32, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, -0.25, -0.5, -1.0, -3.0, -6.0,
        ];

        for &v in &vals {
            let q = Fp4E2M1::from(v);
            let r: f32 = q.into();
            let err = (r - v).abs();

            println!("Fp4 roundtrip: v={v:.8e}, q={q}, r={r:.8e}, err={err:.3e}");

            // Very coarse format; allow big relative error for tiny |v|.
            let tol = 0.5 * v.abs().max(1.0);
            assert!(
                err <= tol,
                "Fp4 roundtrip too large: v={v}, r={r}, err={err}, tol={tol}"
            );
        }
    }

    #[test]
    fn fp4e2m1_division_preserves_scale() {
        let one = Fp4E2M1::from(1.0);
        let half = Fp4E2M1::from(0.5);

        assert_eq!(f32::from(one / one), 1.0);
        assert_eq!(f32::from(one / half), 2.0);
        assert_eq!(f32::from(half / one), 0.5);
        assert!(Fp4E2M1::from(0.0) < half);
        assert!(-half < Fp4E2M1::from(0.0));
    }
}
