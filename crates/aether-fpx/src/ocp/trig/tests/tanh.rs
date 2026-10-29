#[cfg(test)]
mod tests {
    use super::*;
    use crate::fp8::{Fp8E4M3, Fp8E5M2};
    use crate::ocp::trig::*;

    #[inline]
    fn abs_f32(x: f32) -> f32 {
        if x < 0.0 {
            -x
        } else {
            x
        }
    }

    #[test]
    fn tanh_e5m2_special_cases() {
        let pinf = Fp8E5M2::POS_INF;
        let ninf = Fp8E5M2::NEG_INF;
        let nan = Fp8E5M2::CANONICAL_NAN;

        assert_eq!(tanh_e5m2(pinf).to_bits(), Fp8E5M2::from(1.0f32).to_bits());
        assert_eq!(tanh_e5m2(ninf).to_bits(), Fp8E5M2::from(-1.0f32).to_bits());
        assert_eq!(tanh_e5m2(nan).to_bits(), Fp8E5M2::CANONICAL_NAN.to_bits());
    }

    #[test]
    fn tanh_e4m3_special_cases() {
        let nan_pos = Fp8E4M3::CANONICAL_NAN;
        let nan_neg = Fp8E4M3::from_bits(0xFF);

        assert_eq!(
            tanh_e4m3(nan_pos).to_bits(),
            Fp8E4M3::CANONICAL_NAN.to_bits()
        );
        assert_eq!(
            tanh_e4m3(nan_neg).to_bits(),
            Fp8E4M3::CANONICAL_NAN.to_bits()
        );
    }

    #[test]
    fn tanh_e5m2_threshold_behavior() {
        for bits in 0u8..=255 {
            let x = Fp8E5M2::from_bits(bits);
            let xf = f32::from(x);
            if !xf.is_finite() {
                continue;
            }

            let y = tanh_e5m2(x);
            let ax = abs_f32(xf);
            if ax < E5M2_T1 {
                assert_eq!(
                    y.to_bits(),
                    x.to_bits(),
                    "expected passthrough for bits=0x{bits:02X}"
                );
            }
            if ax > E5M2_T2 {
                let expected = if xf.is_sign_negative() {
                    Fp8E5M2::from(-1.0f32)
                } else {
                    Fp8E5M2::from(1.0f32)
                };
                assert_eq!(
                    y.to_bits(),
                    expected.to_bits(),
                    "expected saturation for bits=0x{bits:02X}"
                );
            }
        }
    }

    #[test]
    fn tanh_e4m3_threshold_behavior() {
        for bits in 0u8..=255 {
            let x = Fp8E4M3::from_bits(bits);
            let xf = f32::from(x);
            if xf.is_nan() {
                continue;
            }

            let y = tanh_e4m3(x);
            let ax = abs_f32(xf);
            if ax < E4M3_T1 {
                assert_eq!(
                    y.to_bits(),
                    x.to_bits(),
                    "expected passthrough for bits=0x{bits:02X}"
                );
            }
            if ax > E4M3_T2 {
                let expected = if xf.is_sign_negative() {
                    Fp8E4M3::from(-1.0f32)
                } else {
                    Fp8E4M3::from(1.0f32)
                };
                assert_eq!(
                    y.to_bits(),
                    expected.to_bits(),
                    "expected saturation for bits=0x{bits:02X}"
                );
            }
        }
    }

    #[test]
    fn tanh_e5m2_is_odd_and_bounded() {
        for bits in 0u8..=255 {
            let x = Fp8E5M2::from_bits(bits);
            let xf = f32::from(x);
            if !xf.is_finite() {
                continue;
            }

            let y = tanh_e5m2(x);
            let yv = f32::from(y);
            assert!(yv.is_finite());
            assert!(abs_f32(yv) <= 1.0 + 1e-6);

            let x_neg = Fp8E5M2::from_bits(bits ^ 0x80);
            let x_neg_f = f32::from(x_neg);
            if !x_neg_f.is_finite() {
                continue;
            }
            let y_neg = tanh_e5m2(x_neg);

            let y_mag = y.to_bits() & 0x7F;
            let y_neg_mag = y_neg.to_bits() & 0x7F;
            assert_eq!(
                y_mag, y_neg_mag,
                "odd symmetry magnitude mismatch at bits=0x{bits:02X}"
            );
            assert_eq!((y.to_bits() >> 7) ^ (y_neg.to_bits() >> 7), 1);
        }
    }

    #[test]
    fn tanh_e4m3_is_odd_and_bounded() {
        for bits in 0u8..=255 {
            let x = Fp8E4M3::from_bits(bits);
            let xf = f32::from(x);
            if xf.is_nan() {
                continue;
            }

            let y = tanh_e4m3(x);
            let yv = f32::from(y);
            if yv.is_nan() {
                continue;
            }
            assert!(abs_f32(yv) <= 1.0 + 1e-6);

            let x_neg = Fp8E4M3::from_bits(bits ^ 0x80);
            let x_neg_f = f32::from(x_neg);
            if x_neg_f.is_nan() {
                continue;
            }
            let y_neg = tanh_e4m3(x_neg);
            let y_neg_f = f32::from(y_neg);
            if y_neg_f.is_nan() {
                continue;
            }

            let y_mag = y.to_bits() & 0x7F;
            let y_neg_mag = y_neg.to_bits() & 0x7F;
            assert_eq!(
                y_mag, y_neg_mag,
                "odd symmetry magnitude mismatch at bits=0x{bits:02X}"
            );
            assert_eq!((y.to_bits() >> 7) ^ (y_neg.to_bits() >> 7), 1);
        }
    }

    #[test]
    fn tanh_e5m2_approximation_is_reasonable() {
        for bits in 0u8..=255 {
            let x = Fp8E5M2::from_bits(bits);
            let xf = f32::from(x);
            if !xf.is_finite() {
                continue;
            }

            let y = f32::from(tanh_e5m2(x));
            let ref_y = xf.tanh();
            let err = abs_f32(y - ref_y);
            assert!(
                err <= 0.35,
                "err={err} bits=0x{bits:02X} x={xf} y={y} ref={ref_y}"
            );
        }
    }

    #[test]
    fn tanh_e4m3_approximation_is_reasonable() {
        for bits in 0u8..=255 {
            let x = Fp8E4M3::from_bits(bits);
            let xf = f32::from(x);
            if xf.is_nan() {
                continue;
            }

            let y = f32::from(tanh_e4m3(x));
            if y.is_nan() {
                continue;
            }

            let ref_y = xf.tanh();
            let err = abs_f32(y - ref_y);
            assert!(
                err <= 0.5,
                "err={err} bits=0x{bits:02X} x={xf} y={y} ref={ref_y}"
            );
        }
    }
}
