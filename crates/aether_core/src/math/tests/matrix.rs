#[cfg(test)]
mod tests {
    use crate::{math::Matrix, math::Vector};

    #[test]
    fn test_matrix_add_sub_neg() {
        let a = Matrix {
            data: [[1.0_f64, 2.0], [3.0, 4.0]],
        };
        let b = Matrix {
            data: [[4.0_f64, 3.0], [2.0, 1.0]],
        };

        let c = a + b;
        assert_eq!(c.data, [[5.0, 5.0], [5.0, 5.0]]);

        let d = c - a;
        assert_eq!(d.data, b.data);

        let e = -a;
        assert_eq!(e.data, [[-1.0, -2.0], [-3.0, -4.0]]);
    }

    #[test]
    fn test_matrix_scalar_mul_div() {
        let m = Matrix {
            data: [[1.0_f64, -2.0], [3.0, -4.0]],
        };
        let s = 3.0;

        let mul = m * s;
        assert_eq!(mul.data, [[3.0, -6.0], [9.0, -12.0]]);

        let div = mul / s;
        assert_eq!(div.data, m.data);
    }

    #[test]
    fn test_matrix_mul() {
        let a = Matrix {
            data: [[1.0_f64, 2.0], [3.0, 4.0]],
        };
        let b = Matrix {
            data: [[2.0_f64, 0.0], [1.0, 2.0]],
        };

        let c = a * b;
        // c = [[1*2 + 2*1, 1*0 + 2*2], [3*2 + 4*1, 3*0 + 4*2]]
        //   = [[4, 4], [10, 8]]
        assert_eq!(c.data, [[4.0, 4.0], [10.0, 8.0]]);

        let c_ref = &a * &b;
        assert_eq!(c_ref.data, [[4.0, 4.0], [10.0, 8.0]]);
    }

    #[test]
    fn test_matrix_mul_4x4_fast_path() {
        let a = Matrix {
            data: [
                [1.0_f64, 2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0, 8.0],
                [9.0, 10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0, 16.0],
            ],
        };
        let b = Matrix {
            data: [
                [16.0_f64, 15.0, 14.0, 13.0],
                [12.0, 11.0, 10.0, 9.0],
                [8.0, 7.0, 6.0, 5.0],
                [4.0, 3.0, 2.0, 1.0],
            ],
        };

        let c = a * b;
        assert_eq!(
            c.data,
            [
                [80.0, 70.0, 60.0, 50.0],
                [240.0, 214.0, 188.0, 162.0],
                [400.0, 358.0, 316.0, 274.0],
                [560.0, 502.0, 444.0, 386.0],
            ]
        );
    }


    #[test]
    fn test_matrix_vector_mul() {
        let m = Matrix {
            data: [[1.0_f64, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
        };
        let v = Vector {
            data: [1.0_f64, 0.0, -1.0],
        };

        let r = m * v;
        // r = [
        //   1*1 + 2*0 + 3*(-1) = 1 + 0 - 3 = -2
        //   4*1 + 5*0 + 6*(-1) = 4 + 0 - 6 = -2
        //   7*1 + 8*0 + 9*(-1) = 7 + 0 - 9 = -2
        // ]
        assert_eq!(r.data, [-2.0, -2.0, -2.0]);

        let r_ref = &m * &v;
        assert_eq!(r_ref.data, [-2.0, -2.0, -2.0]);
    }

    #[test]
    fn test_matrix_vector_mul_4_fast_path() {
        let m = Matrix {
            data: [
                [1.0_f64, 2.0, 3.0, 4.0],
                [4.0, 3.0, 2.0, 1.0],
                [0.0, 1.0, 0.0, 1.0],
                [2.0, 0.0, 2.0, 0.0],
            ],
        };
        let v = Vector {
            data: [1.0_f64, -1.0, 2.0, 0.5],
        };

        let r = m * v;
        assert_eq!(r.data, [7.0, 5.5, -0.5, 6.0]);
    }


    #[cfg(feature = "bincode")]
    #[test]
    fn test_matrix_bincode_roundtrip() {
        let m = Matrix {
            data: [[1.0_f64, 2.5, -3.0], [4.0, 0.0, 6.25]],
        };

        let config = bincode::config::standard();
        let bytes = bincode::encode_to_vec(m, config).unwrap();
        let (decoded, _len): (Matrix<f64, 2, 3>, usize) =
            bincode::decode_from_slice(&bytes, config).unwrap();

        assert_eq!(decoded.data, m.data);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_matrix_serde_roundtrip() {
        let m = Matrix {
            data: [[1.0_f64, 2.5, -3.0], [4.0, 0.0, 6.25]],
        };

        let json = serde_json_core::to_string::<_, 256>(&m).unwrap();
        let (decoded, _): (Matrix<f64, 2, 3>, _) =
            serde_json_core::from_str(json.as_str()).unwrap();

        assert_eq!(decoded.data, m.data);
    }
}
