#[cfg(test)]
mod tests {
    use crate::{math::Vector, math::Matrix};

    #[test]
    fn test_matrix_add_sub_neg() {
        let a = Matrix {
            data: [
                [1.0_f64, 2.0],
                [3.0, 4.0],
            ],
        };
        let b = Matrix {
            data: [
                [4.0_f64, 3.0],
                [2.0, 1.0],
            ],
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
            data: [
                [1.0_f64, -2.0],
                [3.0, -4.0],
            ],
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
            data: [
                [1.0_f64, 2.0],
                [3.0, 4.0],
            ],
        };
        let b = Matrix {
            data: [
                [2.0_f64, 0.0],
                [1.0, 2.0],
            ],
        };

        let c = a * b;
        // c = [[1*2 + 2*1, 1*0 + 2*2], [3*2 + 4*1, 3*0 + 4*2]]
        //   = [[4, 4], [10, 8]]
        assert_eq!(c.data, [[4.0, 4.0], [10.0, 8.0]]);
    }

    #[test]
    fn test_matrix_vector_mul() {
        let m = Matrix {
            data: [
                [1.0_f64, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
            ],
        };
        let v = Vector { data: [1.0_f64, 0.0, -1.0] };

        let r = m * v;
        // r = [
        //   1*1 + 2*0 + 3*(-1) = 1 + 0 - 3 = -2
        //   4*1 + 5*0 + 6*(-1) = 4 + 0 - 6 = -2
        //   7*1 + 8*0 + 9*(-1) = 7 + 0 - 9 = -2
        // ]
        assert_eq!(r.data, [-2.0, -2.0, -2.0]);
    }
}
