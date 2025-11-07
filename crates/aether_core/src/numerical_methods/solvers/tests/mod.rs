#[cfg(test)]
mod lu_tests {
    use crate::*;
    use crate::math::{Matrix, Vector};

    #[inline]
    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() <= eps
    }

    #[test]
    fn lu_solve_single_rhs_3x3() {
        // A x = b
        let a: Matrix<f64, 3, 3> = Matrix::new([
            [ 3.0,  2.0, -1.0],
            [ 2.0, -2.0,  4.0],
            [-1.0,  0.5, -1.0],
        ]);
        let b = Vector { data: [1.0, -2.0, 0.0] };

        let x = a.solve(&b).expect("LU failed (singular?)");
        let expected = [1.0, -2.0, -2.0];
        for i in 0..3 {
            assert!(approx_eq(x[i], expected[i], 1e-12), "x[{i}] = {}, expected {}", x[i], expected[i]);
        }
        println!("Single-RHS LU solve passed");
    }

    #[test]
    fn lu_solve_multi_rhs() {
        // A X = B (two RHS)
        let a: Matrix<f64, 3, 3> = Matrix::new([
            [ 3.0,  2.0, -1.0],
            [ 2.0, -2.0,  4.0],
            [-1.0,  0.5, -1.0],
        ]);
        let b: Matrix<f64, 3, 2> = Matrix::new([
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ]);

        let x = a.solve_multi(&b).expect("LU failed (singular?)");

        // Expected from numpy.linalg.solve:
        // [[-2.0,     -2.5    ],
        //  [ 5.333..., 6.0    ],
        //  [ 3.666..., 4.5    ]]
        let expected = [
            [-2.0,        -2.5],
            [ 5.333333333333333, 6.0],
            [ 3.6666666666666665, 4.5],
        ];

        for r in 0..3 {
            for c in 0..2 {
                assert!(
                    approx_eq(x[r][c], expected[r][c], 1e-12),
                    "x[{r}][{c}] = {}, expected {}", x[r][c], expected[r][c]
                );
            }
        }
        println!("Multi-RHS LU solve passed");
    }

    #[test]
    fn lu_pivoting_works() {
        // Pivot required (zero top-left). A = [[0,1],[1,0]]
        let a: Matrix<f64, 2, 2> = Matrix::new([[0.0, 1.0],
                                                [1.0, 0.0]]);
        let b = Vector { data: [1.0, 2.0] };

        let x = a.solve(&b).expect("LU failed");
        let expected = [2.0, 1.0]; // since Ax = b => x = [2,1]
        assert!((x[0] - expected[0]).abs() < 1e-12);
        assert!((x[1] - expected[1]).abs() < 1e-12);
        println!("x = [{}, {}], expected [{}, {}]", x[0], x[1], expected[0], expected[1]);
    }
}

#[cfg(test)]
mod cholesky_tests {
    use crate::*;
    use crate::math::{Matrix, Vector};

    #[inline]
    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() <= eps
    }

    #[test]
    fn cholesky_solve_spd() {
        // SPD matrix
        let q: Matrix<f64, 3, 3> = Matrix::new([
            [4.0, 1.0, 0.0],
            [1.0, 3.0, 0.0],
            [0.0, 0.0, 2.0],
        ]);
        let rhs = Vector { data: [1.0, 2.0, 3.0] };

        let x = q.cholesky_solve(&rhs).expect("Cholesky failed (not SPD?)");

        // Expected from numpy.linalg.solve:
        // [0.09090909, 0.63636364, 1.5]
        let expected = [0.090_909_090_909_090_91, 0.636_363_636_363_636_4, 1.5];

        for i in 0..3 {
            assert!(approx_eq(x[i], expected[i], 1e-12), "x[{i}] = {}, expected {}", x[i], expected[i]);
            println!("x[{i}] = {}, expected {}", x[i], expected[i]);
        }
    }

    #[test]
    fn cholesky_rejects_non_spd() {
        // Symmetric but not positive definite: eigenvalues are 3 and -1
        let a: Matrix<f64, 2, 2> = Matrix::new([
            [1.0, 2.0],
            [2.0, 1.0],
        ]);

        assert!(a.cholesky_lower().is_none(), "Should reject non-SPD matrix");
        let b = Vector { data: [1.0, 1.0] };
        assert!(a.cholesky_solve(&b).is_none(), "Solve should fail for non-SPD");
        println!("Non-SPD matrix correctly rejected");
    }
}

#[cfg(test)]
mod newton_raphson_tests {
    use crate::numerical_methods::solvers::newton::NewtonRaphson;
    #[test]
    fn test_newton_raphson() {
        let solver = NewtonRaphson::new(1e-7, 100);
        // Test: f(x) = x^2 - 2, root at sqrt(2)
        let f = |x: f64| x * x - 2.0;
        let df = |x: f64| 2.0 * x;
        let initial_guess = 1.0;
        let root = solver.solve(initial_guess, f, df).expect("Failed to find root");
        let expected = 2f64.sqrt();
        assert!((root - expected).abs() < 1e-7, "Root found: {}, expected: {}", root, expected);
        println!("Newton-Raphson test passed: root = {}, expected = {}", root, expected);
    }
}