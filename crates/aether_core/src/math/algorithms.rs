use crate::math::matrix::Matrix;
use crate::math::vector::Vector;
use crate::real::Real;

pub trait MatrixAlgorithms<T: Real, const M: usize, const N: usize> {
    fn mul_matrix_strided<const P: usize>(&self, rhs: &Matrix<T, N, P>) -> Matrix<T, M, P>;
    fn mul_vector_unrolled(&self, rhs: &Vector<T, N>) -> Vector<T, M>;
}

impl<T: Real, const M: usize, const N: usize> MatrixAlgorithms<T, M, N> for Matrix<T, M, N> {
    #[inline(always)]
    fn mul_matrix_strided<const P: usize>(&self, rhs: &Matrix<T, N, P>) -> Matrix<T, M, P> {
        const BLOCK: usize = 16;

        let mut result = Matrix {
            data: [[T::ZERO; P]; M],
        };

        for i in 0..M {
            let mut k0 = 0;
            while k0 < N {
                let k_max = (k0 + BLOCK).min(N);

                let mut j0 = 0;
                while j0 < P {
                    let j_max = (j0 + BLOCK).min(P);

                    for k in k0..k_max {
                        let a = self.data[i][k];
                        for j in j0..j_max {
                            result.data[i][j] = result.data[i][j] + a * rhs.data[k][j];
                        }
                    }

                    j0 += BLOCK;
                }

                k0 += BLOCK;
            }
        }

        result
    }

    #[inline(always)]
    fn mul_vector_unrolled(&self, rhs: &Vector<T, N>) -> Vector<T, M> {
        let mut result = Vector {
            data: [T::ZERO; M],
        };

        for i in 0..M {
            let row = &self.data[i];
            let mut sum = T::ZERO;
            for j in 0..N {
                sum = sum + row[j] * rhs.data[j];
            }
            result.data[i] = sum;
        }

        result
    }
}

pub trait VectorAlgorithms<T: Real, const N: usize> {
    fn dot_generic(&self, rhs: &Vector<T, N>) -> T;
    fn norm_generic(&self) -> T;
}

impl<T: Real, const N: usize> VectorAlgorithms<T, N> for Vector<T, N> {
    #[inline]
    fn dot_generic(&self, rhs: &Vector<T, N>) -> T {
        let mut acc = T::ZERO;
        for i in 0..N {
            acc = acc + self.data[i] * rhs.data[i];
        }
        acc
    }

    #[inline]
    fn norm_generic(&self) -> T {
        let mut acc = T::ZERO;
        for &x in &self.data {
            acc = acc + x * x;
        }
        acc.sqrt()
    }
}
