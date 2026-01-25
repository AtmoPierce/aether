// Matrix
/// Construct a statically typed matrix from literals
/// Example: `matrix![1.0, 2.0; 3.0, 4.0]` produces a Matrix<_, 2, 2>
#[macro_export]
macro_rules! matrix {
    // Accept input like: matrix![1.0, 2.0; 3.0, 4.0]
    ( $( [ $( $x:expr ),+ $(,)? ] ),+ $(,)? ) => {{
        use $crate::math::Matrix;
        Matrix {
            data: [
                $(
                    [ $( $x ),+ ],
                )+
            ]
        }
    }};
    // Alternative: Accept without inner brackets
    ( $( $( $x:expr ),+ );+ $(;)? ) => {{
        matrix![
            $(
                [ $( $x ),+ ]
            ),+
        ]
    }};
}
pub(crate) use matrix;

#[macro_export]
macro_rules! eye {
    ($n:expr) => {
        Matrix::identity($n)
    };
}
pub(crate) use eye;

#[macro_export]
macro_rules! zeros {
    ($m:expr, $n:expr) => {
        Matrix::zeros($m, $n)
    };
}
pub(crate) use zeros;

#[macro_export]
macro_rules! ones {
    ($m:expr, $n:expr) => {
        Matrix::ones($m, $n)
    };
}
pub(crate) use ones;

// Vectors
use crate::math::Vector;
use num_traits::Float;
impl<T, const N: usize> Vector<T, N>
where
    T: Float + Copy,
{
    pub fn zeros() -> Self {
        Self {
            data: [T::zero(); N],
        }
    }
    pub fn ones() -> Self {
        Self {
            data: [T::one(); N],
        }
    }
}
