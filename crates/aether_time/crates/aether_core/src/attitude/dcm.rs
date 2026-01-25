use super::euler::Euler;
use super::quaternion::Quaternion;
use crate::{coordinate::Cartesian, math::Matrix, reference_frame::ReferenceFrame};
use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Neg, Sub};
use num_traits::Float;

#[derive(Debug, Clone, Copy, Default)]
pub struct DirectionCosineMatrix<T: Float, From: ReferenceFrame, To: ReferenceFrame> {
    data: Matrix<T, 3, 3>, // still private
    _from: PhantomData<From>,
    _to: PhantomData<To>,
}

impl<T: Float, From: ReferenceFrame, To: ReferenceFrame> DirectionCosineMatrix<T, From, To> {
    pub fn new(m11: T, m12: T, m13: T, m21: T, m22: T, m23: T, m31: T, m32: T, m33: T) -> Self {
        let data = Matrix {
            data: [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]],
        };

        Self {
            data,
            _from: PhantomData,
            _to: PhantomData,
        }
    }

    pub fn from_matrix(data: Matrix<T, 3, 3>) -> Self {
        Self {
            data,
            _from: PhantomData::<From>,
            _to: PhantomData::<To>,
        }
    }

    pub fn as_matrix(&self) -> &Matrix<T, 3, 3> {
        &self.data
    }

    pub fn rotate_x(angle: T) -> Self {
        let (c, s) = (angle.cos(), angle.sin());
        let data: [[T; 3]; 3] = [
            [T::one(), T::zero(), T::zero()],
            [T::zero(), c, s],
            [T::zero(), -s, c],
        ];
        Self {
            data: Matrix { data },
            _from: PhantomData,
            _to: PhantomData,
        }
    }

    pub fn rotate_y(angle: T) -> Self {
        let (c, s) = (angle.cos(), angle.sin());
        let data = [
            [c, T::zero(), -s],
            [T::zero(), T::one(), T::zero()],
            [s, T::zero(), c],
        ];
        Self {
            data: Matrix { data },
            _from: PhantomData,
            _to: PhantomData,
        }
    }

    pub fn rotate_z(angle: T) -> Self {
        let (c, s) = (angle.cos(), angle.sin());

        let data = [
            [c, s, T::zero()],
            [-s, c, T::zero()],
            [T::zero(), T::zero(), T::one()],
        ];
        Self {
            data: Matrix { data },
            _from: PhantomData,
            _to: PhantomData,
        }
    }
}

impl<T, A, B, C> Mul<DirectionCosineMatrix<T, A, B>> for DirectionCosineMatrix<T, B, C>
where
    T: Float,
    A: ReferenceFrame,
    B: ReferenceFrame,
    C: ReferenceFrame,
{
    type Output = DirectionCosineMatrix<T, A, C>;

    fn mul(self, rhs: DirectionCosineMatrix<T, A, B>) -> Self::Output {
        let lhs: Matrix<T, 3, 3> = self.data; // B → C
        let rhs: Matrix<T, 3, 3> = rhs.data; // A → B
        let data: Matrix<T, 3, 3> = lhs * rhs; // A → C
        DirectionCosineMatrix::from_matrix(data)
    }
}

// From
impl<T, A: ReferenceFrame, B: ReferenceFrame> From<Euler<T>> for DirectionCosineMatrix<T, A, B>
where
    T: Float + Mul<Output = T> + Add<Output = T> + Copy + Default,
{
    fn from(euler: Euler<T>) -> Self {
        // ZYX Convention
        let [phi, theta, psi] = euler.data.data;
        let rx: DirectionCosineMatrix<T, _, B> = DirectionCosineMatrix::rotate_x(phi);
        let ry: DirectionCosineMatrix<T, _, B> = DirectionCosineMatrix::rotate_y(theta);
        let rz: DirectionCosineMatrix<T, _, B> = DirectionCosineMatrix::rotate_z(psi);
        let c = rz * ry * rx;

        return c;
    }
}

impl<T: Float, A: ReferenceFrame, B: ReferenceFrame> From<Quaternion<T>>
    for DirectionCosineMatrix<T, A, B>
{
    fn from(q: Quaternion<T>) -> Self {
        let s = q.data.data[0];
        let i = q.data.data[1];
        let j = q.data.data[2];
        let k = q.data.data[3];

        let two = T::one() + T::one();
        let sp2 = s.powf(two);
        let ip2 = i.powf(two);
        let jp2 = j.powf(two);
        let kp2 = k.powf(two);

        let data = [
            [
                sp2 + ip2 - jp2 - kp2,
                two * (i * j + s * k),
                two * (i * k - s * j),
            ],
            [
                two * (i * j - s * k),
                sp2 - ip2 + jp2 - kp2,
                two * (j * k + s * i),
            ],
            [
                two * (i * k + s * j),
                two * (j * k - s * i),
                sp2 - ip2 - jp2 + kp2,
            ],
        ];
        Self {
            data: Matrix { data },
            _from: PhantomData,
            _to: PhantomData,
        }
    }
}

impl<T: Float, A: ReferenceFrame, B: ReferenceFrame> From<&Matrix<T, 3, 3>>
    for DirectionCosineMatrix<T, A, B>
{
    fn from(m: &Matrix<T, 3, 3>) -> Self {
        Self {
            data: *m,
            _from: PhantomData,
            _to: PhantomData,
        }
    }
}

// Cartesian Handling
impl<T: Float, From: ReferenceFrame, To: ReferenceFrame> Mul<Cartesian<T, From>>
    for DirectionCosineMatrix<T, From, To>
{
    type Output = Cartesian<T, To>;

    fn mul(self, rhs: Cartesian<T, From>) -> Self::Output {
        Cartesian {
            data: self.data * rhs.data,
            _reference_frame: PhantomData,
        }
    }
}

// Behavior
impl<T: Float + Copy, From: ReferenceFrame, To: ReferenceFrame> DirectionCosineMatrix<T, From, To> {
    pub fn m11(&self) -> T {
        self.data[(0, 0)]
    }
    pub fn m12(&self) -> T {
        self.data[(0, 1)]
    }
    pub fn m13(&self) -> T {
        self.data[(0, 2)]
    }
    pub fn m21(&self) -> T {
        self.data[(1, 0)]
    }
    pub fn m22(&self) -> T {
        self.data[(1, 1)]
    }
    pub fn m23(&self) -> T {
        self.data[(1, 2)]
    }
    pub fn m31(&self) -> T {
        self.data[(2, 0)]
    }
    pub fn m32(&self) -> T {
        self.data[(2, 1)]
    }
    pub fn m33(&self) -> T {
        self.data[(2, 2)]
    }
    pub fn transpose(&self) -> DirectionCosineMatrix<T, To, From> {
        let mut result = Matrix::<T, 3, 3>::zeros();
        for i in 0..3 {
            for j in 0..3 {
                result[(i, j)] = self.data[(j, i)];
            }
        }
        DirectionCosineMatrix {
            data: result,
            _from: PhantomData,
            _to: PhantomData,
        }
    }
}

#[cfg(feature = "std")]
impl<T: Float + std::fmt::Display, A: ReferenceFrame, B: ReferenceFrame> std::fmt::Display
    for DirectionCosineMatrix<T, A, B>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.data)
    }
}

// Tests
#[cfg(test)]
#[path = "tests/dcm_tests.rs"]
mod dcm_tests;
