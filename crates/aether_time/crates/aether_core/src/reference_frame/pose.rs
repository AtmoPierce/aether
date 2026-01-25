use crate::attitude::DirectionCosineMatrix;
use crate::coordinate::Cartesian;
use crate::math::Matrix;
use crate::reference_frame::{Assembly, ReferenceFrame};
use num_traits::Float;

#[derive(Debug, Clone, Copy)]
pub struct Pose<T: Float, From: ReferenceFrame> {
    /// DCM from the solid's Local frame to the Assembly frame (Local → Assembly)
    pub R_local_to_assembly: DirectionCosineMatrix<T, From, Assembly<T>>,
    /// Center-of-mass position in the Assembly frame
    pub r_cm_in_assembly: Cartesian<T, Assembly<T>>,
}

impl<T: Float, From: ReferenceFrame> Pose<T, From> {
    pub fn identity() -> Self {
        Self {
            R_local_to_assembly: DirectionCosineMatrix::from_matrix(Matrix::<T, 3, 3>::identity()),
            r_cm_in_assembly: Cartesian::zero(),
        }
    }
}
