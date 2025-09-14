pub mod euler;
pub mod dcm;
pub mod quaternion;
pub mod rotation;
pub mod determination;
pub use euler::Euler;
pub use dcm::DirectionCosineMatrix;
pub use quaternion::Quaternion;
pub use rotation::Rotation;
#[cfg(test)]
#[path = "tests/mod.rs"]
pub mod tests;

pub mod attitude{
    use num_traits::Float;
    use crate::reference_frame::ReferenceFrame;

    use super::dcm::DirectionCosineMatrix;
    use super::euler::Euler;
    use super::quaternion::Quaternion;
    pub enum Rotation<T, From: ReferenceFrame, To: ReferenceFrame> 
        where T: Float
            {
            DCMRotation(DirectionCosineMatrix<T, From, To>),
            EulerRotation(Euler<T>),
            QuaternionRotation(Quaternion<T>),
    }
}