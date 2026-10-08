pub mod dcm;
pub mod determination;
pub mod euler;
pub mod quaternion;
pub mod rotation;
pub use dcm::DirectionCosineMatrix;
pub use euler::Euler;
pub use quaternion::Quaternion;
pub use rotation::Rotation;
#[cfg(test)]
#[path = "tests/mod.rs"]
pub mod tests;

pub mod attitude {
    use super::dcm::DirectionCosineMatrix;
    use super::euler::Euler;
    use super::quaternion::Quaternion;
    use crate::real::Real;
    use crate::reference_frame::ReferenceFrame;
    pub enum Rotation<T, From: ReferenceFrame, To: ReferenceFrame>
    where
        T: Real,
    {
        DCMRotation(DirectionCosineMatrix<T, From, To>),
        EulerRotation(Euler<T, From, To>),
        QuaternionRotation(Quaternion<T, From, To>),
    }
}
