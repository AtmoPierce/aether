use crate::coordinate::Cartesian;
use crate::reference_frame::ReferenceFrame;
use num_traits::Float;

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct ICRF<T: Float> {
    pub epoch: T,
}

impl<T: Float> ReferenceFrame for ICRF<T> {}
