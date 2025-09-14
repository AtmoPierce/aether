use crate::reference_frame::{FixedFrame, ReferenceFrame};
use crate::coordinate::Cartesian;
use num_traits::Float;


#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct NED<T: Float>{
    pub origin: T
}

impl<T: Float> ReferenceFrame for NED<T>{}

impl<T: Float> FixedFrame<T> for NED<T>{}