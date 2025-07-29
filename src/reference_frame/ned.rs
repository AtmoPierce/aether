use crate::reference_frame::{FixedFrame, ReferenceFrame};
use crate::coordinate::Cartesian;
use num_traits::Float;


#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NED{
}

impl ReferenceFrame for NED {}

impl<T: Float> FixedFrame<T> for NED{}