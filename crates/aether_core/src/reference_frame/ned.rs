use crate::coordinate::Cartesian;
use crate::reference_frame::{FixedFrame, ReferenceFrame};
use crate::real::Real;
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct NED<T: Real> {
    pub origin: T,
}

impl<T: Real> ReferenceFrame for NED<T> {}

impl<T: Real> FixedFrame<T> for NED<T> {}
