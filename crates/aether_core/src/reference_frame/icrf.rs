use crate::coordinate::Cartesian;
use crate::reference_frame::ReferenceFrame;
use crate::real::Real;
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct ICRF<T: Real> {
    pub epoch: T,
}

impl<T: Real> ReferenceFrame for ICRF<T> {}
