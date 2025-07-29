use crate::reference_frame::{RotatingFrame, ReferenceFrame};
use crate::coordinate::Cartesian;
use num_traits::Float;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ICRF{
}

impl ReferenceFrame for ICRF {}
