use crate::reference_frame::{RotatingFrame, ReferenceFrame};
use crate::coordinate::Cartesian;
use num_traits::Float;

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct ICRF<T: Float>{
    pub epoch: T
}

impl<T: Float> ReferenceFrame for ICRF<T> {

}
