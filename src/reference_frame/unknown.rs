use crate::reference_frame::{RotatingFrame, ReferenceFrame};

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Default, Clone, Copy)]
pub struct Unknown {

}

impl ReferenceFrame for Unknown {

}
