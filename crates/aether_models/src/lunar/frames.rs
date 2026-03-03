use crate::reference_frame::ReferenceFrame;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct LIRF;

impl ReferenceFrame for LIRF {}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct LTRF;

impl ReferenceFrame for LTRF {}
