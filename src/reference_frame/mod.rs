#![allow(non_snake_case)]

use core::marker::PhantomData;

pub mod transforms;
mod traits;
mod icrf;
// mod gcrf;
mod itrf;
mod ned;
mod body;
// mod unknown;

pub use traits::{ReferenceFrame, FixedFrame, RotatingFrame};
pub use icrf::ICRF;
pub use itrf::ITRF;
pub use ned::NED;
pub use body::Body;
// pub use unkown::Unknown;

/// Enum for checking frame categories at runtime if needed
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameKind {
    Fixed,
    Rotating,
    Unknown,
}