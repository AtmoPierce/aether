#![allow(non_snake_case)]

use core::marker::PhantomData;
mod asm;
mod icrf;
mod traits;
pub mod transforms;
// mod gcrf;
mod body;
mod itrf;
mod ned;
mod pose;
mod unknown;

pub use asm::Assembly;
pub use body::Body;
pub use icrf::ICRF;
pub use itrf::ITRF;
pub use ned::NED;
pub use pose::Pose;
pub use traits::{FixedFrame, ReferenceFrame, RotatingFrame};
pub use unknown::Unknown;

/// Enum for checking frame categories at runtime if needed
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameKind {
    Fixed,
    Rotating,
    Unknown,
}
