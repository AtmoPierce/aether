//! Body reference frame
//!
//! In this context, the `Body` frame is assumed to be a fixed frame relative to the object it is attached to.
//! It moves with the object (e.g., a spacecraft or aircraft), but it does not rotate independently.
//! Instead, its orientation is typically defined by the simulation state, such as quaternions or DCMs,
//! and evolves externally as the object moves through space.

use crate::reference_frame::{FixedFrame, ReferenceFrame};
use crate::coordinate::Cartesian;
use num_traits::Float;


#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Body<T: Float>{
    /// Vector from origin to center of mass in body coordinates
    pub mass_center: Cartesian<T, Body<T>>,
    /// Vector from origin to sensor/actuator location (e.g., IMU offset)
    pub lever_arm: Cartesian<T, Body<T>>,
}

impl<T:Float> Body<T>{
    pub fn new(mass_center: Cartesian<T, Body<T>>, lever_arm: Cartesian<T, Body<T>>)->Self{
        Self { mass_center: mass_center, lever_arm: lever_arm}
    }
}

impl<T: Float + Default> Default for Body<T> {
    fn default() -> Self {
        Self { mass_center: Cartesian::default(), lever_arm: Cartesian::default() }
    }
}

// Implement ReferenceFrame for Body<T>
impl<T: Float> ReferenceFrame for Body<T> {}

// Implement FixedFrame for Body<T>
impl<T: Float> FixedFrame<T> for Body<T> {}