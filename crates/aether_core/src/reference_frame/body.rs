//! The `Body` frame moves with the object (e.g., spacecraft), and its orientation
//! is defined externally (e.g., from a quaternion in the simulation state).
//! It is treated as a `FixedFrame` in this context, meaning the transformation to
//! inertial or world frames must be provided by the surrounding system.
use crate::coordinate::Cartesian;
use crate::reference_frame::{FixedFrame, ReferenceFrame};
use num_traits::Float;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Body<T: Float> {
    /// Vector from origin to center of mass in body coordinates
    pub mass_center: Cartesian<T, Body<T>>,
    /// Vector from origin to sensor/actuator location (e.g., IMU offset)
    pub lever_arm: Cartesian<T, Body<T>>,
}

impl<T: Float> Body<T> {
    pub fn new(mass_center: Cartesian<T, Body<T>>, lever_arm: Cartesian<T, Body<T>>) -> Self {
        Self {
            mass_center: mass_center,
            lever_arm: lever_arm,
        }
    }
}

impl<T: Float + Default> Default for Body<T> {
    fn default() -> Self {
        Self {
            mass_center: Cartesian::default(),
            lever_arm: Cartesian::default(),
        }
    }
}

// Implement ReferenceFrame for Body<T>
impl<T: Float> ReferenceFrame for Body<T> {}

// Implement FixedFrame for Body<T>
impl<T: Float> FixedFrame<T> for Body<T> {}

#[cfg(feature = "serde")]
use serde::{Deserialize, Deserializer, Serialize, Serializer};

#[cfg(feature = "serde")]
impl<T> Serialize for Body<T>
where
    T: Float + Serialize,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        (&self.mass_center, &self.lever_arm).serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T> Deserialize<'de> for Body<T>
where
    T: Float + Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let (mass_center, lever_arm): (Cartesian<T, Body<T>>, Cartesian<T, Body<T>>) =
            Deserialize::deserialize(deserializer)?;
        Ok(Body {
            mass_center,
            lever_arm,
        })
    }
}
