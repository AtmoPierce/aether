//! The `Assembly` frame moves with the object (e.g., spacecraft), and its orientation
//! is defined externally (e.g., from a quaternion in the simulation state).
//! It is treated as a `FixedFrame` in this context, meaning the transformation to
//! inertial or world frames must be provided by the surrounding system.
use crate::coordinate::Cartesian;
use crate::reference_frame::{FixedFrame, ReferenceFrame};
use crate::real::Real;
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Assembly<T: Real> {
    /// Vector from origin to center of mass in Assembly coordinates
    pub mass_center: Cartesian<T, Assembly<T>>,
    /// Vector from origin to sensor/actuator location (e.g., IMU offset)
    pub lever_arm: Cartesian<T, Assembly<T>>,
}

impl<T: Real> Assembly<T> {
    pub fn new(
        mass_center: Cartesian<T, Assembly<T>>,
        lever_arm: Cartesian<T, Assembly<T>>,
    ) -> Self {
        Self {
            mass_center: mass_center,
            lever_arm: lever_arm,
        }
    }
}

impl<T: Real + Default> Default for Assembly<T> {
    fn default() -> Self {
        Self {
            mass_center: Cartesian::default(),
            lever_arm: Cartesian::default(),
        }
    }
}

// Implement ReferenceFrame for Assembly<T>
impl<T: Real> ReferenceFrame for Assembly<T> {}

// Implement FixedFrame for Assembly<T>
impl<T: Real> FixedFrame<T> for Assembly<T> {}

#[cfg(feature = "serde")]
use serde::{Deserialize, Deserializer, Serialize, Serializer};

#[cfg(feature = "serde")]
impl<T> Serialize for Assembly<T>
where
    T: Real + Serialize,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        (&self.mass_center, &self.lever_arm).serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T> Deserialize<'de> for Assembly<T>
where
    T: Real + Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let (mass_center, lever_arm): (Cartesian<T, Assembly<T>>, Cartesian<T, Assembly<T>>) =
            Deserialize::deserialize(deserializer)?;
        Ok(Assembly {
            mass_center,
            lever_arm,
        })
    }
}
