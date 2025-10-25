use crate::coordinate::Cartesian;
use crate::reference_frame::{Body, ICRF};
use num_traits::{Float, FromPrimitive};

pub struct QuestObservation<T: Float, Body, ICRF> {
    pub body: Cartesian<T, Body>,
    pub inertial: Cartesian<T, ICRF>,
    pub weight: T,
}

// pub fn quest<Body, ICRF, T: Float>(observation: &[QuestObservation<T, Body, ICRF>]) -> Option<Cartesian<T, Body>> where T: Float + FromPrimitive {
//     if observation.len() < 2 {
//         return None;
//     }

// }
