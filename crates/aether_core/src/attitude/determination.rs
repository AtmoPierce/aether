use crate::coordinate::Cartesian;
use crate::reference_frame::{Body, ICRF};
use crate::real::Real;
pub struct QuestObservation<T: Real, Body, ICRF> {
    pub body: Cartesian<T, Body>,
    pub inertial: Cartesian<T, ICRF>,
    pub weight: T,
}

// pub fn quest<Body, ICRF, T: Real>(observation: &[QuestObservation<T, Body, ICRF>]) -> Option<Cartesian<T, Body>> where T: Real + FromPrimitive {
//     if observation.len() < 2 {
//         return None;
//     }

// }
