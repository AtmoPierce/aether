pub mod cartesian;
pub mod cylindrical;
pub mod spherical;
pub use cartesian::Cartesian;
pub use cylindrical::Cylindrical;
pub use spherical::Spherical;

pub mod coordinate {
    use crate::{
        coordinate::{Cartesian, Cylindrical, Spherical},
        reference_frame::ReferenceFrame,
    };
    use num_traits::Float;

    pub enum Coordinate<T: Float, F: ReferenceFrame> {
        CartesianValue(Cartesian<T, F>),
        CylindricalValue(Cylindrical<T>),
        SphericalValue(Spherical<T>),
    }
}

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
