pub mod lerp;
pub mod polynomial;
pub mod spline;

pub use lerp::{inverse_lerp, lerp, remap};
pub use polynomial::PolynomialInterpolator;
pub use spline::CubicSpline;
