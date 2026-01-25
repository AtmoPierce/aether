use num_traits::{Float, Pow};

use super::cartesian::Cartesian;
use super::cylindrical::Cylindrical;
use crate::{math::Vector, reference_frame::ReferenceFrame};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spherical<T: Float> {
    /// r: radial distance
    /// phi: azimuthal angle (radians, from +x in xy-plane)
    /// theta: inclination angle (radians (z-axis))
    pub data: Vector<T, 3>, // [r, phi, theta]
}

impl<T: Float> Spherical<T> {
    pub fn new(r: T, phi: T, theta: T) -> Self {
        Self {
            data: Vector {
                data: [r, phi, theta],
            },
        }
    }
    pub fn r(&self) -> T {
        self.data.data[0]
    }
    pub fn phi(&self) -> T {
        self.data.data[1]
    }
    pub fn theta(&self) -> T {
        self.data.data[2]
    }
}

impl<T: Float, RF: ReferenceFrame> From<&Cartesian<T, RF>> for Spherical<T> {
    fn from(cart: &Cartesian<T, RF>) -> Self {
        let x = cart.x();
        let y = cart.y();
        let z = cart.z();
        let rho = (x * x + y * y + z * z).sqrt();
        let theta = y.atan2(x); // azimuth angle in xy-plane
        let phi = if rho == T::zero() {
            T::zero()
        } else {
            (z / rho).acos()
        }; // inclination from z-axis
        Spherical::new(rho, theta, phi)
    }
}

impl<T: Float> From<&Cylindrical<T>> for Spherical<T> {
    fn from(c: &Cylindrical<T>) -> Self {
        let r_cyl = c.r();
        let theta = c.theta(); // azimuthal angle in x-y plane
        let z = c.z();

        let rho = (r_cyl * r_cyl + z * z).sqrt(); // spherical radius
        let phi = if rho != T::zero() {
            (z / rho).acos() // inclination from z-axis
        } else {
            T::zero()
        };
        Spherical::new(rho, theta, phi)
    }
}
