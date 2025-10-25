use num_traits::Float;

use super::cartesian::Cartesian;
use super::spherical::Spherical;
use crate::math::Vector;
use crate::reference_frame::ReferenceFrame;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylindrical<T: Float> {
    pub data: Vector<T, 3>, // [r, theta, z]
}
impl<T: Float> Cylindrical<T> {
    pub fn new(r: T, theta: T, z: T) -> Self {
        Self {
            data: Vector {
                data: [r, theta, z],
            },
        }
    }
    pub fn r(&self) -> T {
        self.data.data[0]
    }
    pub fn theta(&self) -> T {
        self.data.data[1]
    }
    pub fn z(&self) -> T {
        self.data.data[2]
    }
}

impl<T: Float, RF: ReferenceFrame> From<&Cartesian<T, RF>> for Cylindrical<T> {
    fn from(cart: &Cartesian<T, RF>) -> Self {
        let x = cart.x();
        let y = cart.y();
        let z = cart.z();
        let r = (x * x + y * y).sqrt();
        let theta = y.atan2(x);
        Cylindrical::new(r, theta, z)
    }
}

impl<T: Float> From<&Spherical<T>> for Cylindrical<T> {
    fn from(s: &Spherical<T>) -> Self {
        let r = s.r() * s.theta().sin();
        let z = s.r() * s.theta().cos();
        let theta = s.phi();
        return Cylindrical::new(r, theta, z);
    }
}
