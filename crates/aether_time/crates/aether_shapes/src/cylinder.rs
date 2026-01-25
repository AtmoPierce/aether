use crate::attributes::Solid;
use aether_core::{math::Vector, reference_frame::ReferenceFrame};
use num_traits::{Float, FloatConst};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderLocal;
impl ReferenceFrame for CylinderLocal {}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylinder<F: Float> {
    pub r: F,
    pub h: F,
}
impl<F: Float + FloatConst> Solid<F> for Cylinder<F> {
    type Local = CylinderLocal;
    fn volume(&self) -> F {
        F::PI() * self.r * self.r * self.h
    }
    fn inertia_principal_cm(&self, m: F) -> Vector<F, 3> {
        let half = (F::one() + F::one()).recip();
        let twelve = F::one()
            + F::one()
            + F::one()
            + F::one()
            + F::one()
            + F::one()
            + F::one()
            + F::one()
            + F::one()
            + F::one()
            + F::one()
            + F::one();
        let r2 = self.r * self.r;
        let h2 = self.h * self.h;
        let iz = m * r2 * half;
        let ixy = m * (F::one() + F::one() + F::one()) * r2 + m * h2; // simplify with constants you prefer
        Vector::new([ixy / twelve, ixy / twelve, iz])
    }
}
