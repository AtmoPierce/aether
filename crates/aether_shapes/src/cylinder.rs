use crate::attributes::Solid;
use aether_core::{math::Vector, reference_frame::ReferenceFrame};
use aether_core::real::Real;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderLocal;
impl ReferenceFrame for CylinderLocal {}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylinder<F: Real> {
    pub r: F,
    pub h: F,
}
impl<F: Real> Solid<F> for Cylinder<F> {
    type Local = CylinderLocal;
    fn volume(&self) -> F {
        F::PI * self.r * self.r * self.h
    }
    fn inertia_principal_cm(&self, m: F) -> Vector<F, 3> {
        let half = F::from_f64(1.0) / (F::ONE + F::ONE);
        let twelve = F::ONE
            + F::ONE
            + F::ONE
            + F::ONE
            + F::ONE
            + F::ONE
            + F::ONE
            + F::ONE
            + F::ONE
            + F::ONE
            + F::ONE
            + F::ONE;
        let r2 = self.r * self.r;
        let h2 = self.h * self.h;
        let iz = m * r2 * half;
        let ixy = m * (F::ONE + F::ONE + F::ONE) * r2 + m * h2; // simplify with constants you prefer
        Vector::new([ixy / twelve, ixy / twelve, iz])
    }
}
