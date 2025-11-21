use crate::attributes::Solid;
use aether_core::{math::Vector, reference_frame::ReferenceFrame};
use aether_core::real::{Real};

// Sphere
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereLocal;
impl ReferenceFrame for SphereLocal {}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere<F: Real> {
    pub r: F,
}
impl<F: Real> Solid<F> for Sphere<F> {
    type Local = SphereLocal;
    fn volume(&self) -> F {
        (F::ONE + F::ONE + F::ONE + F::ONE) / (F::ONE + F::ONE + F::ONE)
            * F::PI
            * self.r
            * self.r
            * self.r
    }
    fn inertia_principal_cm(&self, m: F) -> Vector<F, 3> {
        let two_fifths =
            (F::ONE + F::ONE) / (F::ONE + F::ONE + F::ONE + F::ONE + F::ONE);
        let i = two_fifths * m * self.r * self.r;
        Vector::new([i, i, i])
    }
}
