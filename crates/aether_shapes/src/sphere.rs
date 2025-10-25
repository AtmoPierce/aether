use crate::attributes::Solid;
use aether_core::{math::Vector, reference_frame::ReferenceFrame};
use num_traits::{Float, FloatConst};

// Sphere
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereLocal;
impl ReferenceFrame for SphereLocal {}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere<F: Float> {
    pub r: F,
}
impl<F: Float + FloatConst> Solid<F> for Sphere<F> {
    type Local = SphereLocal;
    fn volume(&self) -> F {
        (F::one() + F::one() + F::one() + F::one()) / (F::one() + F::one() + F::one())
            * F::PI()
            * self.r
            * self.r
            * self.r
    }
    fn inertia_principal_cm(&self, m: F) -> Vector<F, 3> {
        let two_fifths =
            (F::one() + F::one()) / (F::one() + F::one() + F::one() + F::one() + F::one());
        let i = two_fifths * m * self.r * self.r;
        Vector::new([i, i, i])
    }
}
