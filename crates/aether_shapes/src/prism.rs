use crate::attributes::Solid;
use aether_core::{math::Vector, reference_frame::ReferenceFrame};
use aether_core::real::Real;

// Prism
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrismLocal;
impl ReferenceFrame for PrismLocal {}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RectangularPrism<F: Real> {
    pub a: F,
    pub b: F,
    pub c: F,
}
impl<F: Real> Solid<F> for RectangularPrism<F> {
    type Local = PrismLocal;
    fn volume(&self) -> F {
        self.a * self.b * self.c
    }
    fn inertia_principal_cm(&self, m: F) -> Vector<F, 3> {
        let (a, b, c) = (self.a, self.b, self.c);
        let t12 = F::ONE + F::ONE + F::ONE + F::ONE + F::ONE + F::ONE; // 6? ignore; use constants you already have if any
        // Use your existing FromPrimitive helpers; shown directly:
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
        Vector::new([
            m * (b * b + c * c) / twelve,
            m * (a * a + c * c) / twelve,
            m * (a * a + b * b) / twelve,
        ])
    }
}
