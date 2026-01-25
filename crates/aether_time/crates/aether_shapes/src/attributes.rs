use aether_core::{
    math::{Matrix, Vector},
    reference_frame::{Pose, ReferenceFrame},
};
use core::ops::{Add, Sub};
use num_traits::Float;

/// Geometry-only contract (principal moments are in the shape's Local frame).
pub trait Solid<F: Float> {
    type Local: ReferenceFrame;
    fn volume(&self) -> F;
    fn inertia_principal_cm(&self, mass: F) -> Vector<F, 3>;
}

/// A realized solid: shape + density + pose (Local → Assembly, CM position in Assembly)
pub struct MassModel<T: Float, S: Solid<T>> {
    pub solid: S,
    pub density: T,
    pub pose: Pose<T, <S as Solid<T>>::Local>,
}

impl<T: Float, S: Solid<T>> MassModel<T, S> {
    pub fn new(solid: S, density: T, pose: Pose<T, <S as Solid<T>>::Local>) -> Self {
        Self {
            solid,
            density,
            pose,
        }
    }
}

/// Resulting mass properties about the **combined CM**, expressed in the Assembly frame.
#[derive(Debug, Clone, Copy)]
pub struct MassProps<T: Float> {
    pub volume: T,
    pub density: T,
    pub mass: T,
    pub com: Vector<T, 3>,           // CM in Assembly frame
    pub inertia_cm: Matrix<T, 3, 3>, // inertia about combined CM, in Assembly frame
}

impl<T: Float, SO: Solid<T>, SI: Solid<T>> Sub<MassModel<T, SI>> for MassModel<T, SO> {
    type Output = Option<MassProps<T>>;

    fn sub(self, inner: MassModel<T, SI>) -> Self::Output {
        // Volumes & masses
        let v_o = self.solid.volume();
        let v_i = inner.solid.volume();
        if v_o <= T::zero() || v_i <= T::zero() {
            return None;
        }
        let m_o = self.density * v_o;
        let m_i = inner.density * v_i;

        // Principal inertia (Local) → rotate into Assembly (still at each CM)
        let Io_prin = Matrix::<T, 3, 3>::diag_from_vector(&self.solid.inertia_principal_cm(m_o));
        let Ii_prin = Matrix::<T, 3, 3>::diag_from_vector(&inner.solid.inertia_principal_cm(m_i));

        let Ro = *self.pose.R_local_to_assembly.as_matrix(); // Local_o → Assembly
        let Ri = *inner.pose.R_local_to_assembly.as_matrix(); // Local_i → Assembly

        let Io_cm = Ro * Io_prin * Ro.transpose();
        let Ii_cm = Ri * Ii_prin * Ri.transpose();

        // Shift each tensor from its CM to the Assembly origin
        let r_o = self.pose.r_cm_in_assembly.data; // Vector<T,3>
        let r_i = inner.pose.r_cm_in_assembly.data; // Vector<T,3>
        let Io_O = parallel_axis(Io_cm, m_o, r_o);
        let Ii_O = parallel_axis(Ii_cm, m_i, r_i);

        // Net scalars
        let mass = m_o - m_i;
        let volume = v_o - v_i;
        if mass <= T::zero() || volume <= T::zero() {
            return None;
        }
        let density = mass / volume;

        // Combined COM (Assembly)
        let com = (r_o * m_o - r_i * m_i) / mass;

        // Inertia about Assembly origin → shift back to combined CM
        let I_O = Io_O - Ii_O;
        let I_cm = I_O - parallel_axis(Matrix::<T, 3, 3>::zeros(), mass, com);

        Some(MassProps {
            volume,
            density,
            mass,
            com,
            inertia_cm: I_cm,
        })
    }
}

impl<T: Float, SO: Solid<T>, SI: Solid<T>> Add<MassModel<T, SI>> for MassModel<T, SO> {
    type Output = Option<MassProps<T>>;
    fn add(self, other: MassModel<T, SI>) -> Self::Output {
        // Volumes & masses
        let v_a = self.solid.volume();
        let v_b = other.solid.volume();
        if v_a <= T::zero() || v_b <= T::zero() {
            return None;
        }
        let m_a = self.density * v_a;
        let m_b = other.density * v_b;

        // Principal inertia (Local) → rotate into Assembly (still at each CM)
        let Ia_prin = Matrix::<T, 3, 3>::diag_from_vector(&self.solid.inertia_principal_cm(m_a));
        let Ib_prin = Matrix::<T, 3, 3>::diag_from_vector(&other.solid.inertia_principal_cm(m_b));

        let Ra = *self.pose.R_local_to_assembly.as_matrix(); // Local_a → Assembly
        let Rb = *other.pose.R_local_to_assembly.as_matrix(); // Local_b → Assembly

        let Ia_cm = Ra * Ia_prin * Ra.transpose();
        let Ib_cm = Rb * Ib_prin * Rb.transpose();

        // Shift each tensor from its CM to the Assembly origin
        let r_a = self.pose.r_cm_in_assembly.data; // Vector<T,3>
        let r_b = other.pose.r_cm_in_assembly.data; // Vector<T,3>
        let Ia_O = parallel_axis(Ia_cm, m_a, r_a);
        let Ib_O = parallel_axis(Ib_cm, m_b, r_b);

        // Net scalars
        let mass = m_a + m_b;
        let volume = v_a + v_b;
        if mass <= T::zero() || volume <= T::zero() {
            return None;
        }
        let density = mass / volume;

        // Combined COM (Assembly)
        let com = (r_a * m_a + r_b * m_b) / mass;

        // Inertia about Assembly origin → shift back to combined CM
        let I_O = Ia_O + Ib_O;
        let I_cm = I_O - parallel_axis(Matrix::<T, 3, 3>::zeros(), mass, com);

        Some(MassProps {
            volume,
            density,
            mass,
            com,
            inertia_cm: I_cm,
        })
    }
}

#[inline]
fn parallel_axis<T: Float>(i_cm: Matrix<T, 3, 3>, m: T, r: Vector<T, 3>) -> Matrix<T, 3, 3> {
    // I_P = I_CM + m (‖r‖² I − r rᵀ)
    let r2 = r.dot(&r);
    let eye = Matrix::<T, 3, 3>::identity();
    let rrT = Matrix::<T, 3, 3>::outer(&r, &r);
    i_cm + (eye * r2 - rrT) * m
}

/* ----------------------------------------------------------------------------
Optional: keep the "matrix-first" ergonomic API for direct use in tests
---------------------------------------------------------------------------- */

pub trait SolidOps<F: Float>: Solid<F> {
    fn subtract_with<I: Solid<F>>(
        &self,
        inner: &I,
        rho_outer: F,
        rho_inner: F,
        R_o: Matrix<F, 3, 3>,
        r_o: Vector<F, 3>,
        R_i: Matrix<F, 3, 3>,
        r_i: Vector<F, 3>,
    ) -> Option<MassProps<F>>;
    fn add_with<I: Solid<F>>(
        &self,
        other: &I,
        rho_self: F,
        rho_other: F,
        R_self: Matrix<F, 3, 3>,
        r_self: Vector<F, 3>,
        R_other: Matrix<F, 3, 3>,
        r_other: Vector<F, 3>,
    ) -> Option<MassProps<F>>;
}

impl<F: Float, S: Solid<F>> SolidOps<F> for S {
    fn subtract_with<I: Solid<F>>(
        &self,
        inner: &I,
        rho_outer: F,
        rho_inner: F,
        R_o: Matrix<F, 3, 3>,
        r_o: Vector<F, 3>,
        R_i: Matrix<F, 3, 3>,
        r_i: Vector<F, 3>,
    ) -> Option<MassProps<F>> {
        let v_o = self.volume();
        let v_i = inner.volume();
        if v_o <= F::zero() || v_i <= F::zero() {
            return None;
        }
        let m_o = rho_outer * v_o;
        let m_i = rho_inner * v_i;

        let Io_prin = Matrix::<F, 3, 3>::diag_from_vector(&self.inertia_principal_cm(m_o));
        let Ii_prin = Matrix::<F, 3, 3>::diag_from_vector(&inner.inertia_principal_cm(m_i));

        let Io_cm = R_o * Io_prin * R_o.transpose();
        let Ii_cm = R_i * Ii_prin * R_i.transpose();

        let Io_O = parallel_axis(Io_cm, m_o, r_o);
        let Ii_O = parallel_axis(Ii_cm, m_i, r_i);

        let mass = m_o - m_i;
        let volume = v_o - v_i;
        if mass <= F::zero() || volume <= F::zero() {
            return None;
        }
        let density = mass / volume;

        let com = (r_o * m_o - r_i * m_i) / mass;

        let I_O = Io_O - Ii_O;
        let I_cm = I_O - parallel_axis(Matrix::<F, 3, 3>::zeros(), mass, com);

        Some(MassProps {
            volume,
            density,
            mass,
            com,
            inertia_cm: I_cm,
        })
    }

    fn add_with<I: Solid<F>>(
        &self,
        other: &I,
        rho_self: F,
        rho_other: F,
        R_self: Matrix<F, 3, 3>,
        r_self: Vector<F, 3>,
        R_other: Matrix<F, 3, 3>,
        r_other: Vector<F, 3>,
    ) -> Option<MassProps<F>> {
        let v_a = self.volume();
        let v_b = other.volume();
        if v_a <= F::zero() || v_b <= F::zero() {
            return None;
        }
        let m_a = rho_self * v_a;
        let m_b = rho_other * v_b;

        let Ia_prin = Matrix::<F, 3, 3>::diag_from_vector(&self.inertia_principal_cm(m_a));
        let Ib_prin = Matrix::<F, 3, 3>::diag_from_vector(&other.inertia_principal_cm(m_b));

        let Ia_cm = R_self * Ia_prin * R_self.transpose();
        let Ib_cm = R_other * Ib_prin * R_other.transpose();

        let Ia_O = parallel_axis(Ia_cm, m_a, r_self);
        let Ib_O = parallel_axis(Ib_cm, m_b, r_other);

        let mass = m_a + m_b;
        let volume = v_a + v_b;
        if mass <= F::zero() || volume <= F::zero() {
            return None;
        }
        let density = mass / volume;

        let com = (r_self * m_a + r_other * m_b) / mass;

        let I_O = Ia_O + Ib_O;
        let I_cm = I_O - parallel_axis(Matrix::<F, 3, 3>::zeros(), mass, com);

        Some(MassProps {
            volume,
            density,
            mass,
            com,
            inertia_cm: I_cm,
        })
    }
}
