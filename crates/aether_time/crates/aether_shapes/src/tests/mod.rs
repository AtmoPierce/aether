#[cfg(test)]
mod tests_add_sub {
    use super::*;
    use aether_core::attitude::DirectionCosineMatrix;
    use aether_core::coordinate::Cartesian;
    use aether_core::math::{Matrix, Vector};
    use aether_core::reference_frame::{Assembly, Pose, ReferenceFrame};

    use crate::attributes::{MassModel, Solid};
    use crate::cylinder::Cylinder;
    use crate::prism::RectangularPrism;
    use crate::sphere::Sphere;

    type F = f64;

    // --- helpers ---
    fn approx(a: F, b: F, eps: F) {
        assert!(
            (a - b).abs() <= eps,
            "a={} b={} |a-b|={}",
            a,
            b,
            (a - b).abs()
        );
    }
    fn approx_vec3(a: Vector<F, 3>, b: Vector<F, 3>, eps: F) {
        for i in 0..3 {
            approx(a[i], b[i], eps);
        }
    }
    fn approx_diag3(m: Matrix<F, 3, 3>, d: Vector<F, 3>, eps: F) {
        approx(m[(0, 0)], d[0], eps);
        approx(m[(1, 1)], d[1], eps);
        approx(m[(2, 2)], d[2], eps);
        // off-diagonals ~ 0 for aligned principal frames
        approx(m[(0, 1)], 0.0, eps);
        approx(m[(0, 2)], 0.0, eps);
        approx(m[(1, 0)], 0.0, eps);
        approx(m[(1, 2)], 0.0, eps);
        approx(m[(2, 0)], 0.0, eps);
        approx(m[(2, 1)], 0.0, eps);
    }

    // Identity pose (Local -> Assembly is identity; CM at origin)
    fn eye_pose<Local: ReferenceFrame>() -> Pose<F, Local> {
        Pose {
            R_local_to_assembly: DirectionCosineMatrix::from_matrix(Matrix::<F, 3, 3>::identity()),
            r_cm_in_assembly: Cartesian::<F, Assembly<F>>::default(),
        }
    }

    // Build some friendly test geometries
    fn outer_prism() -> RectangularPrism<F> {
        RectangularPrism {
            a: 2.0,
            b: 1.2,
            c: 0.8,
        }
    }
    fn inner_prism() -> RectangularPrism<F> {
        RectangularPrism {
            a: 1.0,
            b: 0.6,
            c: 0.4,
        }
    } // 0.24 m^3

    // Cylinder axis is z; choose height <= outer_prism.height
    fn outer_cyl() -> Cylinder<F> {
        Cylinder { r: 0.50, h: 0.80 }
    } // ~0.628 m^3
    fn inner_cyl() -> Cylinder<F> {
        Cylinder { r: 0.30, h: 0.70 }
    } // ~0.198 m^3

    fn outer_sphere() -> Sphere<F> {
        Sphere { r: 0.55 }
    } // ~0.696 m^3
    fn inner_sphere() -> Sphere<F> {
        Sphere { r: 0.35 }
    } // ~0.179 m^3

    // One density for all
    const RHO: F = 2700.0;
    const EPS: F = 1e-9;

    // -------------------------
    // Subtraction (existing)
    // -------------------------
    fn check_aligned_subtraction<O, I>(outer_shape: &O, inner_shape: &I)
    where
        O: Solid<F> + Copy,
        I: Solid<F> + Copy,
    {
        // Realize as MassModel with identity pose
        let outer = MassModel::new(*outer_shape, RHO, eye_pose::<O::Local>());
        let inner = MassModel::new(*inner_shape, RHO, eye_pose::<I::Local>());

        // Do the subtraction via operator
        let props = (outer - inner).expect("valid subtract");

        // Expected scalars
        let v_o = outer_shape.volume();
        let v_i = inner_shape.volume();
        assert!(
            v_o > v_i,
            "test geometry invalid: inner volume ({v_i}) >= outer volume ({v_o})"
        );
        let v_net = v_o - v_i;
        let m_o = RHO * v_o;
        let m_i = RHO * v_i;
        let m_net = m_o - m_i;

        approx(props.volume, v_net, EPS);
        approx(props.mass, m_net, EPS);
        approx(props.density, RHO, EPS); // same material → density stays RHO

        // Expected principal inertia: subtract
        let io_prin = outer_shape.inertia_principal_cm(m_o);
        let ii_prin = inner_shape.inertia_principal_cm(m_i);
        let expected_prin = io_prin - ii_prin;

        // Inertia should be diagonal with those principal values
        approx_diag3(props.inertia_cm, expected_prin, 1e-8);

        // COM stays at origin for concentric aligned bodies
        approx_vec3(props.com, Vector::new([0.0, 0.0, 0.0]), EPS);

        // Rotate back into the outer local frame and check diagonal again
        let r_local_to_assembly = *eye_pose::<O::Local>().R_local_to_assembly.as_matrix();
        let assembly_to_local = r_local_to_assembly.transpose();
        let i_local = assembly_to_local * props.inertia_cm * r_local_to_assembly;
        approx_diag3(i_local, expected_prin, 1e-8);
    }

    // -------------------------
    // Addition (new)
    // -------------------------
    fn check_aligned_addition<A, B>(a_shape: &A, b_shape: &B)
    where
        A: Solid<F> + Copy,
        B: Solid<F> + Copy,
    {
        let a = MassModel::new(*a_shape, RHO, eye_pose::<A::Local>());
        let b = MassModel::new(*b_shape, RHO, eye_pose::<B::Local>());

        // Do the addition via operator
        let props = (a + b).expect("valid add");

        // Expected scalars
        let v_a = a_shape.volume();
        let v_b = b_shape.volume();
        let v_sum = v_a + v_b;
        let m_a = RHO * v_a;
        let m_b = RHO * v_b;
        let m_sum = m_a + m_b;

        approx(props.volume, v_sum, EPS);
        approx(props.mass, m_sum, EPS);
        approx(props.density, RHO, EPS); // same material → density stays RHO

        // Expected principal inertia: add
        let ia_prin = a_shape.inertia_principal_cm(m_a);
        let ib_prin = b_shape.inertia_principal_cm(m_b);
        let expected_prin = ia_prin + ib_prin;

        // Inertia should be diagonal with those principal values
        approx_diag3(props.inertia_cm, expected_prin, 1e-8);

        // COM stays at origin for concentric aligned bodies
        approx_vec3(props.com, Vector::new([0.0, 0.0, 0.0]), EPS);

        // Rotate back into A's local frame and check diagonal again
        let r_local_to_assembly = *eye_pose::<A::Local>().R_local_to_assembly.as_matrix();
        let assembly_to_local = r_local_to_assembly.transpose();
        let i_local = assembly_to_local * props.inertia_cm * r_local_to_assembly;
        approx_diag3(i_local, expected_prin, 1e-8);
    }

    // --- Prism as outer ---
    #[test]
    fn prism_minus_prism() {
        check_aligned_subtraction(&outer_prism(), &inner_prism());
    }
    #[test]
    fn prism_minus_cylinder() {
        check_aligned_subtraction(&outer_prism(), &inner_cyl());
    }
    #[test]
    fn prism_minus_sphere() {
        check_aligned_subtraction(&outer_prism(), &inner_sphere());
    }

    // --- Cylinder as outer ---
    #[test]
    fn cylinder_minus_prism() {
        check_aligned_subtraction(&outer_cyl(), &inner_prism());
    }
    #[test]
    fn cylinder_minus_cylinder() {
        check_aligned_subtraction(&outer_cyl(), &inner_cyl());
    }
    #[test]
    fn cylinder_minus_sphere() {
        check_aligned_subtraction(&outer_cyl(), &inner_sphere());
    }

    // --- Sphere as outer ---
    #[test]
    fn sphere_minus_prism() {
        check_aligned_subtraction(&outer_sphere(), &inner_prism());
    }
    #[test]
    fn sphere_minus_cylinder() {
        check_aligned_subtraction(&outer_sphere(), &inner_cyl());
    }
    #[test]
    fn sphere_minus_sphere() {
        check_aligned_subtraction(&outer_sphere(), &inner_sphere());
    }

    // --- Addition mirrors the same 3×3 grid ---

    // Prism as A
    #[test]
    fn prism_plus_prism() {
        check_aligned_addition(&outer_prism(), &inner_prism());
    }
    #[test]
    fn prism_plus_cylinder() {
        check_aligned_addition(&outer_prism(), &inner_cyl());
    }
    #[test]
    fn prism_plus_sphere() {
        check_aligned_addition(&outer_prism(), &inner_sphere());
    }

    // Cylinder as A
    #[test]
    fn cylinder_plus_prism() {
        check_aligned_addition(&outer_cyl(), &inner_prism());
    }
    #[test]
    fn cylinder_plus_cylinder() {
        check_aligned_addition(&outer_cyl(), &inner_cyl());
    }
    #[test]
    fn cylinder_plus_sphere() {
        check_aligned_addition(&outer_cyl(), &inner_sphere());
    }

    // Sphere as A
    #[test]
    fn sphere_plus_prism() {
        check_aligned_addition(&outer_sphere(), &inner_prism());
    }
    #[test]
    fn sphere_plus_cylinder() {
        check_aligned_addition(&outer_sphere(), &inner_cyl());
    }
    #[test]
    fn sphere_plus_sphere() {
        check_aligned_addition(&outer_sphere(), &inner_sphere());
    }
}
#[cfg(test)]
mod tests_rot_offset {
    use super::*;
    use aether_core::attitude::DirectionCosineMatrix;
    use aether_core::coordinate::Cartesian;
    use aether_core::math::{Matrix, Vector};
    use aether_core::reference_frame::{Assembly, Pose, ReferenceFrame};

    use crate::attributes::{MassModel, Solid};
    use crate::cylinder::Cylinder;
    use crate::prism::RectangularPrism;
    use crate::sphere::Sphere;

    type F = f64;
    const RHO: F = 2700.0;
    const EPS: F = 1e-9;

    /* ---------------- helpers ---------------- */

    fn approx(a: F, b: F, eps: F) {
        assert!(
            (a - b).abs() <= eps,
            "a={} b={} |a-b|={}",
            a,
            b,
            (a - b).abs()
        );
    }
    fn approx_vec3(a: Vector<F, 3>, b: Vector<F, 3>, eps: F) {
        for i in 0..3 {
            approx(a[i], b[i], eps);
        }
    }
    fn approx_sym(m: Matrix<F, 3, 3>, eps: F) {
        approx(m[(0, 1)], m[(1, 0)], 10.0 * eps);
        approx(m[(0, 2)], m[(2, 0)], 10.0 * eps);
        approx(m[(1, 2)], m[(2, 1)], 10.0 * eps);
    }

    #[inline]
    fn deg(d: F) -> F {
        d.to_radians()
    }

    // Construct a ZYX Euler rotation (Rz*Ry*Rx) as a raw Matrix, then wrap as DCM<Local, Assembly>
    fn dcm_zyx<Local: ReferenceFrame>(
        roll_x: F,
        pitch_y: F,
        yaw_z: F,
    ) -> DirectionCosineMatrix<F, Local, Assembly<F>> {
        let (cx, sx) = (roll_x.cos(), roll_x.sin());
        let (cy, sy) = (pitch_y.cos(), pitch_y.sin());
        let (cz, sz) = (yaw_z.cos(), yaw_z.sin());

        // R = Rz * Ry * Rx
        let r = Matrix::<F, 3, 3>::new([
            [cz * cy, cz * sy * sx - sz * cx, cz * sy * cx + sz * sx],
            [sz * cy, sz * sy * sx + cz * cx, sz * sy * cx - cz * sx],
            [-sy, cy * sx, cy * cx],
        ]);

        DirectionCosineMatrix::from(&r)
    }

    // Build a Pose from Euler ZYX (deg) and a translation vector
    fn pose_euler_xyz<Local: ReferenceFrame>(
        roll_deg: F,
        pitch_deg: F,
        yaw_deg: F,
        r: [F; 3],
    ) -> Pose<F, Local> {
        Pose {
            R_local_to_assembly: dcm_zyx::<Local>(deg(roll_deg), deg(pitch_deg), deg(yaw_deg)),
            r_cm_in_assembly: Cartesian::<F, Assembly<F>> {
                data: Vector::new(r),
                _reference_frame: core::marker::PhantomData,
            },
        }
    }

    #[inline]
    fn parallel_axis(i_cm: Matrix<F, 3, 3>, m: F, r: Vector<F, 3>) -> Matrix<F, 3, 3> {
        let r2 = r.dot(&r);
        let eye = Matrix::<F, 3, 3>::identity();
        let rrT = Matrix::<F, 3, 3>::outer(&r, &r);
        i_cm + (eye * r2 - rrT) * m
    }

    // Expected subtraction: A − B (arbitrary R, r). Returns (mass, volume, com, I_cm)
    fn expected_sub<O: Solid<F>, I: Solid<F>>(
        a: &O,
        rho_a: F,
        Ra: Matrix<F, 3, 3>,
        ra: Vector<F, 3>,
        b: &I,
        rho_b: F,
        Rb: Matrix<F, 3, 3>,
        rb: Vector<F, 3>,
    ) -> Option<(F, F, Vector<F, 3>, Matrix<F, 3, 3>)> {
        let (va, vb) = (a.volume(), b.volume());
        if va <= 0.0 || vb <= 0.0 {
            return None;
        }
        let (ma, mb) = (rho_a * va, rho_b * vb);

        let Ia_prin = Matrix::<F, 3, 3>::diag_from_vector(&a.inertia_principal_cm(ma));
        let Ib_prin = Matrix::<F, 3, 3>::diag_from_vector(&b.inertia_principal_cm(mb));
        let Ia_cm = Ra * Ia_prin * Ra.transpose();
        let Ib_cm = Rb * Ib_prin * Rb.transpose();

        let Ia_O = parallel_axis(Ia_cm, ma, ra);
        let Ib_O = parallel_axis(Ib_cm, mb, rb);

        let mass = ma - mb;
        let volume = va - vb;
        if mass <= 0.0 || volume <= 0.0 {
            return None;
        }

        let com = (ra * ma - rb * mb) / mass;
        let I_O = Ia_O - Ib_O;
        let I_cm = I_O - parallel_axis(Matrix::<F, 3, 3>::zeros(), mass, com);

        Some((mass, volume, com, I_cm))
    }

    // Expected addition: A + B (arbitrary R, r).
    fn expected_add<A: Solid<F>, B: Solid<F>>(
        a: &A,
        rho_a: F,
        Ra: Matrix<F, 3, 3>,
        ra: Vector<F, 3>,
        b: &B,
        rho_b: F,
        Rb: Matrix<F, 3, 3>,
        rb: Vector<F, 3>,
    ) -> Option<(F, F, Vector<F, 3>, Matrix<F, 3, 3>)> {
        let (va, vb) = (a.volume(), b.volume());
        if va <= 0.0 || vb <= 0.0 {
            return None;
        }
        let (ma, mb) = (rho_a * va, rho_b * vb);

        let Ia_prin = Matrix::<F, 3, 3>::diag_from_vector(&a.inertia_principal_cm(ma));
        let Ib_prin = Matrix::<F, 3, 3>::diag_from_vector(&b.inertia_principal_cm(mb));
        let Ia_cm = Ra * Ia_prin * Ra.transpose();
        let Ib_cm = Rb * Ib_prin * Rb.transpose();

        let Ia_O = parallel_axis(Ia_cm, ma, ra);
        let Ib_O = parallel_axis(Ib_cm, mb, rb);

        let mass = ma + mb;
        let volume = va + vb;
        if mass <= 0.0 || volume <= 0.0 {
            return None;
        }

        let com = (ra * ma + rb * mb) / mass;
        let I_O = Ia_O + Ib_O;
        let I_cm = I_O - parallel_axis(Matrix::<F, 3, 3>::zeros(), mass, com);

        Some((mass, volume, com, I_cm))
    }

    /* -------------- shapes -------------- */

    fn prism() -> RectangularPrism<F> {
        RectangularPrism {
            a: 2.0,
            b: 1.2,
            c: 0.8,
        }
    }
    fn cyl() -> Cylinder<F> {
        Cylinder { r: 0.50, h: 0.80 }
    }
    fn sphere() -> Sphere<F> {
        Sphere { r: 0.55 }
    }

    /* -------------- subtraction tests -------------- */

    #[test]
    fn subtraction_with_rotation_and_offset() {
        // Outer: prism at origin, rotated 20° about Z
        let pose_o = pose_euler_xyz::<<RectangularPrism<F> as Solid<F>>::Local>(
            0.0,
            0.0,
            20.0,
            [0.0, 0.0, 0.0],
        );
        // Inner: cylinder rotated 35° about X and 10° about Z, offset in XY
        let pose_i = pose_euler_xyz::<<Cylinder<F> as Solid<F>>::Local>(
            35.0,
            0.0,
            10.0,
            [0.12, -0.06, 0.02],
        );

        let outer = MassModel::new(prism(), RHO, pose_o);
        let inner = MassModel::new(cyl(), RHO, pose_i);

        let got = (outer - inner).expect("valid subtract");

        // Expected via hand composition
        let Ro = *pose_o.R_local_to_assembly.as_matrix();
        let Ri = *pose_i.R_local_to_assembly.as_matrix();
        let ro = pose_o.r_cm_in_assembly.data;
        let ri = pose_i.r_cm_in_assembly.data;

        let (mass, volume, com, I_cm) =
            expected_sub(&prism(), RHO, Ro, ro, &cyl(), RHO, Ri, ri).expect("expected");

        approx(got.mass, mass, 1e-8);
        approx(got.volume, volume, 1e-12);
        approx_vec3(got.com, com, 1e-12);
        approx_sym(got.inertia_cm, 1e-10);
        // Compare all tensor elements (looser tol due to trig)
        for i in 0..3 {
            for j in 0..3 {
                approx(got.inertia_cm[(i, j)], I_cm[(i, j)], 1e-8);
            }
        }
        // Density stays RHO for same-material difference (ρ*(Vo-Vi)/(Vo-Vi) == ρ)
        approx(got.density, RHO, 1e-12);
    }

    /* -------------- addition tests -------------- */

    #[test]
    fn addition_with_two_rotated_offset_parts() {
        // A: sphere rotated 15° about Y, offset +x,+z
        let pose_a =
            pose_euler_xyz::<<Sphere<F> as Solid<F>>::Local>(0.0, 15.0, 0.0, [0.08, 0.00, 0.05]);
        // B: prism rotated 10° about X and 25° about Z, offset -y
        let pose_b = pose_euler_xyz::<<RectangularPrism<F> as Solid<F>>::Local>(
            10.0,
            0.0,
            25.0,
            [0.00, -0.07, 0.00],
        );

        let a = MassModel::new(sphere(), RHO, pose_a);
        let b = MassModel::new(prism(), RHO, pose_b);

        let got = (a + b).expect("valid add");

        // Expected via hand composition
        let Ra = *pose_a.R_local_to_assembly.as_matrix();
        let Rb = *pose_b.R_local_to_assembly.as_matrix();
        let ra = pose_a.r_cm_in_assembly.data;
        let rb = pose_b.r_cm_in_assembly.data;

        let (mass, volume, com, I_cm) =
            expected_add(&sphere(), RHO, Ra, ra, &prism(), RHO, Rb, rb).expect("expected");

        approx(got.mass, mass, 1e-8);
        approx(got.volume, volume, 1e-12);
        approx_vec3(got.com, com, 1e-12);
        approx_sym(got.inertia_cm, 1e-10);
        for i in 0..3 {
            for j in 0..3 {
                approx(got.inertia_cm[(i, j)], I_cm[(i, j)], 1e-8);
            }
        }
        approx(got.density, RHO, 1e-12);
    }
}
