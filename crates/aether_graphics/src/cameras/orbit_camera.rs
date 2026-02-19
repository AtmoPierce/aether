use aether_core::attitude::Quaternion;
use aether_core::math::Vector;
use aether_core::reference_frame::ReferenceFrame;
use aether_core::real::Real;

// 1) Define a phantom camera frame
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Camera<T: Real>(core::marker::PhantomData<T>);
impl<T: Real> ReferenceFrame for Camera<T> {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpenGL<T: Real>(core::marker::PhantomData<T>);
impl<T: Real> ReferenceFrame for OpenGL<T> {}

type CamQuat<T> = Quaternion<T, OpenGL<T>, Camera<T>>;

#[derive(Debug, Clone)]
pub struct OrbitCamera<T: Real = f64> {
    pub target: Vector<T, 3>, // look-at point
    pub distance: T,          // meters
    pub yaw: T,               // rad
    pub pitch: T,             // rad
    orientation: CamQuat<T>,
    pub min_dist: T,
    pub max_dist: T,
    pub orbit_sens: T, // rad per pixel
    pub pan_sens: T,   // meters per pixel (scaled with distance)
    pub zoom_sens: T,  // multiplicative per wheel notch (~120)
}

impl<T: Real> Default for OrbitCamera<T> {
    fn default() -> Self {
        Self {
            target: Vector::new([T::ZERO, T::ZERO, T::ZERO]),
            distance: T::ZERO,
            yaw: T::ZERO,
            pitch: T::ZERO,
            orientation: CamQuat::<T>::identity(),
            min_dist: T::ZERO,
            max_dist: T::from_f64(2.0e12),
            orbit_sens: T::from_f64(0.010),
            pan_sens: T::from_f64(1.25),
            zoom_sens: T::from_f64(1.10),
        }
    }
}

impl<T: Real> OrbitCamera<T> {
    /// Sets yaw/pitch and rebuilds the internal quaternion orientation.
    ///
    /// Yaw is applied about world +Z. Pitch is applied about the camera's right
    /// axis after yaw, which matches the previous Z-up spherical mapping.
    pub fn set_angles(&mut self, yaw: T, pitch: T) {
        self.yaw = yaw;
        self.pitch = pitch;

        let world_up = Vector::new([T::ZERO, T::ZERO, T::ONE]);
        let q_yaw = CamQuat::<T>::from_axis_angle(world_up, yaw);

        // Compute the post-yaw forward and right axis.
        let base_forward = Vector::new([T::ONE, T::ZERO, T::ZERO]);
        let fwd = q_yaw.rotate_vec(&base_forward).normalize();
        let mut right = world_up.cross(&fwd);
        let n2 = right[0] * right[0] + right[1] * right[1] + right[2] * right[2];
        if n2 < T::from_f64(1.0e-16) {
            right = Vector::new([T::ONE, T::ZERO, T::ZERO]).cross(&fwd);
        }
        right = right.normalize();

        let q_pitch = CamQuat::<T>::from_axis_angle(right, pitch);
        self.orientation = q_pitch.mul_hamilton_same(&q_yaw).normalized();
    }

    /// Applies an orbit update using mouse deltas in pixels.
    ///
    /// This updates the internal quaternion orientation (pole-safe) and
    /// refreshes `yaw`/`pitch` from the resulting forward vector.
    pub fn apply_orbit_pixels(&mut self, dx_pixels: T, dy_pixels: T) {
        let dyaw = -dx_pixels * self.orbit_sens;
        let dpitch = -dy_pixels * self.orbit_sens;
        self.apply_orbit_angles(dyaw, dpitch);
    }

    pub fn apply_orbit_angles(&mut self, dyaw: T, dpitch: T) {
        let world_up = Vector::new([T::ZERO, T::ZERO, T::ONE]);
        let base_forward = Vector::new([T::ONE, T::ZERO, T::ZERO]);

        // 1) yaw about world up
        let q_yaw = CamQuat::<T>::from_axis_angle(world_up, dyaw);
        let tmp = q_yaw.mul_hamilton_same(&self.orientation).normalized();

        // 2) pitch about camera right after yaw
        let fwd_tmp = tmp.rotate_vec(&base_forward).normalize();
        let mut right = world_up.cross(&fwd_tmp);
        let n2 = right[0] * right[0] + right[1] * right[1] + right[2] * right[2];
        if n2 < T::from_f64(1.0e-16) {
            right = Vector::new([T::ONE, T::ZERO, T::ZERO]).cross(&fwd_tmp);
        }
        right = right.normalize();
        let q_pitch = CamQuat::<T>::from_axis_angle(right, dpitch);

        self.orientation = q_pitch.mul_hamilton_same(&tmp).normalized();

        // Best-effort yaw/pitch for telemetry/UI (does not drive orientation updates).
        let fwd = self.forward();
        self.yaw = fwd[1].atan2(fwd[0]);
        self.pitch = fwd[2].max(-T::ONE).min(T::ONE).asin();
    }

    /// Applies an orbit update by rotating the camera orientation about the **world** X and Y axes.
    ///
    /// This is useful for UIs that want explicit “rotate around X/Y” controls.
    /// Note: the camera is still Z-up *in the sense that most callers use only `position()` +
    /// a fixed world-up in their look-at; roll is not represented by `right()`/`cam_up()`.
    pub fn apply_orbit_world_angles(&mut self, droll_x: T, dpitch_y: T) {
        let world_x = Vector::new([T::ONE, T::ZERO, T::ZERO]);
        let world_y = Vector::new([T::ZERO, T::ONE, T::ZERO]);

        let qx = CamQuat::<T>::from_axis_angle(world_x, droll_x);
        let qy = CamQuat::<T>::from_axis_angle(world_y, dpitch_y);

        // Apply as world-frame increments.
        let qxy = qx.mul_hamilton_same(&qy).normalized();
        self.orientation = qxy.mul_hamilton_same(&self.orientation).normalized();

        // Refresh best-effort yaw/pitch telemetry.
        let fwd = self.forward();
        self.yaw = fwd[1].atan2(fwd[0]);
        self.pitch = fwd[2].max(-T::ONE).min(T::ONE).asin();
    }

    /// Applies a world-axis orbit update using mouse deltas in pixels.
    ///
    /// Vertical drag rotates about world +X, horizontal drag rotates about world +Y.
    pub fn apply_orbit_world_pixels(&mut self, dx_pixels: T, dy_pixels: T) {
        let d_x = -dy_pixels * self.orbit_sens;
        let d_y = -dx_pixels * self.orbit_sens;
        self.apply_orbit_world_angles(d_x, d_y);
    }

    // Z-up spherical mapping: forward unit vector
    #[inline]
    pub fn forward(&self) -> Vector<T, 3> {
        let base_forward = Vector::new([T::ONE, T::ZERO, T::ZERO]);
        self.orientation.rotate_vec(&base_forward).normalize()
    }

    #[inline]
    pub fn right(&self) -> Vector<T, 3> {
        // right = world_up × forward (Z-up locked), with robust fallback near poles.
        let fwd = self.forward();
        let world_up = Vector::new([T::ZERO, T::ZERO, T::ONE]);
        let mut right = world_up.cross(&fwd);

        let n2 = right[0] * right[0] + right[1] * right[1] + right[2] * right[2];
        if n2 < T::from_f64(1.0e-16) {
            right = Vector::new([T::ONE, T::ZERO, T::ZERO]).cross(&fwd);
        }

        right.normalize()
    }

    #[inline]
    pub fn cam_up(&self) -> Vector<T, 3> {
        // up' = forward × right, locked to world-up hemisphere (+Z)
        let fwd = self.forward();
        let right = self.right();
        let up = fwd.cross(&right);
        let n2 = up[0] * up[0] + up[1] * up[1] + up[2] * up[2];
        if n2 < T::from_f64(1.0e-16) {
            Vector::new([T::ZERO, T::ZERO, T::ONE])
        } else {
            let up_n = up.normalize();
            if up_n[2] < T::ZERO {
                up_n * (-T::ONE)
            } else {
                up_n
            }
        }
    }

    #[inline]
    pub fn position(&self) -> Vector<T, 3> {
        // pos = target - fwd * distance
        self.target - self.forward() * self.distance
    }
}
