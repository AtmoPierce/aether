use aether_core::math::Vector;
use aether_core::reference_frame::ReferenceFrame;
use aether_core::real::Real;

#[derive(Debug, Clone, Copy)]
struct Quat {
    w: f64,
    x: f64,
    y: f64,
    z: f64,
}

impl Quat {
    #[inline]
    fn identity() -> Self {
        Self {
            w: 1.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    #[inline]
    fn from_axis_angle(axis: &Vector<f64, 3>, angle_rad: f64) -> Self {
        let half = 0.5 * angle_rad;
        let (s, c) = half.sin_cos();
        let a = axis.normalize();
        Self {
            w: c,
            x: a[0] * s,
            y: a[1] * s,
            z: a[2] * s,
        }
        .normalized()
    }

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self {
            w: self.w * rhs.w - self.x * rhs.x - self.y * rhs.y - self.z * rhs.z,
            x: self.w * rhs.x + self.x * rhs.w + self.y * rhs.z - self.z * rhs.y,
            y: self.w * rhs.y - self.x * rhs.z + self.y * rhs.w + self.z * rhs.x,
            z: self.w * rhs.z + self.x * rhs.y - self.y * rhs.x + self.z * rhs.w,
        }
    }

    #[inline]
    fn normalized(self) -> Self {
        let n2 = self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z;
        if n2 > 0.0 {
            let inv = 1.0 / n2.sqrt();
            Self {
                w: self.w * inv,
                x: self.x * inv,
                y: self.y * inv,
                z: self.z * inv,
            }
        } else {
            Self::identity()
        }
    }

    #[inline]
    fn rotate_vec3(self, v: &Vector<f64, 3>) -> Vector<f64, 3> {
        // v' = v + 2*q.xyz×(q.xyz×v + w*v)
        // More stable form:
        // t = 2 * cross(q.xyz, v)
        // v' = v + w*t + cross(q.xyz, t)
        let qv = Vector::new([self.x, self.y, self.z]);
        let t = qv.cross(v) * 2.0;
        *v + t * self.w + qv.cross(&t)
    }
}

// 1) Define a phantom camera frame
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Camera<T: Real>(core::marker::PhantomData<T>);
impl<T: Real> ReferenceFrame for Camera<T> {}

#[derive(Debug, Clone)]
pub struct OrbitCamera {
    pub target: Vector<f64, 3>, // look-at point
    pub distance: f64,          // meters
    pub yaw: f64,               // rad
    pub pitch: f64,             // rad
    orientation: Quat,
    pub min_dist: f64,
    pub max_dist: f64,
    pub orbit_sens: f64, // rad per pixel
    pub pan_sens: f64,   // meters per pixel (scaled with distance)
    pub zoom_sens: f64,  // multiplicative per wheel notch (~120)
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self {
            target: Vector::new([0.0, 0.0, 0.0]),
            distance: 0.0,
            yaw: 0.0,
            pitch: 0.0,
            orientation: Quat::identity(),
            min_dist: 0.0,
            max_dist: 2.0e12,
            orbit_sens: 0.010,
            pan_sens: 1.25,
            zoom_sens: 1.10,
        }
    }
}

impl OrbitCamera {
    /// Sets yaw/pitch and rebuilds the internal quaternion orientation.
    ///
    /// Yaw is applied about world +Z. Pitch is applied about the camera's right
    /// axis after yaw, which matches the previous Z-up spherical mapping.
    pub fn set_angles(&mut self, yaw: f64, pitch: f64) {
        self.yaw = yaw;
        self.pitch = pitch;

        let world_up = Vector::new([0.0, 0.0, 1.0]);
        let q_yaw = Quat::from_axis_angle(&world_up, yaw);

        // Compute the post-yaw forward and right axis.
        let base_forward = Vector::new([1.0, 0.0, 0.0]);
        let fwd = q_yaw.rotate_vec3(&base_forward).normalize();
        let mut right = world_up.cross(&fwd);
        let n2 = right[0] * right[0] + right[1] * right[1] + right[2] * right[2];
        if n2 < 1.0e-16 {
            right = Vector::new([1.0, 0.0, 0.0]).cross(&fwd);
        }
        right = right.normalize();

        let q_pitch = Quat::from_axis_angle(&right, pitch);
        self.orientation = q_pitch.mul(q_yaw).normalized();
    }

    /// Applies an orbit update using mouse deltas in pixels.
    ///
    /// This updates the internal quaternion orientation (pole-safe) and
    /// refreshes `yaw`/`pitch` from the resulting forward vector.
    pub fn apply_orbit_pixels(&mut self, dx_pixels: f64, dy_pixels: f64) {
        let dyaw = -dx_pixels * self.orbit_sens;
        let dpitch = -dy_pixels * self.orbit_sens;
        self.apply_orbit_angles(dyaw, dpitch);
    }

    pub fn apply_orbit_angles(&mut self, dyaw: f64, dpitch: f64) {
        let world_up = Vector::new([0.0, 0.0, 1.0]);
        let base_forward = Vector::new([1.0, 0.0, 0.0]);

        // 1) yaw about world up
        let q_yaw = Quat::from_axis_angle(&world_up, dyaw);
        let tmp = q_yaw.mul(self.orientation).normalized();

        // 2) pitch about camera right after yaw
        let fwd_tmp = tmp.rotate_vec3(&base_forward).normalize();
        let mut right = world_up.cross(&fwd_tmp);
        let n2 = right[0] * right[0] + right[1] * right[1] + right[2] * right[2];
        if n2 < 1.0e-16 {
            right = Vector::new([1.0, 0.0, 0.0]).cross(&fwd_tmp);
        }
        right = right.normalize();
        let q_pitch = Quat::from_axis_angle(&right, dpitch);

        self.orientation = q_pitch.mul(tmp).normalized();

        // Best-effort yaw/pitch for telemetry/UI (does not drive orientation updates).
        let fwd = self.forward();
        self.yaw = fwd[1].atan2(fwd[0]);
        self.pitch = fwd[2].clamp(-1.0, 1.0).asin();
    }

    /// Applies an orbit update by rotating the camera orientation about the **world** X and Y axes.
    ///
    /// This is useful for UIs that want explicit “rotate around X/Y” controls.
    /// Note: the camera is still Z-up *in the sense that most callers use only `position()` +
    /// a fixed world-up in their look-at; roll is not represented by `right()`/`cam_up()`.
    pub fn apply_orbit_world_angles(&mut self, droll_x: f64, dpitch_y: f64) {
        let world_x = Vector::new([1.0, 0.0, 0.0]);
        let world_y = Vector::new([0.0, 1.0, 0.0]);

        let qx = Quat::from_axis_angle(&world_x, droll_x);
        let qy = Quat::from_axis_angle(&world_y, dpitch_y);

        // Apply as world-frame increments.
        self.orientation = qx.mul(qy).mul(self.orientation).normalized();

        // Refresh best-effort yaw/pitch telemetry.
        let fwd = self.forward();
        self.yaw = fwd[1].atan2(fwd[0]);
        self.pitch = fwd[2].clamp(-1.0, 1.0).asin();
    }

    /// Applies a world-axis orbit update using mouse deltas in pixels.
    ///
    /// Vertical drag rotates about world +X, horizontal drag rotates about world +Y.
    pub fn apply_orbit_world_pixels(&mut self, dx_pixels: f64, dy_pixels: f64) {
        let d_x = -dy_pixels * self.orbit_sens;
        let d_y = -dx_pixels * self.orbit_sens;
        self.apply_orbit_world_angles(d_x, d_y);
    }

    // Z-up spherical mapping: forward unit vector
    #[inline]
    pub fn forward(&self) -> Vector<f64, 3> {
        let base_forward = Vector::new([1.0, 0.0, 0.0]);
        self.orientation.rotate_vec3(&base_forward).normalize()
    }

    #[inline]
    pub fn right(&self) -> Vector<f64, 3> {
        // right = world_up × forward (Z-up locked), with robust fallback near poles.
        let fwd = self.forward();
        let world_up = Vector::new([0.0, 0.0, 1.0]);
        let mut right = world_up.cross(&fwd);

        let n2 = right[0] * right[0] + right[1] * right[1] + right[2] * right[2];
        if n2 < 1.0e-16 {
            right = Vector::new([1.0, 0.0, 0.0]).cross(&fwd);
        }

        right.normalize()
    }

    #[inline]
    pub fn cam_up(&self) -> Vector<f64, 3> {
        // up' = forward × right, locked to world-up hemisphere (+Z)
        let fwd = self.forward();
        let right = self.right();
        let up = fwd.cross(&right);
        let n2 = up[0] * up[0] + up[1] * up[1] + up[2] * up[2];
        if n2 < 1.0e-16 {
            Vector::new([0.0, 0.0, 1.0])
        } else {
            let up_n = up.normalize();
            if up_n[2] < 0.0 {
                up_n * -1.0
            } else {
                up_n
            }
        }
    }

    #[inline]
    pub fn position(&self) -> Vector<f64, 3> {
        // pos = target - fwd * distance
        self.target - self.forward() * self.distance
    }
}
