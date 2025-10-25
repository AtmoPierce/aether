use aether_core::math::Vector;
use aether_core::reference_frame::ReferenceFrame;

// 1) Define a phantom camera frame
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Camera<T: num_traits::Float>(core::marker::PhantomData<T>);
impl<T: num_traits::Float> ReferenceFrame for Camera<T> {}

#[derive(Debug, Clone)]
pub struct OrbitCamera {
    pub target: Vector<f64, 3>, // look-at point
    pub distance: f64,          // meters
    pub yaw: f64,               // rad
    pub pitch: f64,             // rad
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
            distance: 5_000.0,
            yaw: 0.0,
            pitch: 0.3,
            min_dist: 0.1,
            max_dist: 2.0e12,
            orbit_sens: 0.010,
            pan_sens: 1.25,
            zoom_sens: 1.10,
        }
    }
}

impl OrbitCamera {
    // Z-up spherical mapping: forward unit vector
    #[inline]
    pub fn forward(&self) -> Vector<f64, 3> {
        let cp = self.pitch.cos();
        Vector::new([self.yaw.cos() * cp, self.yaw.sin() * cp, self.pitch.sin()]).normalize()
    }

    #[inline]
    pub fn right(&self) -> Vector<f64, 3> {
        // right = forward × up (Z-up)
        self.forward()
            .cross(&Vector::new([0.0, 0.0, 1.0]))
            .normalize()
    }

    #[inline]
    pub fn cam_up(&self) -> Vector<f64, 3> {
        // up' = right × forward
        self.right().cross(&self.forward()).normalize()
    }

    #[inline]
    pub fn position(&self) -> Vector<f64, 3> {
        // pos = target - fwd * distance
        self.target - self.forward() * self.distance
    }
}
