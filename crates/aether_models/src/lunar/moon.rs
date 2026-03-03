use crate::{
    attitude::DirectionCosineMatrix,
    coordinate::Cartesian,
    math::Matrix,
    reference_frame::{ICRF, ReferenceFrame},
};

use super::{
    constants,
    frames::{LIRF, LTRF},
    transforms,
};

#[derive(Clone, Debug)]
pub struct Moon {
    pub mass: f64,
    pub radius: f64,
    pub rotational_rate_rad_s: f64,
    pub gravitational_constant: f64,
    pub epoch_time: f64,
    pub icrf_to_lirf_at_epoch: DirectionCosineMatrix<f64, ICRF<f64>, LIRF>,
}

impl Default for Moon {
    fn default() -> Self {
        Self {
            mass: constants::MOON_MASS,
            radius: constants::MOON_MEAN_RADIUS,
            rotational_rate_rad_s: constants::MOON_SIDEREAL_ROTATION_RATE_RAD_S,
            gravitational_constant: constants::SELENOCENTRIC_GRAVITATIONAL_CONSTANT,
            epoch_time: 0.0,
            icrf_to_lirf_at_epoch: DirectionCosineMatrix::new(
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
            ),
        }
    }
}

impl Moon {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_epoch_orientation(
        epoch_time: f64,
        icrf_to_lirf_at_epoch: DirectionCosineMatrix<f64, ICRF<f64>, LIRF>,
    ) -> Self {
        Self {
            epoch_time,
            icrf_to_lirf_at_epoch,
            ..Self::default()
        }
    }

    pub fn gravitational_acceleration<RF: ReferenceFrame>(
        &self,
        position_from_center: Cartesian<f64, RF>,
    ) -> Cartesian<f64, RF> {
        let radius = position_from_center.norm();
        if radius <= f64::EPSILON {
            return Cartesian::zero();
        }

        let scalar = -self.gravitational_constant / (radius.powi(3));
        position_from_center * scalar
    }

    pub fn solve_gravitational_force<RF: ReferenceFrame>(
        &self,
        position_from_center: Cartesian<f64, RF>,
        mass: f64,
    ) -> Cartesian<f64, RF> {
        self.gravitational_acceleration(position_from_center) * mass
    }

    pub fn solve_gravitational_torque<RF: ReferenceFrame>(
        &self,
        position_from_center: Cartesian<f64, RF>,
        inertia_matrix: Matrix<f64, 3, 3>,
    ) -> Cartesian<f64, RF> {
        let radius = position_from_center.norm();
        if radius <= f64::EPSILON {
            return Cartesian::zero();
        }

        let value = inertia_matrix * &position_from_center;
        let local = position_from_center.cross(&value);
        let scalar = (3.0 * self.gravitational_constant) / radius.powi(5);
        local * scalar
    }

    pub fn lirf_to_ltrf(&self, time: f64) -> DirectionCosineMatrix<f64, LIRF, LTRF> {
        transforms::lirf_to_ltrf(time)
    }

    pub fn ltrf_to_lirf(&self, time: f64) -> DirectionCosineMatrix<f64, LTRF, LIRF> {
        transforms::ltrf_to_lirf(time)
    }

    pub fn icrf_to_lirf(&self) -> DirectionCosineMatrix<f64, ICRF<f64>, LIRF> {
        self.icrf_to_lirf_at_epoch
    }

    pub fn lirf_to_icrf(&self) -> DirectionCosineMatrix<f64, LIRF, ICRF<f64>> {
        self.icrf_to_lirf_at_epoch.transpose()
    }

    pub fn icrf_to_ltrf(&self, time: f64) -> DirectionCosineMatrix<f64, ICRF<f64>, LTRF> {
        transforms::icrf_to_ltrf_at_epoch(time - self.epoch_time, self.icrf_to_lirf_at_epoch)
    }

    pub fn ltrf_to_icrf(&self, time: f64) -> DirectionCosineMatrix<f64, LTRF, ICRF<f64>> {
        transforms::ltrf_to_icrf_at_epoch(time - self.epoch_time, self.icrf_to_lirf_at_epoch)
    }
}

#[cfg(test)]
mod tests {
    use super::Moon;
    use crate::{coordinate::Cartesian, reference_frame::ICRF};

    #[test]
    fn lunar_gravity_points_towards_center() {
        let moon = Moon::new();
        let r = Cartesian::<f64, ICRF<f64>>::new(moon.radius + 100_000.0, 0.0, 0.0);
        let g = moon.gravitational_acceleration(r);

        assert!(g.x() < 0.0);
        assert!(g.y().abs() < 1.0e-12);
        assert!(g.z().abs() < 1.0e-12);
    }

    #[test]
    fn lunar_gravity_magnitude_near_surface() {
        let moon = Moon::new();
        let r = Cartesian::<f64, ICRF<f64>>::new(moon.radius, 0.0, 0.0);
        let g = moon.gravitational_acceleration(r).norm();

        assert!((g - 1.62).abs() < 0.05);
    }
}
