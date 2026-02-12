use super::atmospheres::ussa;
use super::wgs84;
use crate::{
    attitude::{DirectionCosineMatrix, Euler},
    coordinate::Cartesian,
    math::{Matrix, Vector},
    reference_frame::{ICRF, ITRF},
};

#[derive(Clone, Debug)]
pub struct Earth {
    pub mass: f64,
    pub radius: f64,
    pub rotational_velocity: Cartesian<f64, ITRF<f64>>,
    pub gravitational_constant: f64,
    pub atmosphere: ussa::USSA,
}

impl Default for Earth {
    fn default() -> Self {
        Earth {
            mass: 0.0,
            radius: 0.0,
            rotational_velocity: Cartesian::new(0.0, 0.0, 7.292115146706979e-5),
            gravitational_constant: super::iers::constants::geocentric_gravitational_constant,
            atmosphere: ussa::USSA::new(),
        }
    }
}

impl Earth {
    pub fn new() -> Self {
        Earth {
            mass: super::iers::constants::geocentric_gravitational_constant
                / super::iers::constants::gravitational_constant,
            radius: super::iers::constants::earth_equatorial_radius,
            rotational_velocity: Cartesian::new(0.0, 0.0, 7.292115146706979e-5),
            gravitational_constant: super::iers::constants::geocentric_gravitational_constant,
            atmosphere: ussa::USSA::new(),
        }
    }

    // not currently in use as angular acceleration is not being modeled
    // pub fn solve_euler_force(
    //     &self,
    //     position_from_center: Cartesian<f64, ITRF<f64>>,
    //     mass: f64,
    // ) -> Cartesian<f64, ITRF<f64>> {
    //     let euler_force = -(angular_acceleration.cross(&position_from_center) * mass);
    //     return euler_force;
    // }

    pub fn solve_coriolis_force(
        &self,
        velocity_from_center: Cartesian<f64, ITRF<f64>>,
        mass: f64,
    ) -> Cartesian<f64, ITRF<f64>> {
        let coriolis_force = -(self.rotational_velocity.cross(&velocity_from_center) * 2.0 * mass);
        return coriolis_force;
    }

    pub fn solve_centrifugal_force(
        &self,
        position_from_center: Cartesian<f64, ITRF<f64>>,
        mass: f64,
    ) -> Cartesian<f64, ITRF<f64>> {
        let centrifugal_force =
            -(self.rotational_velocity.cross(&(self.rotational_velocity.cross(&position_from_center))) * mass);
        return centrifugal_force;
    }

    pub fn solve_gravitational_force(
        &self,
        position_from_center: Cartesian<f64, ITRF<f64>>,
        mass: f64,
    ) -> Cartesian<f64, ITRF<f64>> {
        let gravitational_force = super::wgs84::gravity::gravity_rectangular(
            position_from_center.x(),
            position_from_center.y(),
            position_from_center.z(),
        ) * mass;
        return gravitational_force;
    }
    pub fn solve_gravitational_torque(
        &self,
        position_from_center: Cartesian<f64, ITRF<f64>>,
        inertia_matrix: Matrix<f64, 3, 3>,
    ) -> Cartesian<f64, ITRF<f64>> {
        let value = inertia_matrix * position_from_center;
        let local = position_from_center.cross(&value);
        let scalar = (3.0 * super::iers::constants::geocentric_gravitational_constant)
            / (position_from_center.norm().powf(5.0));
        let gravitational_torque = local * scalar;
        return gravitational_torque;
    }
    pub fn geocentric_to_ecef(
        &self,
        latitude: f64,
        longitude: f64,
        altitude: f64,
    ) -> Cartesian<f64, ITRF<f64>> {
        return super::wgs84::transforms::geocentric_to_ecef(latitude, longitude, altitude);
    }

    pub fn ecef_to_geocentric_ferrari(&self, x: f64, y: f64, z: f64) -> Vector<f64, 3> {
        return super::wgs84::transforms::ecef_to_geocentric_ferrari(x, y, z);
    }

    pub fn ecef_to_geocentric(&self, x: f64, y: f64, z: f64) -> Vector<f64, 3> {
        return super::wgs84::transforms::ecef_to_geocentric(x, y, z);
    }

    pub fn icrf_to_itrf(&self, time: f64) -> DirectionCosineMatrix<f64, ICRF<f64>, ITRF<f64>> {
        return wgs84::transforms::icrf_to_itrf(time, self.rotational_velocity);
    }
    pub fn itrf_to_icrf(&self, time: f64) -> DirectionCosineMatrix<f64, ITRF<f64>, ICRF<f64>> {
        return wgs84::transforms::itrf_to_icrf(time, self.rotational_velocity);
    }

    pub fn get_temperature(&self, geometric_height: f64) -> Result<f64, &'static str> {
        return Ok(self
            .atmosphere
            .temperature(geometric_height)
            .expect("Could not get temperature."));
    }
    pub fn get_pressure(&self, geometric_height: f64) -> Result<f64, &'static str> {
        return Ok(self
            .atmosphere
            .pressure(geometric_height)
            .expect("Could not get pressure."));
    }
    pub fn get_density(&self, geometric_height: f64) -> Result<f64, &'static str> {
        return Ok(self
            .atmosphere
            .density(geometric_height)
            .expect("Could not get density."));
    }
    pub fn get_speed_of_sound(&self, geometric_height: f64) -> Result<f64, &'static str> {
        return Ok(self
            .atmosphere
            .speed_of_sound(geometric_height)
            .expect("Could not get speed of sound."));
    }
}

mod tests {
    // use hifitime::Epoch;
    // use crate::Earth;

    #[test]
    fn test_update_crs_to_trs_dcm() {
        // let mut earth = Earth::new();
        // let crs_to_trs_dcm = earth.update_crs_to_trs_dcm(Epoch::now().unwrap().to_et_seconds());
        // let theta_x = crs_to_trs_dcm.m32.atan2(crs_to_trs_dcm.m33);
        // let theta_y = -crs_to_trs_dcm.m31.atan2((crs_to_trs_dcm.m32.powf(2.0) + crs_to_trs_dcm.m33.powf(2.0)).sqrt());
        // let theta_z = crs_to_trs_dcm.m21.atan2(crs_to_trs_dcm.m11);
        // println!("Angles: {},{},{}", theta_x.to_degrees(), theta_y.to_degrees(), theta_z.to_degrees());
    }
}
