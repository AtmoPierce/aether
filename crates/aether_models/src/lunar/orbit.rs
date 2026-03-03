use super::constants;

#[derive(Clone, Copy, Debug)]
pub struct MoonOrbitAroundEarth {
    pub semi_major_axis_m: f64,
    pub eccentricity: f64,
    pub inclination_rad: f64,
    pub sidereal_period_s: f64,
    pub mu_earth_m3_s2: f64,
    pub mu_moon_m3_s2: f64,
}

impl Default for MoonOrbitAroundEarth {
    fn default() -> Self {
        Self {
            semi_major_axis_m: constants::MOON_SEMIMAJOR_AXIS_M,
            eccentricity: constants::MOON_ORBIT_ECCENTRICITY,
            inclination_rad: constants::MOON_ORBIT_INCLINATION_RAD,
            sidereal_period_s: constants::MOON_SIDEREAL_ORBIT_PERIOD_S,
            mu_earth_m3_s2: crate::terrestrial::iers::constants::geocentric_gravitational_constant,
            mu_moon_m3_s2: constants::SELENOCENTRIC_GRAVITATIONAL_CONSTANT,
        }
    }
}

impl MoonOrbitAroundEarth {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn mu_total_m3_s2(&self) -> f64 {
        self.mu_earth_m3_s2 + self.mu_moon_m3_s2
    }

    pub fn mean_motion_rad_s(&self) -> f64 {
        core::f64::consts::TAU / self.sidereal_period_s
    }

    pub fn keplerian_period_s(&self) -> f64 {
        core::f64::consts::TAU * (self.semi_major_axis_m.powi(3) / self.mu_total_m3_s2()).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::MoonOrbitAroundEarth;

    #[test]
    fn mean_motion_and_period_are_consistent() {
        let orbit = MoonOrbitAroundEarth::new();
        let reconstructed = core::f64::consts::TAU / orbit.mean_motion_rad_s();
        assert!((reconstructed - orbit.sidereal_period_s).abs() < 1.0e-9);
    }

    #[test]
    fn keplerian_period_matches_mean_period_order() {
        let orbit = MoonOrbitAroundEarth::new();
        let rel_err =
            ((orbit.keplerian_period_s() - orbit.sidereal_period_s) / orbit.sidereal_period_s).abs();
        assert!(rel_err < 0.05);
    }
}
