//! Lunar physical constants used by the `Moon` gravity model.
//!
//! Sources / conventions:
//! - Selenocentric gravitational parameter (`GM`): JPL SSD Astrodynamic Parameters
//!   (DE440), Moon = 4902.800118 km^3/s^2.
//!   https://ssd.jpl.nasa.gov/astro_par.html
//! - Mean radius: 1737.4 km (IAU/NASA standard planetary facts usage).
//! - Sidereal rotation period: 27.321661 days (used to derive angular rate).
//!
//! Notes:
//! - `MOON_MASS` is derived from `GM / G` to remain internally consistent with
//!   the crate's gravitational constant in `terrestrial::iers::constants`.

pub const SELENOCENTRIC_GRAVITATIONAL_CONSTANT: f64 = 4.902_800_118e12; // [m^3/s^2]
pub const MOON_MEAN_RADIUS: f64 = 1_737_400.0; // [m]
pub const MOON_SIDEREAL_ROTATION_PERIOD_DAYS: f64 = 27.321_661; // [day]
pub const MOON_SIDEREAL_ROTATION_RATE_RAD_S: f64 =
	core::f64::consts::TAU / (MOON_SIDEREAL_ROTATION_PERIOD_DAYS * 86_400.0); // [rad/s]
pub const MOON_MASS: f64 =
	SELENOCENTRIC_GRAVITATIONAL_CONSTANT / crate::terrestrial::iers::constants::gravitational_constant; // [kg]

// Mean lunar orbit around Earth (long-term averaged values, not osculating elements)
pub const MOON_SEMIMAJOR_AXIS_M: f64 = 384_400_000.0; // [m]
pub const MOON_ORBIT_ECCENTRICITY: f64 = 0.0549; // [-]
pub const MOON_ORBIT_INCLINATION_DEG: f64 = 5.145; // [deg] to ecliptic
pub const MOON_ORBIT_INCLINATION_RAD: f64 =
	MOON_ORBIT_INCLINATION_DEG * core::f64::consts::PI / 180.0; // [rad]
pub const MOON_SIDEREAL_ORBIT_PERIOD_DAYS: f64 = 27.321_661; // [day]
pub const MOON_SIDEREAL_ORBIT_PERIOD_S: f64 = MOON_SIDEREAL_ORBIT_PERIOD_DAYS * 86_400.0; // [s]
