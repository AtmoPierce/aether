// Constants defined using https: f64 = //iers-conventions.obspm.fr/packaged_versions/iersconventions_v1_0_0.tar.gz
//  Chapter 1
// TODO Add uncertainty values to different constants, can likely make this a trait that automatically samples...
// TODO Traited units attached to these values would be nice.

// Natural defining constants
pub const speed_of_light: f64 = 299792458.0; // [m/s] Speed of Light
                                             // Auxillary Defining Constants
pub const gaussian_gravitational_constant: f64 = 1.720209895e-2; // Gaussian gravitational constant
pub const l_g: f64 = 6.969290134e-10; // 1−d(TT)/d(TCG)
pub const l_b: f64 = 1.550519768e-8; // 1−d(TDB)/d(TCB)
pub const tdb_0: f64 = -6.55e-5; // TDB−TCB at JD 2443144.5 TAI
pub const earth_angular_position_initial: f64 = 0.7790572732640; // [rev] Earth Rotation Angle (ERA) at J2000.0
pub const earth_angular_rate: f64 = 1.00273781191135448; // [rev/UT1day] Rate of Advance of Earth
                                                         // Natural measurable constant
pub const gravitational_constant: f64 = 6.67428e-11; // [m3/kg*s2] Constant of Gravitation
                                                     //Body Constants
pub const heliocentric_gravitational_constant: f64 = 1.32712442099e20; // [m3/s2] Heliocentric gravitational constant
pub const sun_dynamical_form_factor: f64 = 2.0e-7; // Dynamical form factor of the Sun
pub const moon_earth_mass_ratio: f64 = 0.0123000371; // Moon-Earth Mass Ratio
                                                     // Earth Constants
pub const geocentric_gravitational_constant: f64 = 3.986004418e14; // [m3/s2] Geocentric gravitational constant
pub const earth_equatorial_radius: f64 = 6378136.6; // [m] Equatorial radius of the earth
pub const earth_dynamical_form_factor: f64 = 1.0826359e-3; // Dynamical form factor of the Earth
pub const earth_flattening_factor: f64 = 298.25642; // Flattening factor of the Earth
pub const earth_mean_equatorial_gravity: f64 = 9.7803278; // [m/s2] Mean equatorial gravity
pub const earth_geoid_potential: f64 = 62636856.0; // [m2/s2] Potential of the geoid
pub const earth_geoid_potential_scale_factor: f64 = 6363672.6; // [m] Geopotential scale factor (GM_0/W_O)
pub const earth_dynamical_flattening: f64 = 3273795.0e-9; // Dynamical Flateening
                                                          // Initial value at J2000.0
pub const ecliptic_obliquity: f64 = 84381.406; // Obliquity of the ecliptic
                                               // Other Constants
pub const astronomical_unit: f64 = 1.49597870700e11; // [m] Astronmical unit
pub const l_c: f64 = 1.48082686741e-8; // Average value of 1-d(TCG)/d(TCB) (IDK what this is)
