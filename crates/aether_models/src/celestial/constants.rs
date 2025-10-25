pub const DEGREES_TO_RADIANS: f64 = core::f64::consts::PI / 180.0;
pub const RADIANS_TO_DEGREES: f64 = 180.0 / core::f64::consts::PI;
pub const DEGREES_TO_ARCMINUTES: f64 = 60.0;
pub const ARCMINUTES_TO_DEGREES: f64 = 1.0 / 60.0;
pub const DEGREES_TO_ARCSECONDS: f64 = 3600.0;
pub const ARCSECONDS_TO_DEGREES: f64 = 1.0 / 3600.0;
pub const DEGREES_TO_MILLIARCSECONDS: f64 = 3600000.0;
pub const MILLIARCSECONDS_TO_DEGREES: f64 = 1.0 / 3600000.0;
pub const DEGREES_TO_MICROARCSECONDS: f64 = 3600000000.0;
pub const MICROARCSECONDS_TO_DEGREES: f64 = 1.0 / 3600000000.0;

/* Pi */
pub const IAU_DPI: f64 = 3.141592653589793238462643;

/* 2Pi */
pub const IAU_D2PI: f64 = 6.283185307179586476925287;

/* Radians to degrees */
pub const IAU_DR2D: f64 = 57.29577951308232087679815;

/* Degrees to radians */
pub const IAU_DD2R: f64 = 1.745329251994329576923691e-2;

/* Radians to arcseconds */
pub const IAU_DR2AS: f64 = 206264.8062470963551564734;

/* Arcseconds to radians */
pub const IAU_DAS2R: f64 = 4.848136811095359935899141e-6;

/* Seconds of time to radians */
pub const IAU_DS2R: f64 = 7.272205216643039903848712e-5;

/* Arcseconds in a full circle */
pub const IAU_TURNAS: f64 = 1296000.0;

/* Milliarcseconds to radians */
pub const IAU_DMAS2R: f64 = IAU_DAS2R / 1e3;

/* Length of tropical year B1900 days */
pub const IAU_DTY: f64 = 365.242198781;

/* Seconds per day. */
pub const IAU_DAYSEC: f64 = 86400.0;

/* Days per Julian year */
pub const IAU_DJY: f64 = 365.25;

/* Days per Julian century */
pub const IAU_DJC: f64 = 36525.0;

/* Days per Julian millennium */
pub const IAU_DJM: f64 = 365250.0;

/* Reference epoch J2000.0, Julian Date */
pub const IAU_DJ00: f64 = 2451545.0;

/* Julian Date of Modified Julian Date zero */
pub const IAU_DJM0: f64 = 2400000.5;

/* Reference epoch J2000.0, Modified Julian Date */
pub const IAU_DJM00: f64 = 51544.5;

/* 1977 Jan 1.0 as MJD */
pub const IAU_DJM77: f64 = 43144.0;

/* TT minus TAI s */
pub const IAU_TTMTAI: f64 = 32.184;

/* Astronomical unit m, IAU 2012 */
pub const IAU_DAU: f64 = 149597870.7e3;

/* Speed of light m/s */
pub const IAU_CMPS: f64 = 299792458.0;

/* Light time for 1 au s */
pub const IAU_AULT: f64 = IAU_DAU / IAU_CMPS;

/* Speed of light au per day */
pub const IAU_DC: f64 = IAU_DAYSEC / IAU_AULT;

/* L_G = 1 - dTT/dTCG */
pub const IAU_ELG: f64 = 6.969290134e-10;

/* L_B = 1 - dTDB/dTCB, and TDB s at TAI 1977/1/1.0 */
pub const IAU_ELB: f64 = 1.550519768e-8;
pub const IAU_TDB0: f64 = -6.55e-5;

/* Schwarzschild radius of the Sun au */
/* = 2 * 1.32712440041e20 / 2.99792458e8^2 / 1.49597870700e11 */
pub const IAU_SRS: f64 = 1.97412574336e-8;
