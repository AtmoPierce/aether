use aether_core::{coordinate::Cartesian, reference_frame::NED};
use aether_models::terrestrial::wgs84::constants::{a as WGS84_SEMI_MAJOR_AXIS_M, f as WGS84_FLATTENING};
use core::fmt;

#[cfg(all(feature = "no_std", not(feature = "std")))]
use aether_core::real::Real;

#[cfg(feature = "std")]
use std::path::PathBuf;

pub const WMM_EXTRACTED_RELATIVE_DIR: &str = "catalog/wmm";
pub const WMM_MAX_DEGREE: usize = 12;
pub const WMM_RECORD_COUNT: usize = 90;
pub const WMM_TRIANGULAR_COUNT: usize = (WMM_MAX_DEGREE + 1) * (WMM_MAX_DEGREE + 2) / 2;
pub const WMM_ARRAY_DIM: usize = WMM_MAX_DEGREE + 1;
pub const WMM_REFERENCE_RADIUS_KM: f64 = 6_371.2;
pub const DEFAULT_POLE_EPSILON: f64 = 1.0e-12;

pub const fn wmm_reference_radius_m() -> f64 {
	WMM_REFERENCE_RADIUS_KM * 1_000.0
}

pub fn wgs84_mean_radius_m() -> f64 {
	WGS84_SEMI_MAJOR_AXIS_M * (1.0 - WGS84_FLATTENING / 3.0)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WmmHeader {
	pub epoch: f64,
	pub model_name: &'static str,
	pub release_date: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WmmCoefficientRecord {
	pub degree: usize,
	pub order: usize,
	pub g_nm_nt: f64,
	pub h_nm_nt: f64,
	pub g_dot_nm_nt_per_year: f64,
	pub h_dot_nm_nt_per_year: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeodeticPoint {
	pub latitude_rad: f64,
	pub longitude_rad: f64,
	pub height_m: f64,
}

impl GeodeticPoint {
	pub const fn new(latitude_rad: f64, longitude_rad: f64, height_m: f64) -> Self {
		Self {
			latitude_rad,
			longitude_rad,
			height_m,
		}
	}

	pub fn from_degrees(latitude_deg: f64, longitude_deg: f64, height_m: f64) -> Self {
		Self::new(latitude_deg.to_radians(), longitude_deg.to_radians(), height_m)
	}
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeocentricPoint {
	pub radius_m: f64,
	pub latitude_rad: f64,
	pub longitude_rad: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MagneticElements {
	pub geodetic: GeodeticPoint,
	pub geocentric: GeocentricPoint,
	pub decimal_year: f64,
	pub field_ned_nt: Cartesian<f64, NED<f64>>,
	pub secular_variation_ned_nt_per_year: Cartesian<f64, NED<f64>>,
	pub horizontal_intensity_nt: f64,
	pub total_intensity_nt: f64,
	pub declination_rad: f64,
	pub inclination_rad: f64,
	pub grid_variation_rad: f64,
	pub horizontal_intensity_rate_nt_per_year: f64,
	pub total_intensity_rate_nt_per_year: f64,
	pub declination_rate_rad_per_year: f64,
	pub inclination_rate_rad_per_year: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WmmError {
	#[cfg(feature = "std")]
	Io,
	PoleSingularity {
		latitude_rad: f64,
	},
	UndefinedHorizontalField,
}

impl fmt::Display for WmmError {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		match self {
			#[cfg(feature = "std")]
			Self::Io => write!(f, "failed to access WMM coefficient source"),
			Self::PoleSingularity { latitude_rad } => write!(
				f,
				"WMM spherical-harmonic evaluation is singular at geocentric latitude {latitude_rad} rad"
			),
			Self::UndefinedHorizontalField => write!(f, "horizontal magnetic field magnitude is zero"),
		}
	}
}

#[cfg(feature = "std")]
impl std::error::Error for WmmError {}

#[derive(Debug, Clone, Copy, PartialEq)]
struct TriangularTable {
	data: [f64; WMM_TRIANGULAR_COUNT],
}

impl TriangularTable {
	const fn from_array(data: [f64; WMM_TRIANGULAR_COUNT]) -> Self {
		Self { data }
	}

	const fn index(degree: usize, order: usize) -> usize {
		degree * (degree + 1) / 2 + order
	}

	fn set(&mut self, degree: usize, order: usize, value: f64) {
		self.data[Self::index(degree, order)] = value;
	}

	fn get(&self, degree: usize, order: usize) -> f64 {
		if degree > WMM_MAX_DEGREE || order > degree {
			0.0
		} else {
			self.data[Self::index(degree, order)]
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimeAdjustedCoefficients {
	pub decimal_year: f64,
	pub max_degree: usize,
	g: TriangularTable,
	h: TriangularTable,
	g_dot: TriangularTable,
	h_dot: TriangularTable,
}

impl TimeAdjustedCoefficients {
	pub fn g(&self, degree: usize, order: usize) -> f64 {
		self.g.get(degree, order)
	}

	pub fn h(&self, degree: usize, order: usize) -> f64 {
		self.h.get(degree, order)
	}

	pub fn g_dot(&self, degree: usize, order: usize) -> f64 {
		self.g_dot.get(degree, order)
	}

	pub fn h_dot(&self, degree: usize, order: usize) -> f64 {
		self.h_dot.get(degree, order)
	}
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WorldMagneticModel {
	pub header: WmmHeader,
	pub max_degree: usize,
	records: [WmmCoefficientRecord; WMM_RECORD_COUNT],
	g: TriangularTable,
	h: TriangularTable,
	g_dot: TriangularTable,
	h_dot: TriangularTable,
}

include!("wmm_generated.rs");

pub const WMM_2025: WorldMagneticModel = WorldMagneticModel {
	header: GENERATED_WMM_HEADER,
	max_degree: WMM_MAX_DEGREE,
	records: GENERATED_WMM_RECORDS,
	g: TriangularTable::from_array(GENERATED_WMM_G),
	h: TriangularTable::from_array(GENERATED_WMM_H),
	g_dot: TriangularTable::from_array(GENERATED_WMM_G_DOT),
	h_dot: TriangularTable::from_array(GENERATED_WMM_H_DOT),
};

impl Default for WorldMagneticModel {
	fn default() -> Self {
		WMM_2025
	}
}

impl WorldMagneticModel {
	pub const fn embedded() -> Self {
		WMM_2025
	}

	pub const fn records(&self) -> &[WmmCoefficientRecord; WMM_RECORD_COUNT] {
		&self.records
	}

	pub fn coefficients_at(&self, decimal_year: f64) -> TimeAdjustedCoefficients {
		let delta_years = decimal_year - self.header.epoch;
		let mut g = self.g;
		let mut h = self.h;

		for record in self.records.iter() {
			g.set(
				record.degree,
				record.order,
				record.g_nm_nt + delta_years * record.g_dot_nm_nt_per_year,
			);
			h.set(
				record.degree,
				record.order,
				record.h_nm_nt + delta_years * record.h_dot_nm_nt_per_year,
			);
		}

		TimeAdjustedCoefficients {
			decimal_year,
			max_degree: self.max_degree,
			g,
			h,
			g_dot: self.g_dot,
			h_dot: self.h_dot,
		}
	}

	pub fn evaluate(&self, decimal_year: f64, point: GeodeticPoint) -> Result<MagneticElements, WmmError> {
		let coefficients = self.coefficients_at(decimal_year);
		self.evaluate_with_coefficients(&coefficients, point)
	}

	pub fn evaluate_with_coefficients(
		&self,
		coefficients: &TimeAdjustedCoefficients,
		point: GeodeticPoint,
	) -> Result<MagneticElements, WmmError> {
		let geocentric = geodetic_to_geocentric(point);
		let sin_phi_prime = geocentric.latitude_rad.sin();
		let cos_phi_prime = geocentric.latitude_rad.cos();
		if cos_phi_prime.abs() <= DEFAULT_POLE_EPSILON {
			return Err(WmmError::PoleSingularity {
				latitude_rad: geocentric.latitude_rad,
			});
		}

		let (p_bar, dp_bar_dphi) = schmidt_semi_normalized_legendre(sin_phi_prime, cos_phi_prime);
		let (cos_m_lambda, sin_m_lambda) = longitude_harmonics(point.longitude_rad);
		let delta_lat = geocentric.latitude_rad - point.latitude_rad;

		let mut x_prime = 0.0;
		let mut y_prime = 0.0;
		let mut z_prime = 0.0;
		let mut x_prime_dot = 0.0;
		let mut y_prime_dot = 0.0;
		let mut z_prime_dot = 0.0;

		for degree in 1..=coefficients.max_degree {
			let radial_scale = (wmm_reference_radius_m() / geocentric.radius_m).powi((degree + 2) as i32);
			for order in 0..=degree {
				let cos_term = cos_m_lambda[order];
				let sin_term = sin_m_lambda[order];
				let g = coefficients.g(degree, order);
				let h = coefficients.h(degree, order);
				let g_dot = coefficients.g_dot(degree, order);
				let h_dot = coefficients.h_dot(degree, order);

				let combined = g * cos_term + h * sin_term;
				let combined_dot = g_dot * cos_term + h_dot * sin_term;
				let azimuthal = order as f64 * (g * sin_term - h * cos_term);
				let azimuthal_dot = order as f64 * (g_dot * sin_term - h_dot * cos_term);

				x_prime -= radial_scale * combined * dp_bar_dphi[degree][order];
				y_prime += radial_scale * azimuthal * p_bar[degree][order] / cos_phi_prime;
				z_prime -= radial_scale * (degree as f64 + 1.0) * combined * p_bar[degree][order];

				x_prime_dot -= radial_scale * combined_dot * dp_bar_dphi[degree][order];
				y_prime_dot += radial_scale * azimuthal_dot * p_bar[degree][order] / cos_phi_prime;
				z_prime_dot -= radial_scale * (degree as f64 + 1.0) * combined_dot * p_bar[degree][order];
			}
		}

		let cos_delta = delta_lat.cos();
		let sin_delta = delta_lat.sin();

		let x = x_prime * cos_delta - z_prime * sin_delta;
		let y = y_prime;
		let z = x_prime * sin_delta + z_prime * cos_delta;

		let x_dot = x_prime_dot * cos_delta - z_prime_dot * sin_delta;
		let y_dot = y_prime_dot;
		let z_dot = x_prime_dot * sin_delta + z_prime_dot * cos_delta;

		let horizontal_intensity_nt = (x * x + y * y).sqrt();
		if horizontal_intensity_nt == 0.0 {
			return Err(WmmError::UndefinedHorizontalField);
		}
		let total_intensity_nt = (horizontal_intensity_nt * horizontal_intensity_nt + z * z).sqrt();

		let declination_rad = y.atan2(x);
		let inclination_rad = z.atan2(horizontal_intensity_nt);
		let horizontal_intensity_rate_nt_per_year = (x * x_dot + y * y_dot) / horizontal_intensity_nt;
		let total_intensity_rate_nt_per_year =
			(x * x_dot + y * y_dot + z * z_dot) / total_intensity_nt;
		let declination_rate_rad_per_year =
			(x * y_dot - y * x_dot) / (horizontal_intensity_nt * horizontal_intensity_nt);
		let inclination_rate_rad_per_year =
			(horizontal_intensity_nt * z_dot - z * horizontal_intensity_rate_nt_per_year)
				/ (total_intensity_nt * total_intensity_nt);

		Ok(MagneticElements {
			geodetic: point,
			geocentric,
			decimal_year: coefficients.decimal_year,
			field_ned_nt: Cartesian::new(x, y, z),
			secular_variation_ned_nt_per_year: Cartesian::new(x_dot, y_dot, z_dot),
			horizontal_intensity_nt,
			total_intensity_nt,
			declination_rad,
			inclination_rad,
			grid_variation_rad: declination_rad,
			horizontal_intensity_rate_nt_per_year,
			total_intensity_rate_nt_per_year,
			declination_rate_rad_per_year,
			inclination_rate_rad_per_year,
		})
	}
}

#[cfg(feature = "std")]
const COEFFICIENT_FILE_CANDIDATES: [&str; 4] = [
	"catalog/wmm/WMM2025COF/WMM.COF",
	"catalog/wmm/WMM2025COF/WMM2025.COF",
	"catalog/wmm/WMM.COF",
	"catalog/wmm/WMM2025.COF",
];

#[cfg(feature = "std")]
pub fn default_coefficients_path() -> Option<PathBuf> {
	let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
	for relative in COEFFICIENT_FILE_CANDIDATES {
		let path = manifest_dir.join(relative);
		if path.is_file() {
			return Some(path);
		}
	}
	None
}

#[cfg(feature = "std")]
pub fn has_extracted_catalog() -> bool {
	default_coefficients_path().is_some()
}

fn geodetic_to_geocentric(point: GeodeticPoint) -> GeocentricPoint {
	let sin_lat = point.latitude_rad.sin();
	let cos_lat = point.latitude_rad.cos();
	let eccentricity_sq = WGS84_FLATTENING * (2.0 - WGS84_FLATTENING);
	let prime_vertical_radius =
		WGS84_SEMI_MAJOR_AXIS_M / (1.0 - eccentricity_sq * sin_lat * sin_lat).sqrt();
	let p = (prime_vertical_radius + point.height_m) * cos_lat;
	let z = (prime_vertical_radius * (1.0 - eccentricity_sq) + point.height_m) * sin_lat;
	let radius_m = (p * p + z * z).sqrt();
	let latitude_rad = (z / radius_m).asin();

	GeocentricPoint {
		radius_m,
		latitude_rad,
		longitude_rad: point.longitude_rad,
	}
}

fn schmidt_semi_normalized_legendre(
	sin_lat: f64,
	cos_lat: f64,
) -> (
	[[f64; WMM_ARRAY_DIM]; WMM_ARRAY_DIM],
	[[f64; WMM_ARRAY_DIM]; WMM_ARRAY_DIM],
) {
	let mut p = [[0.0; WMM_ARRAY_DIM]; WMM_ARRAY_DIM];
	p[0][0] = 1.0;

	for order in 1..=WMM_MAX_DEGREE {
		p[order][order] = ((2 * order - 1) as f64) * cos_lat * p[order - 1][order - 1];
	}

	for order in 0..WMM_MAX_DEGREE {
		p[order + 1][order] = ((2 * order + 1) as f64) * sin_lat * p[order][order];
	}

	for order in 0..=WMM_MAX_DEGREE {
		for degree in (order + 2)..=WMM_MAX_DEGREE {
			p[degree][order] = (((2 * degree - 1) as f64) * sin_lat * p[degree - 1][order]
				- ((degree + order - 1) as f64) * p[degree - 2][order])
				/ ((degree - order) as f64);
		}
	}

	let mut p_bar = [[0.0; WMM_ARRAY_DIM]; WMM_ARRAY_DIM];
	let mut dp_bar_dphi = [[0.0; WMM_ARRAY_DIM]; WMM_ARRAY_DIM];

	for degree in 0..=WMM_MAX_DEGREE {
		for order in 0..=degree {
			let schmidt = schmidt_factor(degree, order);
			p_bar[degree][order] = schmidt * p[degree][order];

			let previous = if degree > 0 && order < degree {
				p[degree - 1][order]
			} else {
				0.0
			};

			let derivative_wrt_x = if cos_lat.abs() <= DEFAULT_POLE_EPSILON {
				0.0
			} else {
				((degree as f64) * sin_lat * p[degree][order]
					- ((degree + order) as f64) * previous)
					/ (sin_lat * sin_lat - 1.0)
			};

			dp_bar_dphi[degree][order] = schmidt * cos_lat * derivative_wrt_x;
		}
	}

	(p_bar, dp_bar_dphi)
}

fn schmidt_factor(degree: usize, order: usize) -> f64 {
	if order == 0 {
		1.0
	} else {
		let mut ratio = 2.0;
		for value in (degree - order + 1)..=(degree + order) {
			ratio /= value as f64;
		}
		ratio.sqrt()
	}
}

fn longitude_harmonics(longitude_rad: f64) -> ([f64; WMM_ARRAY_DIM], [f64; WMM_ARRAY_DIM]) {
	let mut cos_m_lambda = [0.0; WMM_ARRAY_DIM];
	let mut sin_m_lambda = [0.0; WMM_ARRAY_DIM];
	cos_m_lambda[0] = 1.0;
	let cos_lambda = longitude_rad.cos();
	let sin_lambda = longitude_rad.sin();

	for order in 1..=WMM_MAX_DEGREE {
		cos_m_lambda[order] =
			cos_m_lambda[order - 1] * cos_lambda - sin_m_lambda[order - 1] * sin_lambda;
		sin_m_lambda[order] =
			sin_m_lambda[order - 1] * cos_lambda + cos_m_lambda[order - 1] * sin_lambda;
	}

	(cos_m_lambda, sin_m_lambda)
}

#[cfg(test)]
mod tests {
	use super::*;

	const WMM2025_TEST_VALUES: &str = include_str!("../../catalog/wmm/WMM2025COF/WMM2025_TestValues.txt");

	#[derive(Debug, Clone, Copy)]
	struct SampleCase {
		decimal_year: f64,
		height_km: f64,
		latitude_deg: f64,
		longitude_deg: f64,
		declination_deg: f64,
		inclination_deg: f64,
		h_nt: f64,
		x_nt: f64,
		y_nt: f64,
		z_nt: f64,
		f_nt: f64,
		d_declination_deg_per_year: f64,
		d_inclination_deg_per_year: f64,
		d_h_nt_per_year: f64,
		d_x_nt_per_year: f64,
		d_y_nt_per_year: f64,
		d_z_nt_per_year: f64,
		d_f_nt_per_year: f64,
	}

	#[test]
	fn exposes_embedded_wmm2025_coefficients() {
		let model = WorldMagneticModel::embedded();
		assert_eq!(model.header.epoch, 2025.0);
		assert_eq!(model.header.model_name, "WMM-2025");
		assert_eq!(model.header.release_date, "11/13/2024");
		assert_eq!(model.max_degree, WMM_MAX_DEGREE);
		assert_eq!(model.records().len(), WMM_RECORD_COUNT);
		assert!((model.coefficients_at(2025.0).g(1, 0) + 29_351.8).abs() < 1.0e-12);
		assert!((model.coefficients_at(2026.0).g(1, 0) + 29_339.8).abs() < 1.0e-9);
		assert!((wmm_reference_radius_m() - 6_371_200.0).abs() < 1.0e-12);
		assert!((wgs84_mean_radius_m() - 6_371_008.771_415_059).abs() < 1.0e-6);
	}

	#[test]
	fn reproduces_official_wmm2025_test_values() {
		let model = WorldMagneticModel::embedded();
		let cases = parse_test_values(WMM2025_TEST_VALUES);

		for case in &cases {
			let point = GeodeticPoint::from_degrees(
				case.latitude_deg,
				case.longitude_deg,
				case.height_km * 1_000.0,
			);
			let solution = model.evaluate(case.decimal_year, point).expect("evaluate sample");

			assert_close(solution.declination_rad.to_degrees(), case.declination_deg, 0.02, "declination");
			assert_close(solution.inclination_rad.to_degrees(), case.inclination_deg, 0.02, "inclination");
			assert_close(solution.horizontal_intensity_nt, case.h_nt, 0.1, "H");
			assert_close(solution.field_ned_nt.x(), case.x_nt, 0.1, "X");
			assert_close(solution.field_ned_nt.y(), case.y_nt, 0.1, "Y");
			assert_close(solution.field_ned_nt.z(), case.z_nt, 0.1, "Z");
			assert_close(solution.total_intensity_nt, case.f_nt, 0.1, "F");
			assert_close(
				solution.declination_rate_rad_per_year.to_degrees(),
				case.d_declination_deg_per_year,
				0.02,
				"Ddot",
			);
			assert_close(
				solution.inclination_rate_rad_per_year.to_degrees(),
				case.d_inclination_deg_per_year,
				0.02,
				"Idot",
			);
			assert_close(solution.horizontal_intensity_rate_nt_per_year, case.d_h_nt_per_year, 0.1, "Hdot");
			assert_close(solution.secular_variation_ned_nt_per_year.x(), case.d_x_nt_per_year, 0.1, "Xdot");
			assert_close(solution.secular_variation_ned_nt_per_year.y(), case.d_y_nt_per_year, 0.1, "Ydot");
			assert_close(solution.secular_variation_ned_nt_per_year.z(), case.d_z_nt_per_year, 0.1, "Zdot");
			assert_close(solution.total_intensity_rate_nt_per_year, case.d_f_nt_per_year, 0.1, "Fdot");
		}
	}

	fn parse_test_values(contents: &str) -> std::vec::Vec<SampleCase> {
		contents
			.lines()
			.filter_map(|line| {
				let line = line.trim();
				if line.is_empty() || line.starts_with('#') {
					return None;
				}

				let values: std::vec::Vec<f64> = line
					.split_whitespace()
					.map(|token| token.parse::<f64>().expect("numeric test value"))
					.collect();
				assert_eq!(values.len(), 18, "expected 18 columns in sample table");

				Some(SampleCase {
					decimal_year: values[0],
					height_km: values[1],
					latitude_deg: values[2],
					longitude_deg: values[3],
					declination_deg: values[4],
					inclination_deg: values[5],
					h_nt: values[6],
					x_nt: values[7],
					y_nt: values[8],
					z_nt: values[9],
					f_nt: values[10],
					d_declination_deg_per_year: values[11],
					d_inclination_deg_per_year: values[12],
					d_h_nt_per_year: values[13],
					d_x_nt_per_year: values[14],
					d_y_nt_per_year: values[15],
					d_z_nt_per_year: values[16],
					d_f_nt_per_year: values[17],
				})
			})
			.collect()
	}

	fn assert_close(actual: f64, expected: f64, tolerance: f64, label: &str) {
		let error = (actual - expected).abs();
		assert!(
			error <= tolerance,
			"{label} mismatch: actual={actual}, expected={expected}, abs_error={error}, tolerance={tolerance}"
		);
	}
}
