use aether_catalogs::{GeodeticPoint, WorldMagneticModel};
use std::{
	env,
	fs::File,
	io::{BufWriter, Write},
	path::PathBuf,
	process,
};

const DEFAULT_DECIMAL_YEAR: f64 = 2025.0;
const DEFAULT_ALTITUDE_KM: f64 = 450.0;
const DEFAULT_LAT_STEP_DEG: f64 = 5.0;
const DEFAULT_LON_STEP_DEG: f64 = 5.0;
const DEFAULT_MAX_ABS_LAT_DEG: f64 = 89.0;

struct Config {
	output_path: PathBuf,
	decimal_year: f64,
	altitude_km: f64,
	lat_step_deg: f64,
	lon_step_deg: f64,
	max_abs_lat_deg: f64,
}

fn main() {
	let config = match parse_args(env::args().skip(1).collect()) {
		Ok(config) => config,
		Err(message) => {
			eprintln!("{message}");
			print_usage();
			process::exit(2);
		}
	};

	if let Err(error) = generate_csv(&config) {
		eprintln!("failed to generate WMM CSV: {error}");
		process::exit(1);
	}
}

fn parse_args(args: Vec<String>) -> Result<Config, String> {
	if args.is_empty() {
		return Err("missing required output CSV path".to_string());
	}
	if args.iter().any(|arg| arg == "-h" || arg == "--help") {
		print_usage();
		process::exit(0);
	}

	let output_path = PathBuf::from(&args[0]);
	let decimal_year = parse_optional_f64(args.get(1), DEFAULT_DECIMAL_YEAR, "decimal_year")?;
	let altitude_km = parse_optional_f64(args.get(2), DEFAULT_ALTITUDE_KM, "altitude_km")?;
	let lat_step_deg = parse_optional_f64(args.get(3), DEFAULT_LAT_STEP_DEG, "lat_step_deg")?;
	let lon_step_deg = parse_optional_f64(args.get(4), DEFAULT_LON_STEP_DEG, "lon_step_deg")?;
	let max_abs_lat_deg = parse_optional_f64(args.get(5), DEFAULT_MAX_ABS_LAT_DEG, "max_abs_lat_deg")?;

	if lat_step_deg <= 0.0 {
		return Err("lat_step_deg must be positive".to_string());
	}
	if lon_step_deg <= 0.0 {
		return Err("lon_step_deg must be positive".to_string());
	}
	if !(0.0..90.0).contains(&max_abs_lat_deg) {
		return Err("max_abs_lat_deg must be in [0, 90)".to_string());
	}

	Ok(Config {
		output_path,
		decimal_year,
		altitude_km,
		lat_step_deg,
		lon_step_deg,
		max_abs_lat_deg,
	})
}

fn parse_optional_f64(value: Option<&String>, default: f64, label: &str) -> Result<f64, String> {
	match value {
		Some(value) => value
			.parse::<f64>()
			.map_err(|error| format!("invalid {label} '{value}': {error}")),
		None => Ok(default),
	}
}

fn print_usage() {
	eprintln!(
		concat!(
			"Usage:\n",
			"  cargo run -p aether_catalogs --example wmm_grid_csv --features \"std world-magnetic-model\" -- \\\n",
			"      <output.csv> [decimal_year] [altitude_km] [lat_step_deg] [lon_step_deg] [max_abs_lat_deg]\n\n",
			"Defaults:\n",
			"  decimal_year    = 2025.0\n",
			"  altitude_km     = 450.0\n",
			"  lat_step_deg    = 5.0\n",
			"  lon_step_deg    = 5.0\n",
			"  max_abs_lat_deg = 89.0\n"
		)
	);
}

fn generate_csv(config: &Config) -> Result<(), Box<dyn std::error::Error>> {
	if let Some(parent) = config.output_path.parent() {
		if !parent.as_os_str().is_empty() {
			std::fs::create_dir_all(parent)?;
		}
	}

	let file = File::create(&config.output_path)?;
	let mut writer = BufWriter::new(file);
	let model = WorldMagneticModel::embedded();
	let latitudes = grid_values(-config.max_abs_lat_deg, config.max_abs_lat_deg, config.lat_step_deg);
	let longitudes = grid_values(-180.0, 180.0, config.lon_step_deg);

	writeln!(
		writer,
		"decimal_year,altitude_km,latitude_deg,longitude_deg,x_nt,y_nt,z_nt,h_nt,f_nt,declination_deg,inclination_deg,x_dot_nt_per_year,y_dot_nt_per_year,z_dot_nt_per_year,f_dot_nt_per_year"
	)?;

	for latitude_deg in latitudes {
		for longitude_deg in &longitudes {
			let point = GeodeticPoint::from_degrees(latitude_deg, *longitude_deg, config.altitude_km * 1_000.0);
			let field = model.evaluate(config.decimal_year, point)?;
			writeln!(
				writer,
				"{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
				config.decimal_year,
				config.altitude_km,
				latitude_deg,
				longitude_deg,
				field.field_ned_nt.x(),
				field.field_ned_nt.y(),
				field.field_ned_nt.z(),
				field.horizontal_intensity_nt,
				field.total_intensity_nt,
				field.declination_rad.to_degrees(),
				field.inclination_rad.to_degrees(),
				field.secular_variation_ned_nt_per_year.x(),
				field.secular_variation_ned_nt_per_year.y(),
				field.secular_variation_ned_nt_per_year.z(),
				field.total_intensity_rate_nt_per_year,
			)?;
		}
	}

	writer.flush()?;
	println!("wrote {}", config.output_path.display());
	Ok(())
}

fn grid_values(min_value: f64, max_value: f64, step: f64) -> Vec<f64> {
	let mut values = Vec::new();
	let mut value = min_value;
	while value <= max_value + 1.0e-9 {
		values.push((value * 1_000_000.0).round() / 1_000_000.0);
		value += step;
	}

	if let Some(last) = values.last_mut() {
		if (*last - max_value).abs() > 1.0e-9 {
			values.push(max_value);
		} else {
			*last = max_value;
		}
	}

	values
}
