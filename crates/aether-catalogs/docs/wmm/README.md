# World Magnetic Model 2025

Source:
- https://www.ncei.noaa.gov/sites/default/files/2024-12/WMM2025COF.zip

## Using WMM for Local Magnetometer Field Calibration

The WMM should be used as a reference for the **expected Earth magnetic field** at a known
location, altitude, and date. It does **not** replace normal magnetometer calibration, because it
does not model local disturbances from motors, wiring, batteries, steel, or nearby structures.

### Recommended process

1. Evaluate the local Earth field with the WMM using geodetic latitude, longitude, altitude, and
	decimal year.
2. Rotate the expected field from NED into the body or sensor frame using a current attitude
	estimate.
3. Compare the expected body-frame field against the magnetometer measurement.
4. Fit hard-iron bias and, if needed, soft-iron scale/misalignment terms over many orientations.
5. Apply declination from the WMM when converting magnetic heading into true heading.

### Typical equations

Expected local field in the navigation frame:

$$
\mathbf{b}^{ned}_{expected} = \operatorname{WMM}(\varphi, \lambda, h, t)
$$

Rotate into the body frame:

$$
\mathbf{b}^{body}_{expected} = R_{body \leftarrow ned}\,\mathbf{b}^{ned}_{expected}
$$

Simple magnetometer measurement model:

$$
\mathbf{m}_{raw} \approx S\,\mathbf{b}^{body}_{expected} + \mathbf{b}
$$

where:

- $\mathbf{b}$ is the hard-iron bias
- $S$ is the soft-iron scale/misalignment matrix

Calibrated reading:

$$
\mathbf{m}_{cal} = S^{-1}(\mathbf{m}_{raw} - \mathbf{b})
$$

True-heading correction using declination $D$:

$$
\psi_{true} = \psi_{mag} + D
$$

### Minimal crate usage

```rust
use aether_catalogs::wmm::{GeodeticPoint, WorldMagneticModel};

let model = WorldMagneticModel::embedded();
let point = GeodeticPoint::from_degrees(latitude_deg, longitude_deg, altitude_m);
let field = model.evaluate(decimal_year, point)?;

let b_ned_nt = field.field_ned_nt;
let declination_rad = field.declination_rad;
let inclination_rad = field.inclination_rad;
let total_intensity_nt = field.total_intensity_nt;
```

### Practical notes

- Use WMM to predict the nominal local field vector and field magnitude.
- Use that prediction together with attitude and a tumble dataset to estimate calibration terms.
- Compare measured field magnitude against `total_intensity_nt` to detect disturbances.
- Compare measured dip angle against `inclination_rad` as a consistency check.
- Apply `declination_rad` to convert magnetic heading to true heading.

### In this crate

- `WorldMagneticModel::embedded()` provides the embedded WMM-2025 coefficients.
- `GeodeticPoint::from_degrees()` constructs the evaluation point.
- `WorldMagneticModel::evaluate()` returns `MagneticElements`, including:
  - `field_ned_nt`
  - `declination_rad`
  - `inclination_rad`
  - `total_intensity_nt`
  - `secular_variation_ned_nt_per_year`

## Generating a CSV Grid and Plotting with Python

If you want a globe or map-style visualization, a simple workflow is:

1. Generate a regular latitude/longitude CSV grid from Rust.
2. Plot that CSV in Python.

### Generate the CSV

An example generator is available at [crates/aether_catalogs/examples/wmm_grid_csv.rs](crates/aether_catalogs/examples/wmm_grid_csv.rs).

Example:

```bash
cargo run -p aether_catalogs --example wmm_grid_csv --features "std world-magnetic-model" -- \
  crates/aether_catalogs/docs/wmm/output/wmm_450km.csv 2025.0 450.0 5.0 5.0 89.0
```

Positional arguments:

1. `output.csv`
2. `decimal_year` (default `2025.0`)
3. `altitude_km` (default `450.0`)
4. `lat_step_deg` (default `5.0`)
5. `lon_step_deg` (default `5.0`)
6. `max_abs_lat_deg` (default `89.0`)

The generated CSV includes:

- `latitude_deg`
- `longitude_deg`
- `x_nt`, `y_nt`, `z_nt`
- `h_nt`, `f_nt`
- `declination_deg`, `inclination_deg`
- secular-variation terms

### Plot the CSV in Python

A plotting script is available at [crates/aether_catalogs/docs/wmm/plot_wmm_csv.py](crates/aether_catalogs/docs/wmm/plot_wmm_csv.py).

Example:

```bash
/usr/bin/python3 crates/aether_catalogs/docs/wmm/plot_wmm_csv.py \
  crates/aether_catalogs/docs/wmm/output/wmm_450km.csv \
  --output crates/aether_catalogs/docs/wmm/output/wmm_450km_3d.png \
  --view 3d \
  --field f_nt \
  --quiver-step 4 \
  --saa-threshold 32000
```

If `--output` is omitted, the script now saves beside the CSV as
`<csv_stem>_<field>_<view>.png`.

Use `--show` to request an interactive window when a display is available.

For a 2D longitude/latitude map, either omit `--view` or pass `--view 2d`.

Notes:

- `--field` selects the scalar coloring field: `f_nt`, `h_nt`, or `z_nt`.
- `--view` selects either `2d` or `3d`. The default is `2d`.
- `--quiver-step` controls arrow density.
- In `3d` view, arrow density is additionally capped so arrows are drawn at most once every `5` latitude samples and `5` longitude samples.
- `--saa-threshold` optionally highlights low-total-field regions such as the South Atlantic Anomaly.
- `--show` opens the figure interactively when running in a GUI session.
- The script plots horizontal field direction arrows on top of the colored map or globe.

### Python dependencies

The plotting script expects:

- `matplotlib`
- `numpy`

## Analyzing Maximum Magnetic-North Delta

If you want to quantify how quickly magnetic north changes over adjacent latitude/longitude cells,
use [crates/aether_catalogs/docs/wmm/analyze_declination_delta.py](crates/aether_catalogs/docs/wmm/analyze_declination_delta.py).

Example:

```bash
./.venv/bin/python crates/aether_catalogs/docs/wmm/analyze_declination_delta.py \
  crates/aether_catalogs/docs/wmm/output/wmm_450km.csv
```

This reports the maximum adjacent change in `declination_deg` for:

- latitude neighbors
- longitude neighbors
- the overall maximum of the two

It also supports filtered reports by minimum horizontal field strength `h_nt`, since declination is
poorly conditioned near magnetic poles where the horizontal field approaches zero.

Example with custom thresholds:

```bash
./.venv/bin/python crates/aether_catalogs/docs/wmm/analyze_declination_delta.py \
  crates/aether_catalogs/docs/wmm/output/wmm_450km.csv \
  --h-min 0 1000 5000 10000
```

To also save a 2D heatmap of adjacent declination deltas:

```bash
./.venv/bin/python crates/aether_catalogs/docs/wmm/analyze_declination_delta.py \
  crates/aether_catalogs/docs/wmm/output/wmm_450km.csv \
  --plot \
  --plot-h-min 10000
```

If `--output` is omitted, the plot is saved as:

`<csv_stem>_declination_delta_hmin_<threshold>.png`

The heatmap includes three panels:

- latitude-adjacent `|ΔD|`
- longitude-adjacent `|ΔD|`
- pointwise maximum adjacent `|ΔD|`