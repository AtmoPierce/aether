use aether_core::math::Vector;

#[cfg(feature = "pro")]
use aether_pro::step::{write_step_mesh_file, write_step_mesh_str, TriangleMesh};

pub type ProfilePoint2 = Vector<f64, 2>;
pub type ProfilePoint3 = Vector<f64, 3>;

#[derive(Debug, Clone, PartialEq)]
pub struct ClosedPolynomialProfile2D {
    pub x_start: f64,
    pub x_end: f64,
    pub upper_coefficients: Vec<f64>,
    pub lower_coefficients: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolynomialProfileSection2D {
    pub x_stations: Vec<f64>,
    pub upper_points: Vec<ProfilePoint2>,
    pub lower_points: Vec<ProfilePoint2>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlanarProjection3D {
    pub origin: ProfilePoint3,
    pub u_axis: ProfilePoint3,
    pub v_axis: ProfilePoint3,
}

#[cfg(feature = "pro")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProfileExtrusion3D {
    pub projection: PlanarProjection3D,
    pub extrusion: ProfilePoint3,
}

impl ClosedPolynomialProfile2D {
    pub fn new(
        x_start: f64,
        x_end: f64,
        upper_coefficients: Vec<f64>,
        lower_coefficients: Vec<f64>,
    ) -> Result<Self, String> {
        if !(x_start.is_finite() && x_end.is_finite()) {
            return Err("Polynomial profile x range must be finite".to_string());
        }
        if x_end <= x_start {
            return Err("Polynomial profile requires x_end > x_start".to_string());
        }
        if upper_coefficients.is_empty() || lower_coefficients.is_empty() {
            return Err(
                "Polynomial profile requires both upper and lower coefficient lists".to_string(),
            );
        }

        Ok(Self {
            x_start,
            x_end,
            upper_coefficients,
            lower_coefficients,
        })
    }

    pub fn evaluate_upper(&self, x: f64) -> f64 {
        evaluate_polynomial(&self.upper_coefficients, x)
    }

    pub fn evaluate_lower(&self, x: f64) -> f64 {
        evaluate_polynomial(&self.lower_coefficients, x)
    }

    pub fn sample_section(
        &self,
        station_count: usize,
    ) -> Result<PolynomialProfileSection2D, String> {
        if station_count < 2 {
            return Err("Polynomial profile sampling requires at least two stations".to_string());
        }

        let mut x_stations = Vec::with_capacity(station_count);
        let mut upper_points = Vec::with_capacity(station_count);
        let mut lower_points = Vec::with_capacity(station_count);

        for index in 0..station_count {
            let fraction = index as f64 / (station_count - 1) as f64;
            let x = self.x_start + fraction * (self.x_end - self.x_start);
            x_stations.push(x);
            upper_points.push(point2(x, self.evaluate_upper(x)));
            lower_points.push(point2(x, self.evaluate_lower(x)));
        }

        Ok(PolynomialProfileSection2D {
            x_stations,
            upper_points,
            lower_points,
        })
    }

    pub fn projected_outline(
        &self,
        station_count: usize,
        projection: PlanarProjection3D,
    ) -> Result<Vec<ProfilePoint3>, String> {
        Ok(self
            .sample_section(station_count)?
            .projected_closed_outline_points(projection))
    }

    #[cfg(feature = "pro")]
    pub fn extruded_mesh(
        &self,
        station_count: usize,
        extrusion: ProfileExtrusion3D,
    ) -> Result<TriangleMesh, String> {
        self.sample_section(station_count)?.extruded_mesh(extrusion)
    }

    #[cfg(feature = "pro")]
    pub fn step_mesh_str(
        &self,
        name: &str,
        station_count: usize,
        extrusion: ProfileExtrusion3D,
    ) -> Result<String, String> {
        let mesh = self.extruded_mesh(station_count, extrusion)?;
        write_step_mesh_str(name, &mesh)
    }

    #[cfg(feature = "pro")]
    pub fn write_step_mesh_file(
        &self,
        path: impl AsRef<std::path::Path>,
        name: &str,
        station_count: usize,
        extrusion: ProfileExtrusion3D,
    ) -> Result<(), String> {
        let mesh = self.extruded_mesh(station_count, extrusion)?;
        write_step_mesh_file(path, name, &mesh)
    }
}

impl PolynomialProfileSection2D {
    pub fn outline_points(&self) -> Vec<ProfilePoint2> {
        let mut outline = self.upper_points.iter().rev().copied().collect::<Vec<_>>();
        let lower_start = if self
            .lower_points
            .first()
            .zip(self.upper_points.first())
            .map(|(lower, upper)| lower == upper)
            .unwrap_or(false)
        {
            1
        } else {
            0
        };
        let lower_end = if self
            .lower_points
            .last()
            .zip(self.upper_points.last())
            .map(|(lower, upper)| lower == upper)
            .unwrap_or(false)
        {
            self.lower_points.len().saturating_sub(1)
        } else {
            self.lower_points.len()
        };
        outline.extend(self.lower_points[lower_start..lower_end].iter().copied());
        outline
    }

    pub fn closed_outline_points(&self) -> Vec<ProfilePoint2> {
        let mut outline = self.outline_points();
        if let Some(first) = outline.first().copied() {
            let needs_closure = outline.last().map(|last| *last != first).unwrap_or(false);
            if needs_closure {
                outline.push(first);
            }
        }
        outline
    }

    pub fn projected_closed_outline_points(
        &self,
        projection: PlanarProjection3D,
    ) -> Vec<ProfilePoint3> {
        self.closed_outline_points()
            .into_iter()
            .map(|point_local| projection.project_point(point_local))
            .collect()
    }

    #[cfg(feature = "pro")]
    pub fn extruded_mesh(&self, extrusion: ProfileExtrusion3D) -> Result<TriangleMesh, String> {
        if self.upper_points.len() != self.lower_points.len() || self.upper_points.len() < 2 {
            return Err(
                "Polynomial section requires matched upper and lower point lists".to_string(),
            );
        }

        let plane_normal = extrusion
            .projection
            .u_axis
            .cross(&extrusion.projection.v_axis);
        if plane_normal.norm() <= 1.0e-12 {
            return Err("Polynomial profile projection axes must span a valid plane".to_string());
        }
        if extrusion.extrusion.norm() <= 1.0e-12 {
            return Err("Polynomial profile extrusion vector must be non-zero".to_string());
        }
        if plane_normal.dot(&extrusion.extrusion).abs() <= 1.0e-12 {
            return Err(
                "Polynomial profile extrusion vector must not lie in the section plane".to_string(),
            );
        }

        let station_count = self.upper_points.len();
        let mut positions = Vec::with_capacity(station_count * 4);

        for point_local in &self.upper_points {
            positions.push(extrusion.projection.project_point(*point_local));
        }
        for point_local in &self.lower_points {
            positions.push(extrusion.projection.project_point(*point_local));
        }
        for point_local in &self.upper_points {
            positions.push(extrusion.projection.project_point(*point_local) + extrusion.extrusion);
        }
        for point_local in &self.lower_points {
            positions.push(extrusion.projection.project_point(*point_local) + extrusion.extrusion);
        }

        let mesh_center = positions
            .iter()
            .copied()
            .fold(point3(0.0, 0.0, 0.0), |accum, position| accum + position)
            / positions.len() as f64;

        let mut indices = Vec::new();
        for segment in 0..(station_count - 1) {
            let upper_front_0 = segment as u32;
            let upper_front_1 = (segment + 1) as u32;
            let lower_front_0 = (station_count + segment) as u32;
            let lower_front_1 = (station_count + segment + 1) as u32;

            let upper_back_0 = (2 * station_count + segment) as u32;
            let upper_back_1 = (2 * station_count + segment + 1) as u32;
            let lower_back_0 = (3 * station_count + segment) as u32;
            let lower_back_1 = (3 * station_count + segment + 1) as u32;

            push_triangle_outward(
                &mut indices,
                &positions,
                mesh_center,
                upper_front_0,
                lower_front_0,
                upper_front_1,
            );
            push_triangle_outward(
                &mut indices,
                &positions,
                mesh_center,
                upper_front_1,
                lower_front_0,
                lower_front_1,
            );

            push_triangle_outward(
                &mut indices,
                &positions,
                mesh_center,
                upper_back_0,
                upper_back_1,
                lower_back_0,
            );
            push_triangle_outward(
                &mut indices,
                &positions,
                mesh_center,
                upper_back_1,
                lower_back_1,
                lower_back_0,
            );
        }

        let front_outline =
            outline_indices_for_points(&self.upper_points, &self.lower_points, 0, station_count);
        let back_outline = outline_indices_for_points(
            &self.upper_points,
            &self.lower_points,
            2 * station_count,
            3 * station_count,
        );
        for edge in 0..front_outline.len() {
            let next = (edge + 1) % front_outline.len();
            let front_a = front_outline[edge];
            let front_b = front_outline[next];
            let back_a = back_outline[edge];
            let back_b = back_outline[next];

            push_triangle_outward(
                &mut indices,
                &positions,
                mesh_center,
                front_a,
                front_b,
                back_b,
            );
            push_triangle_outward(
                &mut indices,
                &positions,
                mesh_center,
                front_a,
                back_b,
                back_a,
            );
        }

        Ok(TriangleMesh { positions, indices })
    }
}

impl PlanarProjection3D {
    pub fn new(origin: ProfilePoint3, u_axis: ProfilePoint3, v_axis: ProfilePoint3) -> Self {
        Self {
            origin,
            u_axis,
            v_axis,
        }
    }

    pub fn project_point(&self, point_local: ProfilePoint2) -> ProfilePoint3 {
        self.origin + self.u_axis * point_local[0] + self.v_axis * point_local[1]
    }
}

impl Default for PlanarProjection3D {
    fn default() -> Self {
        Self {
            origin: point3(0.0, 0.0, 0.0),
            u_axis: point3(1.0, 0.0, 0.0),
            v_axis: point3(0.0, 1.0, 0.0),
        }
    }
}

#[cfg(feature = "pro")]
impl ProfileExtrusion3D {
    pub fn new(projection: PlanarProjection3D, extrusion: ProfilePoint3) -> Self {
        Self {
            projection,
            extrusion,
        }
    }
}

#[cfg(feature = "pro")]
impl Default for ProfileExtrusion3D {
    fn default() -> Self {
        Self {
            projection: PlanarProjection3D::default(),
            extrusion: point3(0.0, 0.0, 1.0),
        }
    }
}

fn evaluate_polynomial(coefficients: &[f64], x: f64) -> f64 {
    coefficients
        .iter()
        .rev()
        .fold(0.0, |accum, coefficient| accum * x + coefficient)
}

#[cfg(feature = "pro")]
fn outline_indices_for_points(
    upper_points: &[ProfilePoint2],
    lower_points: &[ProfilePoint2],
    upper_base: usize,
    lower_base: usize,
) -> Vec<u32> {
    let station_count = upper_points.len();
    let mut outline = Vec::with_capacity(2 * station_count);
    for index in (0..station_count).rev() {
        outline.push((upper_base + index) as u32);
    }

    let lower_start = if lower_points
        .first()
        .zip(upper_points.first())
        .map(|(lower, upper)| lower == upper)
        .unwrap_or(false)
    {
        1
    } else {
        0
    };
    let lower_end = if lower_points
        .last()
        .zip(upper_points.last())
        .map(|(lower, upper)| lower == upper)
        .unwrap_or(false)
    {
        station_count.saturating_sub(1)
    } else {
        station_count
    };

    for index in lower_start..lower_end {
        outline.push((lower_base + index) as u32);
    }
    outline
}

#[cfg(feature = "pro")]
fn push_triangle_outward(
    indices: &mut Vec<u32>,
    positions: &[ProfilePoint3],
    mesh_center: ProfilePoint3,
    a: u32,
    b: u32,
    c: u32,
) {
    let pa = positions[a as usize];
    let pb = positions[b as usize];
    let pc = positions[c as usize];
    let normal = (pb - pa).cross(&(pc - pa));
    let triangle_center = (pa + pb + pc) / 3.0;
    if normal.dot(&(triangle_center - mesh_center)) >= 0.0 {
        indices.extend_from_slice(&[a, b, c]);
    } else {
        indices.extend_from_slice(&[a, c, b]);
    }
}

fn point2(x: f64, y: f64) -> ProfilePoint2 {
    Vector::new([x, y])
}

fn point3(x: f64, y: f64, z: f64) -> ProfilePoint3 {
    Vector::new([x, y, z])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn example_profile() -> ClosedPolynomialProfile2D {
        ClosedPolynomialProfile2D::new(0.0, 1.0, vec![0.08, 0.10, -0.18], vec![-0.08, -0.10, 0.18])
            .unwrap()
    }

    #[test]
    fn samples_closed_outline_from_polynomials() {
        let section = example_profile().sample_section(9).unwrap();
        let outline = section.closed_outline_points();

        assert_eq!(section.upper_points.len(), 9);
        assert_eq!(section.lower_points.len(), 9);
        assert_eq!(outline.len(), 19);
        assert_eq!(outline.first(), outline.last());
    }

    #[test]
    fn projects_outline_to_requested_plane() {
        let projection = PlanarProjection3D::new(
            point3(1.0, 2.0, 3.0),
            point3(1.0, 0.0, 0.0),
            point3(0.0, 0.0, 1.0),
        );
        let projected = example_profile().projected_outline(7, projection).unwrap();

        assert!(projected
            .iter()
            .all(|point| (point[1] - 2.0).abs() < 1.0e-12));
        assert!(projected
            .iter()
            .any(|point| (point[2] - 3.0).abs() > 1.0e-12));
    }

    #[cfg(feature = "pro")]
    #[test]
    fn extruded_mesh_is_watertight() {
        let mesh = example_profile()
            .extruded_mesh(
                17,
                ProfileExtrusion3D::new(PlanarProjection3D::default(), point3(0.0, 0.0, 0.2)),
            )
            .unwrap();

        assert!(!mesh.is_empty());
        assert!(mesh.is_watertight());
        assert!(mesh.surface_area() > 0.0);
    }

    #[cfg(feature = "pro")]
    #[test]
    fn step_export_contains_step_header() {
        let step = example_profile()
            .step_mesh_str(
                "poly",
                17,
                ProfileExtrusion3D::new(PlanarProjection3D::default(), point3(0.0, 0.0, 0.2)),
            )
            .unwrap();

        assert!(step.contains("ISO-10303-21"));
    }
}
