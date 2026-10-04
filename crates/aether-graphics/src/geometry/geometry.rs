#[cfg(feature = "pro")]
use super::polynomial_profile::{ClosedPolynomialProfile2D, ProfileExtrusion3D};
use aether_core::math::{Matrix, Vector};
#[cfg(feature = "pro")]
use aether_pro::step::{
    load_step_assembly_part_file, load_step_part_file, step_length_unit_scale_to_si_file,
    tessellate_step_file, FacePatch, Part, TessellationOptions, TriangleMesh,
};
#[cfg(feature = "pro")]
use bincode::{
    config::standard,
    serde::{decode_from_slice, encode_to_vec},
};
#[cfg(feature = "pro")]
use serde::{Deserialize, Serialize};
use std::io::BufReader;
use std::path::Path;

pub type Vec2<T> = Vector<T, 2>;
pub type Vec3<T> = Vector<T, 3>;
pub type Vec4<T> = Vector<T, 4>;
pub type Mat4<T> = Matrix<T, 4, 4>;

#[cfg(feature = "pro")]
const STEP_GEOMETRY_CACHE_VERSION: u32 = 12;
#[cfg(feature = "pro")]
const LARGE_STEP_PREVIEW_BYTES: u64 = 32 * 1024 * 1024;
#[cfg(not(feature = "pro"))]
const STEP_SUPPORT_ERROR: &str = "STEP loading requires the aether_graphics `pro` feature";

#[cfg(feature = "pro")]
#[derive(Debug, Clone, Serialize, Deserialize)]
struct GeometryDiskCache {
    version: u32,
    source_len: u64,
    source_modified_ms: u128,
    positions: Vec<[f32; 3]>,
    indices: Vec<u32>,
    colors: Vec<[f32; 4]>,
    texture: Vec<[f32; 2]>,
    x_offset: f32,
    y_offset: f32,
    z_offset: f32,
}

#[cfg(feature = "pro")]
impl GeometryDiskCache {
    fn from_geometry(geometry: &Geometry, source_len: u64, source_modified_ms: u128) -> Self {
        Self {
            version: STEP_GEOMETRY_CACHE_VERSION,
            source_len,
            source_modified_ms,
            positions: geometry
                .positions
                .iter()
                .map(|p| [p[0], p[1], p[2]])
                .collect(),
            indices: geometry.indices.clone(),
            colors: geometry
                .colors
                .iter()
                .map(|c| [c[0], c[1], c[2], c[3]])
                .collect(),
            texture: geometry.texture.iter().map(|uv| [uv[0], uv[1]]).collect(),
            x_offset: geometry.x_offset,
            y_offset: geometry.y_offset,
            z_offset: geometry.z_offset,
        }
    }

    fn into_geometry(self) -> Geometry {
        let mut geometry = Geometry {
            positions: self.positions.into_iter().map(Vec3::new).collect(),
            indices: self.indices,
            colors: self.colors.into_iter().map(Vec4::new).collect(),
            texture: self.texture.into_iter().map(Vec2::new).collect(),
            elements: Vec::new(),
            x_offset: self.x_offset,
            y_offset: self.y_offset,
            z_offset: self.z_offset,
        };

        geometry.elements = geometry
            .positions
            .iter()
            .zip(geometry.colors.iter())
            .zip(geometry.texture.iter())
            .map(|((p, c), uv)| {
                Element::new(p[0], p[1], p[2], c[0], c[1], c[2], c[3], uv[0], uv[1])
            })
            .collect();

        geometry
    }
}

#[derive(Copy, Clone, Debug)]
#[repr(C, packed)]
pub struct Element {
    x: f32,
    y: f32,
    z: f32,
    r: f32,
    g: f32,
    b: f32,
    a: f32,
    s: f32,
    t: f32,
}

impl Element {
    #[allow(clippy::too_many_arguments)]
    pub fn new(x: f32, y: f32, z: f32, r: f32, g: f32, b: f32, a: f32, s: f32, t: f32) -> Element {
        Element {
            x,
            y,
            z,
            r,
            g,
            b,
            a,
            s,
            t,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Geometry {
    pub positions: Vec<Vec3<f32>>,
    pub indices: Vec<u32>,
    pub colors: Vec<Vec4<f32>>,
    pub texture: Vec<Vec2<f32>>,
    pub elements: Vec<Element>,
    pub x_offset: f32,
    pub y_offset: f32,
    pub z_offset: f32,
}

impl Default for Geometry {
    fn default() -> Self {
        Geometry {
            positions: Vec::default(),
            indices: Vec::default(),
            colors: Vec::default(),
            elements: Vec::default(),
            texture: Vec::default(),
            x_offset: 0.0,
            y_offset: 0.0,
            z_offset: 0.0,
        }
    }
}

impl Geometry {
    #[cfg(feature = "pro")]
    pub fn from_triangle_mesh(mesh: TriangleMesh) -> Result<Self, String> {
        let mut geometry = Self::default();
        geometry.load_triangle_mesh(mesh)?;
        Ok(geometry)
    }

    #[cfg(feature = "pro")]
    pub fn from_polynomial_profile(
        profile: &ClosedPolynomialProfile2D,
        station_count: usize,
        extrusion: ProfileExtrusion3D,
    ) -> Result<Self, String> {
        Self::from_triangle_mesh(profile.extruded_mesh(station_count, extrusion)?)
    }

    #[cfg(feature = "pro")]
    fn step_preview_tessellation_options(path: &Path) -> Result<TessellationOptions, String> {
        let source_len = std::fs::metadata(path)
            .map_err(|e| format!("Failed to read metadata for {}: {}", path.display(), e))?
            .len();

        if source_len >= LARGE_STEP_PREVIEW_BYTES {
            Ok(TessellationOptions {
                circle_segments: 16,
                line_segments: 2,
            })
        } else {
            Ok(TessellationOptions {
                circle_segments: 24,
                line_segments: 2,
            })
        }
    }

    #[cfg(feature = "pro")]
    fn should_use_fast_step_preview(path: &Path) -> Result<bool, String> {
        Ok(std::fs::metadata(path)
            .map_err(|e| format!("Failed to read metadata for {}: {}", path.display(), e))?
            .len()
            >= LARGE_STEP_PREVIEW_BYTES)
    }

    #[cfg(feature = "pro")]
    fn file_signature(path: &Path) -> Result<(u64, u128), String> {
        let metadata = std::fs::metadata(path)
            .map_err(|e| format!("Failed to read metadata for {}: {}", path.display(), e))?;
        let len = metadata.len();
        let modified = metadata
            .modified()
            .map_err(|e| format!("Failed to read modified time for {}: {}", path.display(), e))?;
        let modified_ms = modified
            .duration_since(std::time::UNIX_EPOCH)
            .map_err(|e| {
                format!(
                    "Failed to convert modified time for {}: {}",
                    path.display(),
                    e
                )
            })?
            .as_millis();
        Ok((len, modified_ms))
    }

    #[cfg(feature = "pro")]
    fn step_geometry_cache_path(path: &Path) -> Result<std::path::PathBuf, String> {
        let file_name = path
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| format!("Invalid cache file name for {}", path.display()))?;
        Ok(path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(format!(".{}.telegraph-stepmesh.bin", file_name)))
    }

    #[cfg(feature = "pro")]
    fn try_load_step_cache(path: &Path) -> Result<Option<Geometry>, String> {
        let cache_path = Self::step_geometry_cache_path(path)?;
        if !cache_path.exists() {
            return Ok(None);
        }

        let bytes = match std::fs::read(&cache_path) {
            Ok(bytes) => bytes,
            Err(_) => return Ok(None),
        };
        let cache: GeometryDiskCache =
            match decode_from_slice::<GeometryDiskCache, _>(&bytes, standard()) {
                Ok((cache, _)) => cache,
                Err(_) => return Ok(None),
            };
        if cache.version != STEP_GEOMETRY_CACHE_VERSION {
            return Ok(None);
        }

        let (source_len, source_modified_ms) = Self::file_signature(path)?;
        if cache.source_len != source_len || cache.source_modified_ms != source_modified_ms {
            return Ok(None);
        }

        Ok(Some(cache.into_geometry()))
    }

    #[cfg(feature = "pro")]
    fn write_step_cache(&self, path: &Path) -> Result<(), String> {
        let cache_path = Self::step_geometry_cache_path(path)?;
        let (source_len, source_modified_ms) = Self::file_signature(path)?;
        let cache = GeometryDiskCache::from_geometry(self, source_len, source_modified_ms);
        let bytes = encode_to_vec(&cache, standard()).map_err(|e| {
            format!(
                "Failed to serialize STEP cache for {}: {}",
                path.display(),
                e
            )
        })?;
        std::fs::write(&cache_path, bytes)
            .map_err(|e| format!("Failed to write STEP cache {}: {}", cache_path.display(), e))
    }

    pub fn tint(&mut self, r: f32, g: f32, b: f32, a: f32) {
        for el in &mut self.elements {
            el.r = r;
            el.g = g;
            el.b = b;
            el.a = a;
        }

        if !self.colors.is_empty() {
            for c in &mut self.colors {
                c[0] = r;
                c[1] = g;
                c[2] = b;
                c[3] = a;
            }
        }
    }

    pub fn flatten_elements(&self) -> Vec<f32> {
        self.elements
            .clone()
            .into_iter()
            .flat_map(|p| {
                core::iter::once(p.x).chain(
                    core::iter::once(p.y).chain(
                        core::iter::once(p.z).chain(
                            core::iter::once(p.r).chain(
                                core::iter::once(p.g).chain(
                                    core::iter::once(p.b).chain(
                                        core::iter::once(p.a).chain(
                                            core::iter::once(p.s).chain(core::iter::once(p.t)),
                                        ),
                                    ),
                                ),
                            ),
                        ),
                    ),
                )
            })
            .collect()
    }

    pub fn load_mesh_file<P: AsRef<Path>>(&mut self, path: P) -> Result<(), String> {
        let path = path.as_ref();
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.to_ascii_lowercase())
            .ok_or_else(|| "Mesh file has no extension".to_string())?;

        match ext.as_str() {
            "stl" => self.load_stl(path),
            "stp" | "step" => self.load_step(path),
            _ => Err(format!(
                "Unsupported mesh extension `.{}`. Supported: .stl, .stp, .step",
                ext
            )),
        }
    }

    fn reset_mesh_buffers(&mut self) {
        self.positions.clear();
        self.indices.clear();
        self.colors.clear();
        self.texture.clear();
        self.elements.clear();
    }

    fn push_vertex(&mut self, x: f32, y: f32, z: f32) {
        self.positions.push(Vec3::new([x, y, z]));
        self.colors.push(Vec4::new([1.0, 1.0, 1.0, 1.0]));
        self.texture.push(Vec2::new([0.0, 0.0]));
        self.elements
            .push(Element::new(x, y, z, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0));
    }

    pub fn recenter_to_origin(&mut self) {
        if self.positions.is_empty() {
            return;
        }

        let mut min_x = self.positions[0][0];
        let mut min_y = self.positions[0][1];
        let mut min_z = self.positions[0][2];
        let mut max_x = self.positions[0][0];
        let mut max_y = self.positions[0][1];
        let mut max_z = self.positions[0][2];

        for p in &self.positions {
            min_x = min_x.min(p[0]);
            min_y = min_y.min(p[1]);
            min_z = min_z.min(p[2]);
            max_x = max_x.max(p[0]);
            max_y = max_y.max(p[1]);
            max_z = max_z.max(p[2]);
        }

        let cx = 0.5 * (min_x + max_x);
        let cy = 0.5 * (min_y + max_y);
        let cz = 0.5 * (min_z + max_z);

        for p in &mut self.positions {
            p[0] -= cx;
            p[1] -= cy;
            p[2] -= cz;
        }

        for el in &mut self.elements {
            el.x -= cx;
            el.y -= cy;
            el.z -= cz;
        }
    }

    fn load_stl(&mut self, path: &Path) -> Result<(), String> {
        let file = std::fs::File::open(path)
            .map_err(|e| format!("Failed to open STL file {}: {}", path.display(), e))?;
        let mut reader = BufReader::new(file);
        let mesh = stl_io::read_stl(&mut reader)
            .map_err(|e| format!("Failed to parse STL {}: {}", path.display(), e))?;

        self.reset_mesh_buffers();

        for vertex in &mesh.vertices {
            self.push_vertex(vertex[0], vertex[1], vertex[2]);
        }

        for face in &mesh.faces {
            for &idx in &face.vertices {
                let index_u32 = u32::try_from(idx)
                    .map_err(|_| format!("Mesh index {} exceeds u32 range", idx))?;
                if index_u32 as usize >= self.positions.len() {
                    return Err(format!(
                        "STL file {} has invalid index {} for {} vertices",
                        path.display(),
                        index_u32,
                        self.positions.len()
                    ));
                }
                self.indices.push(index_u32);
            }
        }

        if self.indices.is_empty() {
            return Err(format!("STL file {} has no triangles", path.display()));
        }

        Ok(())
    }

    #[cfg(feature = "pro")]
    fn load_step(&mut self, path: &Path) -> Result<(), String> {
        if let Some(cached) = Self::try_load_step_cache(path)? {
            *self = cached;
            return Ok(());
        }

        let preview_options = Self::step_preview_tessellation_options(path)?;
        let scale_to_m = step_length_unit_scale_to_si_file(path)
            .map_err(|e| format!("aether_pro::step unit read failed: {e}"))?
            .unwrap_or(1.0);

        if Self::should_use_fast_step_preview(path)? {
            let mut mesh = tessellate_step_file(path, preview_options)
                .map_err(|e| format!("aether_pro::step mesh load failed: {e}"))?;
            if (scale_to_m - 1.0).abs() > 1.0e-12 {
                mesh = mesh.scaled(scale_to_m);
            }
            self.load_triangle_mesh(mesh)?;
        } else {
            let mut part = load_step_part_file(path, preview_options)
                .map_err(|e| format!("aether_pro::step part load failed: {e}"))?;
            if (scale_to_m - 1.0).abs() > 1.0e-12 {
                part.mesh = part.mesh.scaled(scale_to_m);
            }
            self.load_part(part)?;
        }
        let _ = self.write_step_cache(path);
        Ok(())
    }

    #[cfg(not(feature = "pro"))]
    fn load_step(&mut self, path: &Path) -> Result<(), String> {
        let _ = path;
        Err(STEP_SUPPORT_ERROR.to_string())
    }

    #[cfg(feature = "pro")]
    pub fn load_step_assembly_toml(&mut self, path: &Path) -> Result<(), String> {
        let part = load_step_assembly_part_file(path)
            .map_err(|e| format!("aether_pro::step assembly load failed: {e}"))?;
        self.load_part(part)
    }

    #[cfg(not(feature = "pro"))]
    pub fn load_step_assembly_toml(&mut self, path: &Path) -> Result<(), String> {
        let _ = path;
        Err(STEP_SUPPORT_ERROR.to_string())
    }

    #[cfg(feature = "pro")]
    fn load_part(&mut self, part: Part) -> Result<(), String> {
        self.reset_mesh_buffers();

        let default_color = part
            .appearance
            .base_color_rgba
            .unwrap_or([1.0, 1.0, 1.0, 1.0]);

        if part.face_patches.is_empty() {
            self.load_triangle_mesh_with_color(part.mesh, default_color)
        } else {
            self.load_triangle_mesh_patches(part.mesh, &part.face_patches, default_color)
        }
    }

    #[cfg(feature = "pro")]
    fn load_triangle_mesh(&mut self, mesh: TriangleMesh) -> Result<(), String> {
        self.load_triangle_mesh_with_color(mesh, [1.0, 1.0, 1.0, 1.0])
    }

    #[cfg(feature = "pro")]
    fn load_triangle_mesh_with_color(
        &mut self,
        mesh: TriangleMesh,
        color: [f32; 4],
    ) -> Result<(), String> {
        self.reset_mesh_buffers();

        for position in mesh.positions {
            let x = position[0] as f32;
            let y = position[1] as f32;
            let z = position[2] as f32;
            self.positions.push(Vec3::new([x, y, z]));
            self.colors.push(Vec4::new(color));
            self.texture.push(Vec2::new([0.0, 0.0]));
            self.elements.push(Element::new(
                x, y, z, color[0], color[1], color[2], color[3], 0.0, 0.0,
            ));
        }

        for index in mesh.indices {
            if index as usize >= self.positions.len() {
                return Err(format!(
                    "Triangle mesh index {} exceeds vertex count {}",
                    index,
                    self.positions.len()
                ));
            }
            self.indices.push(index);
        }

        if self.indices.is_empty() {
            return Err("Triangle mesh has no triangles".to_string());
        }

        Ok(())
    }

    #[cfg(feature = "pro")]
    fn load_triangle_mesh_patches(
        &mut self,
        mesh: TriangleMesh,
        face_patches: &[FacePatch],
        default_color: [f32; 4],
    ) -> Result<(), String> {
        self.reset_mesh_buffers();

        let triangle_count = mesh.indices.len() / 3;
        let mut patch_for_triangle = vec![default_color; triangle_count];

        for patch in face_patches {
            let color = patch.appearance.base_color_rgba.unwrap_or(default_color);
            let end = (patch.first_triangle + patch.triangle_count).min(triangle_count);
            for tri in patch.first_triangle..end {
                patch_for_triangle[tri] = color;
            }
        }

        for (tri_idx, tri) in mesh.indices.chunks_exact(3).enumerate() {
            let color = patch_for_triangle[tri_idx];
            let base_index = u32::try_from(self.positions.len())
                .map_err(|_| "Expanded mesh vertex count exceeds u32 range".to_string())?;

            for &index in tri {
                let position = *mesh.positions.get(index as usize).ok_or_else(|| {
                    format!(
                        "Triangle mesh index {} exceeds vertex count {}",
                        index,
                        mesh.positions.len()
                    )
                })?;
                let x = position[0] as f32;
                let y = position[1] as f32;
                let z = position[2] as f32;
                self.positions.push(Vec3::new([x, y, z]));
                self.colors.push(Vec4::new(color));
                self.texture.push(Vec2::new([0.0, 0.0]));
                self.elements.push(Element::new(
                    x, y, z, color[0], color[1], color[2], color[3], 0.0, 0.0,
                ));
            }

            self.indices
                .extend_from_slice(&[base_index, base_index + 1, base_index + 2]);
        }

        if self.indices.is_empty() {
            return Err("Triangle mesh has no triangles".to_string());
        }

        Ok(())
    }

    pub fn cone(
        &mut self,
        diameter: f32,
        length: f32,
        length_subdivisions: u32,
        angle_subdivisions: u32,
    ) {
        let radius = diameter / 2.0;
        self.positions.clear();
        for i in 0..length_subdivisions + 1 {
            let x = i as f32 / length_subdivisions as f32;
            for j in 0..angle_subdivisions {
                let angle = 2.0 * std::f32::consts::PI * j as f32 / angle_subdivisions as f32;
                self.positions.push(Vec3::new([
                    radius * x * angle.cos(),
                    length * (radius - x),
                    radius * x * angle.sin(),
                ]));
            }
        }
        self.indices.clear();

        for i in 0..length_subdivisions {
            for j in 0..angle_subdivisions {
                self.indices.push(i * angle_subdivisions + j);
                self.indices
                    .push(i * angle_subdivisions + (j + 1) % angle_subdivisions);
                self.indices
                    .push((i + 1) * angle_subdivisions + (j + 1) % angle_subdivisions);

                self.indices.push(i * angle_subdivisions + j);
                self.indices
                    .push((i + 1) * angle_subdivisions + (j + 1) % angle_subdivisions);
                self.indices.push((i + 1) * angle_subdivisions + j);

                self.colors.push(Vec4::new([1.0, 0.0, 0.0, 1.0]));
                self.colors.push(Vec4::new([0.0, 1.0, 0.0, 1.0]));
                self.colors.push(Vec4::new([0.0, 0.0, 1.0, 1.0]));
                self.colors.push(Vec4::new([1.0, 1.0, 0.0, 1.0]));
                self.colors.push(Vec4::new([0.0, 0.0, 1.0, 1.0]));
                self.colors.push(Vec4::new([1.0, 0.0, 0.0, 1.0]));
            }
        }
    }

    pub fn cone_blunted_spherical(
        &mut self,
        diameter: f32,
        length: f32,
        x_tangency_point: f32,
        length_subdivisions: u32,
        angle_subdivisions: u32,
    ) {
        self.reset_mesh_buffers();

        if diameter <= 0.0
            || length <= 0.0
            || x_tangency_point <= 0.0
            || x_tangency_point >= length
            || length_subdivisions == 0
            || angle_subdivisions < 3
        {
            return;
        }

        let base_radius = diameter * 0.5;
        let tangent_radius = x_tangency_point
            * (base_radius
                + (base_radius * base_radius + length * length
                    - x_tangency_point * x_tangency_point)
                    .sqrt())
            / (length + x_tangency_point);
        let nose_radius = (tangent_radius * tangent_radius + x_tangency_point * x_tangency_point)
            / (2.0 * x_tangency_point);
        let tangent_axial = length - x_tangency_point;
        let sphere_center_axial = length - nose_radius;

        for axial_index in 0..=length_subdivisions {
            let axial = length * axial_index as f32 / length_subdivisions as f32;
            let ring_radius = if axial <= tangent_axial {
                let slope = (base_radius - tangent_radius) / tangent_axial.max(f32::EPSILON);
                base_radius - slope * axial
            } else {
                let dx = axial - sphere_center_axial;
                (nose_radius * nose_radius - dx * dx).max(0.0).sqrt()
            };

            for angle_index in 0..angle_subdivisions {
                let angle =
                    2.0 * std::f32::consts::PI * angle_index as f32 / angle_subdivisions as f32;
                self.push_vertex(
                    ring_radius * angle.cos(),
                    -axial - self.y_offset,
                    ring_radius * angle.sin(),
                );
            }
        }

        for axial_index in 0..length_subdivisions {
            for angle_index in 0..angle_subdivisions {
                self.indices
                    .push(axial_index * angle_subdivisions + angle_index);
                self.indices.push(
                    axial_index * angle_subdivisions + (angle_index + 1) % angle_subdivisions,
                );
                self.indices.push(
                    (axial_index + 1) * angle_subdivisions + (angle_index + 1) % angle_subdivisions,
                );

                self.indices
                    .push(axial_index * angle_subdivisions + angle_index);
                self.indices.push(
                    (axial_index + 1) * angle_subdivisions + (angle_index + 1) % angle_subdivisions,
                );
                self.indices
                    .push((axial_index + 1) * angle_subdivisions + angle_index);
            }
        }
    }

    pub fn cone_blunted_ellipsoidal() {}

    pub fn cylinder(
        &mut self,
        diameter: f32,
        length: f32,
        length_subdivisions: u32,
        angle_subdivisions: u32,
    ) {
        let radius = diameter * 0.5;

        self.positions.clear();
        self.texture.clear();
        self.elements.clear();
        self.indices.clear();

        let sector_step = 2.0 * std::f32::consts::PI / angle_subdivisions as f32;
        let stack_step = 1.0 / length_subdivisions as f32;

        let y_top = -self.y_offset;
        let y_bot = -length - self.y_offset;

        for i in 0..=length_subdivisions {
            let t = i as f32 * stack_step;
            let y = -length * t - self.y_offset;

            for j in 0..angle_subdivisions {
                let a = j as f32 * sector_step;
                let x = radius * f32::cos(a);
                let z = radius * f32::sin(a);
                let s = j as f32 / angle_subdivisions as f32;

                self.positions.push(Vec3::new([x, y, z]));
                self.texture.push(Vec2::new([s, t]));
                self.elements
                    .push(Element::new(x, y, z, 1.0, 1.0, 1.0, 1.0, s, t));
            }
        }

        let row = angle_subdivisions;
        for i in 0..length_subdivisions {
            let k1 = i * row;
            let k2 = k1 + row;

            for j in 0..angle_subdivisions {
                let jn = (j + 1) % angle_subdivisions;

                let i0 = k1 + j;
                let i1 = k2 + j;
                let i2 = k1 + jn;
                let i3 = k2 + jn;

                self.indices.push(i0);
                self.indices.push(i1);
                self.indices.push(i2);

                self.indices.push(i2);
                self.indices.push(i1);
                self.indices.push(i3);
            }
        }

        let mut push_cap_vertex = |this: &mut Geometry, x: f32, y: f32, z: f32| -> u32 {
            let s = 0.5 + x / diameter;
            let t = 0.5 + z / diameter;
            this.positions.push(Vec3::new([x, y, z]));
            this.texture.push(Vec2::new([s, t]));
            this.elements
                .push(Element::new(x, y, z, 1.0, 1.0, 1.0, 1.0, s, t));
            (this.positions.len() - 1) as u32
        };

        let top_center = push_cap_vertex(self, 0.0, y_top, 0.0);
        let mut top_ring: Vec<u32> = Vec::with_capacity(angle_subdivisions as usize);
        for j in 0..angle_subdivisions {
            let a = j as f32 * sector_step;
            top_ring.push(push_cap_vertex(
                self,
                radius * a.cos(),
                y_top,
                radius * a.sin(),
            ));
        }
        for j in 0..angle_subdivisions {
            let jn = (j + 1) % angle_subdivisions;
            self.indices.push(top_center);
            self.indices.push(top_ring[j as usize]);
            self.indices.push(top_ring[jn as usize]);
        }

        let bot_center = push_cap_vertex(self, 0.0, y_bot, 0.0);
        let mut bot_ring: Vec<u32> = Vec::with_capacity(angle_subdivisions as usize);
        for j in 0..angle_subdivisions {
            let a = j as f32 * sector_step;
            bot_ring.push(push_cap_vertex(
                self,
                radius * a.cos(),
                y_bot,
                radius * a.sin(),
            ));
        }
        for j in 0..angle_subdivisions {
            let jn = (j + 1) % angle_subdivisions;
            self.indices.push(bot_center);
            self.indices.push(bot_ring[jn as usize]);
            self.indices.push(bot_ring[j as usize]);
        }
    }

    pub fn sphere(&mut self, radius: f32, stack_count: u32, sector_count: u32) {
        self.reset_mesh_buffers();

        if stack_count < 2 || sector_count < 3 {
            return;
        }

        let sector_step = 2.0 * std::f32::consts::PI / sector_count as f32;
        let stack_step = std::f32::consts::PI / stack_count as f32;
        let white = Vec4::new([1.0, 1.0, 1.0, 1.0]);

        for stack in 0..=stack_count {
            let polar_angle = stack as f32 * stack_step;
            let z = radius * polar_angle.cos();
            let xy_radius = radius * polar_angle.sin();
            let t = 1.0 - stack as f32 / stack_count as f32;

            for sector in 0..=sector_count {
                let azimuth = sector as f32 * sector_step;
                let x = xy_radius * azimuth.cos();
                let y = xy_radius * azimuth.sin();

                let s = if sector == sector_count {
                    1.0 - f32::EPSILON
                } else {
                    sector as f32 / sector_count as f32
                };

                self.positions.push(Vec3::new([x, y, z]));
                self.colors.push(white);
                self.texture.push(Vec2::new([s, t]));
                self.elements.push(Element::new(
                    x, y, z, white[0], white[1], white[2], white[3], s, t,
                ));
            }
        }

        for stack in 0..stack_count {
            let ring_start = stack * (sector_count + 1);
            let next_ring_start = ring_start + sector_count + 1;

            for sector in 0..sector_count {
                let i0 = ring_start + sector;
                let i1 = next_ring_start + sector;
                let i2 = i0 + 1;
                let i3 = i1 + 1;

                if stack != 0 {
                    self.indices.push(i0);
                    self.indices.push(i1);
                    self.indices.push(i2);
                }

                if stack != (stack_count - 1) {
                    self.indices.push(i2);
                    self.indices.push(i1);
                    self.indices.push(i3);
                }
            }
        }
    }

    fn circle(radius: f32, steps: usize) -> Vec<Vec3<f32>> {
        let mut points = vec![];
        if steps > 2 {
            let pi_2 = -1.0_f32.acos() * 2.0;
            for i in 1..steps {
                let a = pi_2 / (steps as f32 * i as f32);
                let x = radius * a.cos();
                let y = radius * a.sin();
                points.push(Vec3::new([x, y, 0.0]));
            }
        }
        points
    }

    pub fn cube(&mut self) {
        self.positions.clear();
        self.texture.clear();
        self.elements.clear();
        self.colors.clear();
        self.indices.clear();

        let c = Vec4::new([1.0, 1.0, 1.0, 1.0]);
        let n = -0.5f32;
        let p = 0.5f32;

        let ia = self.positions.len() as u32;
        self.positions.push(Vec3::new([p, n, n]));
        self.texture.push(Vec2::new([0.0, 0.0]));
        self.elements
            .push(Element::new(p, n, n, c[0], c[1], c[2], c[3], 0.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([p, n, p]));
        self.texture.push(Vec2::new([1.0, 0.0]));
        self.elements
            .push(Element::new(p, n, p, c[0], c[1], c[2], c[3], 1.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([p, p, p]));
        self.texture.push(Vec2::new([1.0, 1.0]));
        self.elements
            .push(Element::new(p, p, p, c[0], c[1], c[2], c[3], 1.0, 1.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([p, p, n]));
        self.texture.push(Vec2::new([0.0, 1.0]));
        self.elements
            .push(Element::new(p, p, n, c[0], c[1], c[2], c[3], 0.0, 1.0));
        self.colors.push(c);
        self.indices
            .extend_from_slice(&[ia, ia + 1, ia + 2, ia, ia + 2, ia + 3]);

        let ia = self.positions.len() as u32;
        self.positions.push(Vec3::new([n, n, p]));
        self.texture.push(Vec2::new([0.0, 0.0]));
        self.elements
            .push(Element::new(n, n, p, c[0], c[1], c[2], c[3], 0.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([n, n, n]));
        self.texture.push(Vec2::new([1.0, 0.0]));
        self.elements
            .push(Element::new(n, n, n, c[0], c[1], c[2], c[3], 1.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([n, p, n]));
        self.texture.push(Vec2::new([1.0, 1.0]));
        self.elements
            .push(Element::new(n, p, n, c[0], c[1], c[2], c[3], 1.0, 1.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([n, p, p]));
        self.texture.push(Vec2::new([0.0, 1.0]));
        self.elements
            .push(Element::new(n, p, p, c[0], c[1], c[2], c[3], 0.0, 1.0));
        self.colors.push(c);
        self.indices
            .extend_from_slice(&[ia, ia + 1, ia + 2, ia, ia + 2, ia + 3]);

        let ia = self.positions.len() as u32;
        self.positions.push(Vec3::new([n, p, n]));
        self.texture.push(Vec2::new([0.0, 0.0]));
        self.elements
            .push(Element::new(n, p, n, c[0], c[1], c[2], c[3], 0.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([p, p, n]));
        self.texture.push(Vec2::new([1.0, 0.0]));
        self.elements
            .push(Element::new(p, p, n, c[0], c[1], c[2], c[3], 1.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([p, p, p]));
        self.texture.push(Vec2::new([1.0, 1.0]));
        self.elements
            .push(Element::new(p, p, p, c[0], c[1], c[2], c[3], 1.0, 1.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([n, p, p]));
        self.texture.push(Vec2::new([0.0, 1.0]));
        self.elements
            .push(Element::new(n, p, p, c[0], c[1], c[2], c[3], 0.0, 1.0));
        self.colors.push(c);
        self.indices
            .extend_from_slice(&[ia, ia + 1, ia + 2, ia, ia + 2, ia + 3]);

        let ia = self.positions.len() as u32;
        self.positions.push(Vec3::new([n, n, p]));
        self.texture.push(Vec2::new([0.0, 0.0]));
        self.elements
            .push(Element::new(n, n, p, c[0], c[1], c[2], c[3], 0.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([p, n, p]));
        self.texture.push(Vec2::new([1.0, 0.0]));
        self.elements
            .push(Element::new(p, n, p, c[0], c[1], c[2], c[3], 1.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([p, n, n]));
        self.texture.push(Vec2::new([1.0, 1.0]));
        self.elements
            .push(Element::new(p, n, n, c[0], c[1], c[2], c[3], 1.0, 1.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([n, n, n]));
        self.texture.push(Vec2::new([0.0, 1.0]));
        self.elements
            .push(Element::new(n, n, n, c[0], c[1], c[2], c[3], 0.0, 1.0));
        self.colors.push(c);
        self.indices
            .extend_from_slice(&[ia, ia + 1, ia + 2, ia, ia + 2, ia + 3]);

        let ia = self.positions.len() as u32;
        self.positions.push(Vec3::new([n, n, p]));
        self.texture.push(Vec2::new([0.0, 0.0]));
        self.elements
            .push(Element::new(n, n, p, c[0], c[1], c[2], c[3], 0.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([n, p, p]));
        self.texture.push(Vec2::new([1.0, 0.0]));
        self.elements
            .push(Element::new(n, p, p, c[0], c[1], c[2], c[3], 1.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([p, p, p]));
        self.texture.push(Vec2::new([1.0, 1.0]));
        self.elements
            .push(Element::new(p, p, p, c[0], c[1], c[2], c[3], 1.0, 1.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([p, n, p]));
        self.texture.push(Vec2::new([0.0, 1.0]));
        self.elements
            .push(Element::new(p, n, p, c[0], c[1], c[2], c[3], 0.0, 1.0));
        self.colors.push(c);
        self.indices
            .extend_from_slice(&[ia, ia + 1, ia + 2, ia, ia + 2, ia + 3]);

        let ia = self.positions.len() as u32;
        self.positions.push(Vec3::new([p, n, n]));
        self.texture.push(Vec2::new([0.0, 0.0]));
        self.elements
            .push(Element::new(p, n, n, c[0], c[1], c[2], c[3], 0.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([p, p, n]));
        self.texture.push(Vec2::new([1.0, 0.0]));
        self.elements
            .push(Element::new(p, p, n, c[0], c[1], c[2], c[3], 1.0, 0.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([n, p, n]));
        self.texture.push(Vec2::new([1.0, 1.0]));
        self.elements
            .push(Element::new(n, p, n, c[0], c[1], c[2], c[3], 1.0, 1.0));
        self.colors.push(c);
        self.positions.push(Vec3::new([n, n, n]));
        self.texture.push(Vec2::new([0.0, 1.0]));
        self.elements
            .push(Element::new(n, n, n, c[0], c[1], c[2], c[3], 0.0, 1.0));
        self.colors.push(c);
        self.indices
            .extend_from_slice(&[ia, ia + 1, ia + 2, ia, ia + 2, ia + 3]);
    }

    pub fn triangular_prism(&mut self, a_side_length: f32, b_side_length: f32, depth: f32) {
        self.positions.clear();
        self.texture.clear();
        self.elements.clear();
        self.colors.clear();
        self.indices.clear();

        let c = Vec4::new([0.0, 1.0, 0.0, 1.0]);

        let v0 = Vec3::new([0.0, 0.0, 0.0]);
        let v1 = Vec3::new([a_side_length, 0.0, 0.0]);
        let v2 = Vec3::new([0.0, 0.0, b_side_length]);

        let v3 = Vec3::new([0.0, depth, 0.0]);
        let v4 = Vec3::new([a_side_length, depth, 0.0]);
        let v5 = Vec3::new([0.0, depth, b_side_length]);

        let verts = [v0, v1, v2, v3, v4, v5];
        for (i, v) in verts.iter().enumerate() {
            let s = if i % 2 == 0 { 0.0 } else { 1.0 };
            let t = if i < 3 { 0.0 } else { 1.0 };
            self.positions.push(*v);
            self.texture.push(Vec2::new([s, t]));
            self.elements
                .push(Element::new(v[0], v[1], v[2], c[0], c[1], c[2], c[3], s, t));
            self.colors.push(c);
        }

        self.indices.extend_from_slice(&[0, 1, 2]);
        self.indices.extend_from_slice(&[3, 5, 4]);
        self.indices.extend_from_slice(&[0, 1, 4, 0, 4, 3]);
        self.indices.extend_from_slice(&[1, 2, 5, 1, 5, 4]);
        self.indices.extend_from_slice(&[2, 0, 3, 2, 3, 5]);
    }

    pub fn ellipsoid() {}

    pub fn set_x_offset(&mut self, new_x_offset: f32) {
        self.x_offset = new_x_offset;
    }

    pub fn set_y_offset(&mut self, new_y_offset: f32) {
        self.y_offset = new_y_offset;
    }

    pub fn set_z_offset(&mut self, new_z_offset: f32) {
        self.z_offset = new_z_offset;
    }

    pub fn triangle(&mut self) {
        self.elements
            .push(Element::new(0.5, 0.5, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0));
        self.elements
            .push(Element::new(0.5, -0.5, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0));
        self.elements
            .push(Element::new(-0.5, -0.5, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0));
        self.elements
            .push(Element::new(-0.5, 0.5, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0));
        self.indices.push(0);
        self.indices.push(1);
        self.indices.push(1);
        self.indices.push(2);
        self.indices.push(3);
    }
}
