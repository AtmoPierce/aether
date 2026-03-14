use rkyv::{Archive, Deserialize as RkyvDeserialize, Serialize as RkyvSerialize};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Read;
use sha2::{Digest, Sha256};
pub mod color;
use color::gaia_color::{bp_rp_to_rgb, srgb_to_linear};

/// GNC Catalog record from USNO GNC v1.1
#[derive(Archive, RkyvDeserialize, RkyvSerialize, Deserialize, Serialize, Debug)]
pub struct GncRecord {
    pub gnc_id: u64,
    pub hip_id: Option<u64>,
    pub gaiadr3_id: Option<u64>,

    pub ra_2016: f64,
    pub de_2016: f64,
    pub e_ra_2016: f32,
    pub e_de_2016: f32,
    pub pmra: Option<f32>,
    pub pmde: Option<f32>,
    pub e_pmra: Option<f32>,
    pub e_pmde: Option<f32>,
    pub plx: Option<f32>,
    pub e_plx: Option<f32>,
    pub cat_2016: Option<String>,

    pub gmag: Option<f32>,
    pub e_gmag: Option<f32>,
    pub bpmag: Option<f32>,
    pub e_bpmag: Option<f32>,
    pub rpmag: Option<f32>,
    pub e_rpmag: Option<f32>,

    pub jmag: Option<f32>,
    pub e_jmag: Option<f32>,
    pub hmag: Option<f32>,
    pub e_hmag: Option<f32>,
    pub kmag: Option<f32>,
    pub e_kmag: Option<f32>,

    pub flag_var: Option<String>,
    pub flag_mult: Option<String>,

    pub neighbor_id: u64,
    pub neighbor_distance_arcsec: f32,
    pub shift_arsec: Option<f32>,

    pub flag_shift: Option<String>,
    pub flag_pos: Option<String>,
    pub flag_pm: Option<String>,
    pub flag_neighbor: Option<String>,
}

pub struct GncCatalogReader {
    pub records: Vec<GncRecord>,
}

impl GncCatalogReader {
    pub const GNC_V1_1_MAR_2_2023_SHA256: &'static str = "d317268726ff5fd28d87d6e3d0a3aa0f26abc0408a43e711f4b87579bb00b640";

    pub fn default_catalog_path() -> std::path::PathBuf {
        std::path::PathBuf::from(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/catalog/gnc_v1_1_mar_2_2023.csv"
        ))
    }

    /// Load the default catalog. Lookup order:
    /// 1. `CARGO_MANIFEST_DIR`/catalog/... (compile-time layout)
    /// 2. `$AETHER_CATALOGS_DIR/gnc_v1_1_mar_2_2023.csv` (runtime override)
    /// If neither exists, return an error instructing the user to download the CSV.
    pub fn from_default_catalog() -> Result<Self, Box<dyn std::error::Error>> {
        let local = Self::default_catalog_path();
        if local.exists() {
            return Self::from_csv_verified_v1_1(local);
        }

        if let Some(env_dir) = std::env::var_os("AETHER_CATALOGS_DIR") {
            let p = std::path::PathBuf::from(env_dir).join("gnc_v1_1_mar_2_2023.csv");
            if p.exists() {
                return Self::from_csv_verified_v1_1(p);
            }
        }

        Err(format!(
            "Catalog not found at {}. Please download the CSV and place it in that path or set `AETHER_CATALOGS_DIR`.",
            local.display()
        )
        .into())
    }

    

    pub fn from_csv_verified_v1_1<P: AsRef<std::path::Path>>(
        path: P,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        Self::from_csv_with_expected_sha256(path, Self::GNC_V1_1_MAR_2_2023_SHA256)
    }

    pub fn from_csv_with_expected_sha256<P: AsRef<std::path::Path>>(
        path: P,
        expected_sha256_hex: &str,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        Self::verify_catalog_checksum(&path, expected_sha256_hex)?;
        Self::from_csv(path)
    }

    pub fn verify_catalog_checksum<P: AsRef<std::path::Path>>(
        path: P,
        expected_sha256_hex: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::open(path.as_ref())?;
        let mut hasher = Sha256::new();
        let mut buffer = [0_u8; 64 * 1024];

        loop {
            let read = file.read(&mut buffer)?;
            if read == 0 {
                break;
            }
            hasher.update(&buffer[..read]);
        }

        let digest_hex = format!("{:x}", hasher.finalize());
        if digest_hex.eq_ignore_ascii_case(expected_sha256_hex) {
            Ok(())
        } else {
            Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Catalog SHA-256 checksum mismatch for {}: expected {}, got {}",
                    path.as_ref().display(),
                    expected_sha256_hex,
                    digest_hex,
                ),
            )
            .into())
        }
    }

    fn propagated_uncertainty(
        t: f64,
        e_ra_2016: f32,
        e_pmra: Option<f32>,
        e_de_2016: f32,
        e_pmde: Option<f32>,
        cat_2016: &Option<String>,
    ) -> (f64, f64) {
        let e_ra_2016 = e_ra_2016 as f64;
        let e_de_2016 = e_de_2016 as f64;
        let e_pmra = e_pmra.unwrap_or(0.0) as f64;
        let e_pmde = e_pmde.unwrap_or(0.0) as f64;

        let (e_ra_t, e_de_t) = match cat_2016 {
            Some(c) if c == "GAIADR3" => {
                let e_ra = (e_ra_2016.powi(2) + e_pmra.powi(2) * (t - 2016.0_f64).powi(2)).sqrt();
                let e_de = (e_de_2016.powi(2) + e_pmde.powi(2) * (t - 2016.0_f64).powi(2)).sqrt();
                (e_ra, e_de)
            }
            _ => {
                let delta_t2 = (t - 1991.25_f64).powi(2) - (2016.0_f64 - 1991.25_f64).powi(2);
                let e_ra = (e_ra_2016.powi(2) + e_pmra.powi(2) * delta_t2).sqrt();
                let e_de = (e_de_2016.powi(2) + e_pmde.powi(2) * delta_t2).sqrt();
                (e_ra, e_de)
            }
        };

        (e_ra_t, e_de_t)
    }

    pub fn from_csv<P: AsRef<std::path::Path>>(
        path: P,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(path)?;
        let mut rdr = csv::ReaderBuilder::new()
            .flexible(true)
            .from_reader(file);
        let mut records = Vec::new();
        let mut skipped = 0usize;

        for result in rdr.deserialize() {
            match result {
                Ok(record) => records.push(record),
                Err(_) => skipped += 1,
            }
        }

        if skipped > 0 {
            eprintln!(
                "[aether_catalogs] Skipped {} malformed GNC catalog rows",
                skipped
            );
        }

        Ok(GncCatalogReader { records })
    }
}

#[derive(Debug, Clone)]
pub struct GncCatalogMotionSample {
    pub epoch: f64,
    pub icrs_position_au: [f64; 3],
    pub icrs_velocity_au_per_yr: [f64; 3],
    pub uncertainty_ra_mas: f64,
    pub uncertainty_dec_mas: f64,
}

pub trait IcrsPropagator {
    fn icrs_sample(&self, epoch: f64) -> GncCatalogMotionSample;
}

impl IcrsPropagator for GncRecord {
    fn icrs_sample(&self, epoch: f64) -> GncCatalogMotionSample {
        let dt = epoch - 2016.0;

        let ra_deg = self.ra_2016 + (self.pmra.unwrap_or(0.0) as f64 * dt / 3_600_000.0);
        let dec_deg = self.de_2016 + (self.pmde.unwrap_or(0.0) as f64 * dt / 3_600_000.0);

        let ra_rad = ra_deg.to_radians();
        let dec_rad = dec_deg.to_radians();

        let cos_dec = dec_rad.cos();
        let sin_dec = dec_rad.sin();
        let cos_ra = ra_rad.cos();
        let sin_ra = ra_rad.sin();

        let unit = [cos_dec * cos_ra, cos_dec * sin_ra, sin_dec];

        let distance_au = self
            .plx
            .filter(|&p| p > 0.0)
            .map(|plx_mas| {
                let plx_arcsec = plx_mas as f64 / 1000.0;
                let distance_pc = 1.0 / plx_arcsec;
                distance_pc * 206_264.806
            })
            .unwrap_or(1.0); // fallback: 1 AU if no parallax or invalid

        let position = unit.map(|v| v * distance_au);

        let velocity = match (self.plx, self.pmra, self.pmde) {
            (Some(plx), Some(pmra), Some(pmde)) if plx > 0.0 => {
                let plx_mas = plx as f64;
                let scale = 4.74047 / plx_mas;
                let mu_ra = pmra as f64;
                let mu_de = pmde as f64;

                let t_ra = [-sin_ra, cos_ra, 0.0];
                let t_de = [-cos_ra * sin_dec, -sin_ra * sin_dec, cos_dec];

                [
                    scale * (mu_ra * t_ra[0] + mu_de * t_de[0]),
                    scale * (mu_ra * t_ra[1] + mu_de * t_de[1]),
                    scale * (mu_ra * t_ra[2] + mu_de * t_de[2]),
                ]
            }
            _ => [0.0, 0.0, 0.0],
        };

        let (e_ra, e_dec) = GncCatalogReader::propagated_uncertainty(
            epoch,
            self.e_ra_2016,
            self.e_pmra,
            self.e_de_2016,
            self.e_pmde,
            &self.cat_2016,
        );

        GncCatalogMotionSample {
            epoch,
            icrs_position_au: position,
            icrs_velocity_au_per_yr: velocity,
            uncertainty_ra_mas: e_ra,
            uncertainty_dec_mas: e_dec,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GncCatalogPhotoSample {
    pub gmag: f32,
    pub bpmag: f32,
    pub rpmag: f32,
    pub rgb: [f32; 3],
}

pub trait PhotometricColor {
    fn photometric_color(&self) -> Option<GncCatalogPhotoSample>;
    fn random_color(&self) -> Option<GncCatalogPhotoSample>;
}

impl PhotometricColor for GncRecord {
    fn photometric_color(&self) -> Option<GncCatalogPhotoSample> {
        // Use explicit pattern match to avoid any accidental early-return None
        if let (Some(bp), Some(rp), Some(g)) = (self.bpmag, self.rpmag, self.gmag) {
            let bp_rp = bp - rp;
            let rgb_srgb = bp_rp_to_rgb(bp_rp, false);
            // If your renderer expects linear (most do), convert:
            let rgb = srgb_to_linear(rgb_srgb);
            return Some(GncCatalogPhotoSample {
                gmag: g,
                bpmag: bp,
                rpmag: rp,
                rgb,
            });
        }
        None
    }

    fn random_color(&self) -> Option<GncCatalogPhotoSample> {
        use rand::Rng;
        let mut rng = rand::rng();
        let r: f32 = rng.random_range(0.0..1.0);
        let g: f32 = rng.random_range(0.0..0.1); // Stars are more red/bluish than greenish (imo)
        let b: f32 = rng.random_range(0.0..1.0);
        Some(GncCatalogPhotoSample {
            gmag: 0.0,
            bpmag: 0.0,
            rpmag: 0.0,
            rgb: [r, g, b],
        })
    }
}

#[derive(Debug, Clone, Copy)]
pub enum GncCatalogInfraredColor {
    From2Mass {
        rgb: [f32; 3],
        j: f32,
        h: f32,
        k: f32,
    },
    Missing,
}

pub trait InfraredColor {
    fn photometric_infrared_rgb(&self) -> Option<GncCatalogInfraredColor>;
}

impl InfraredColor for GncRecord {
    fn photometric_infrared_rgb(&self) -> Option<GncCatalogInfraredColor> {
        match (self.jmag, self.hmag, self.kmag) {
            (Some(j), Some(h), Some(k)) => {
                let clamp = |mag: f32| (1.0 - ((mag - 8.0) / 8.0)).clamp(0.0, 1.0);
                let r = clamp(k);
                let g = clamp(h);
                let b = clamp(j);

                Some(GncCatalogInfraredColor::From2Mass {
                    rgb: [r, g, b],
                    j,
                    h,
                    k,
                })
            }
            _ => Some(GncCatalogInfraredColor::Missing),
        }
    }
}
