//! Remez export + codegen helpers.
pub mod accuracy;
pub mod coefficients;
pub mod config;
pub mod export;
pub mod generators;
pub mod macros;
pub mod metrics;
mod payne_hanek;
pub mod plots;
pub mod report;
pub use macros::*;

pub mod analysis;
pub use analysis::*;

pub use metrics::{
    ABS_TOL_NEAR_ZERO_F32, ABS_TOL_NEAR_ZERO_F64, MAX_ULP_ALLOWED_F128, MAX_ULP_ALLOWED_F32,
    MAX_ULP_ALLOWED_F64,
};

pub use export::{write_coeffs_to_c, write_function_to_rs, RemezExport};

pub use config::load_fpx_remez_toml;
