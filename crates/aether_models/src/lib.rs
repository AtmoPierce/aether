// Re-export core modules so `crate::{math, coordinate, ...}` resolve
pub use aether::{attitude, coordinate, math, numerical_methods, reference_frame, utils};

// Expose this crate's content under `crate::models::...` as expected by modules
pub mod models;

// Root modules for this crate's own content
pub mod celestial;
pub mod terrestial;
