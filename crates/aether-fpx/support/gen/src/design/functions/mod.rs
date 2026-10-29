pub mod sqrt_f;
pub use sqrt_f::*;

// sqrt_generate is not useful to approximate...
// pub mod sqrt_generate;

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
