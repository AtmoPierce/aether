pub mod macros;
pub mod matrix;
pub mod vector;
pub use matrix::Matrix;
pub use vector::Vector;

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
