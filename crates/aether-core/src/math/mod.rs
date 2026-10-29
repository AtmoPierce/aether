pub mod algorithms;
pub mod arch;
pub mod complex;
pub mod macros;
pub mod matrix;
pub mod tensor;
pub mod vector;
pub use complex::{Complex, ComplexField};
pub use matrix::Matrix;
pub use tensor::{Tensor3, Tensor4, TensorView, TensorViewMut};
pub use vector::Vector;

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
