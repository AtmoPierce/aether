pub mod macros;
pub mod matrix;
pub mod vector;
pub mod tensor;
pub use matrix::Matrix;
pub use vector::Vector;
pub use tensor::{Tensor3, Tensor4, TensorView, TensorViewMut};

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
