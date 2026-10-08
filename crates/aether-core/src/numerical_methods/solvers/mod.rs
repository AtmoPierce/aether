pub mod lu;
pub mod newton;
pub mod spd;
pub use lu::*;
pub use newton::*;
pub use spd::*;

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
