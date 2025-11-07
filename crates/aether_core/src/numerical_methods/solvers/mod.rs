pub mod lu;
pub mod spd;
pub mod newton;
pub use lu::*;
pub use spd::*;
pub use newton::*;  

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;