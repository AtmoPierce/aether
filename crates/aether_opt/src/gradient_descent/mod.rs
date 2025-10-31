pub mod basic;
pub use basic::*;

pub mod stochastic;
pub use stochastic::*;

#[cfg(std)]
pub mod stochastic_parallel;
#[cfg(std)]
pub use stochastic_parallel::*;

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
