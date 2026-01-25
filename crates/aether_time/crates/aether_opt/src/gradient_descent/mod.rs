pub mod basic;
pub use basic::*;

pub mod stochastic;
pub use stochastic::*;

pub mod stochastic_parallel;
pub use stochastic_parallel::*;

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
