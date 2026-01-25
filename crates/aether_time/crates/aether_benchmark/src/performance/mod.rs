mod clock;
mod performance;
mod power;
mod timer;
pub use clock::*;
#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
