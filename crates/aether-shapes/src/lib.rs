#![cfg_attr(all(feature = "no_std", not(feature = "std")), no_std)]

pub mod attributes;
pub mod cylinder;
pub mod prism;
pub mod sphere;

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
