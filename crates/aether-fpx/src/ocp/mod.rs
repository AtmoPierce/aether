// Ours
// pub mod real;
// pub use real::*;
pub mod exp2i;
pub mod fp4;
pub mod fp6;
pub mod fp8;
pub mod rep;
pub mod scale;
pub mod trig;

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;

// Macros
#[macro_export]
macro_rules! fp8_e5 {
    ($lit:literal) => {
        $crate::ocp::fp8::Fp8E5M2::from($lit as f32)
    };
}

#[macro_export]
macro_rules! fp8_e4 {
    ($lit:literal) => {
        $crate::ocp::fp8::Fp8E4M3::from($lit as f32)
    };
}

// Macros
#[macro_export]
macro_rules! fp6_e3 {
    ($lit:literal) => {
        $crate::ocp::fp6::Fp6E3M2::from($lit as f32)
    };
}

#[macro_export]
macro_rules! fp4_e2 {
    ($lit:literal) => {
        $crate::ocp::fp4::Fp4E2M1::from($lit as f32)
    };
}
