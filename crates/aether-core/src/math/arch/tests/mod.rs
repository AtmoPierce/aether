#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod x86;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod matrix_simd_x86;

#[cfg(target_arch = "aarch64")]
mod neon;

#[cfg(target_arch = "aarch64")]
mod matrix_simd_arm;

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
mod precision;
