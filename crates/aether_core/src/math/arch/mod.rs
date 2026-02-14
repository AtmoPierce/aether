#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub mod x86;

#[cfg(any(target_arch = "arm", target_arch = "aarch64"))]
pub mod arm;

#[cfg(test)]
pub mod tests;
