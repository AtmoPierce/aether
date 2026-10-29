//! Explicit architecture-specific SIMD polynomial evaluation.
//!
//! Each function evaluates one polynomial at every input value in place. The
//! slice API is shared by all targets; only the packed lane width differs.

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod x86;

#[cfg(target_arch = "aarch64")]
mod aarch64;

/// Evaluate one polynomial for one binary32 input by placing distinct
/// coefficient streams in SIMD lanes and combining them horizontally.
#[inline]
pub fn eval_blocked_horner_f32<const N: usize>(x: f32, coefficients: &[f32; N]) -> f32 {
    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        not(target_os = "none")
    ))]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: both target features were detected immediately above.
        return unsafe { x86::eval_blocked_horner_f32(x, coefficients) };
    }
    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is part of the AArch64 architecture profile.
        return unsafe { aarch64::eval_blocked_horner_f32(x, coefficients) };
    }
    #[allow(unreachable_code)]
    scalar_horner_value_f32(x, coefficients)
}

/// Evaluate one polynomial for one binary32 input with a horizontal SIMD
/// Estrin tree.
#[inline]
pub fn eval_horizontal_estrin_f32<const N: usize>(x: f32, coefficients: &[f32; N]) -> f32 {
    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        not(target_os = "none")
    ))]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: both target features were detected immediately above.
        return unsafe { x86::eval_horizontal_estrin_f32(x, coefficients) };
    }
    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is part of the AArch64 architecture profile.
        return unsafe { aarch64::eval_horizontal_estrin_f32(x, coefficients) };
    }
    #[allow(unreachable_code)]
    crate::trig::common::eval_estrin_f32_slice(x, coefficients)
}

/// Binary64 counterpart to [`eval_blocked_horner_f32`].
#[inline]
pub fn eval_blocked_horner_f64<const N: usize>(x: f64, coefficients: &[f64; N]) -> f64 {
    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        not(target_os = "none")
    ))]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: both target features were detected immediately above.
        return unsafe { x86::eval_blocked_horner_f64(x, coefficients) };
    }
    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is part of the AArch64 architecture profile.
        return unsafe { aarch64::eval_blocked_horner_f64(x, coefficients) };
    }
    #[allow(unreachable_code)]
    scalar_horner_value_f64(x, coefficients)
}

/// Binary64 counterpart to [`eval_horizontal_estrin_f32`].
#[inline]
pub fn eval_horizontal_estrin_f64<const N: usize>(x: f64, coefficients: &[f64; N]) -> f64 {
    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        not(target_os = "none")
    ))]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: both target features were detected immediately above.
        return unsafe { x86::eval_horizontal_estrin_f64(x, coefficients) };
    }
    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is part of the AArch64 architecture profile.
        return unsafe { aarch64::eval_horizontal_estrin_f64(x, coefficients) };
    }
    #[allow(unreachable_code)]
    crate::trig::common::eval_estrin_f64_slice(x, coefficients)
}

#[inline]
pub fn eval_horner_f32_in_place<const N: usize>(values: &mut [f32], coefficients: &[f32; N]) {
    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        not(target_os = "none")
    ))]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: both target features were detected immediately above.
        unsafe { x86::eval_horner_f32(values, coefficients) };
        return;
    }

    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        target_os = "none",
        target_feature = "avx2",
        target_feature = "fma"
    ))]
    {
        // SAFETY: the required features are part of the compilation target.
        unsafe { x86::eval_horner_f32(values, coefficients) };
        return;
    }

    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is part of the AArch64 architecture profile.
        unsafe { aarch64::eval_horner_f32(values, coefficients) };
        return;
    }

    #[allow(unreachable_code)]
    scalar_horner_f32(values, coefficients);
}

#[inline]
pub fn eval_estrin_f32_in_place<const N: usize>(values: &mut [f32], coefficients: &[f32; N]) {
    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        not(target_os = "none")
    ))]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: both target features were detected immediately above.
        unsafe { x86::eval_estrin_f32(values, coefficients) };
        return;
    }

    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        target_os = "none",
        target_feature = "avx2",
        target_feature = "fma"
    ))]
    {
        // SAFETY: the required features are part of the compilation target.
        unsafe { x86::eval_estrin_f32(values, coefficients) };
        return;
    }

    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is part of the AArch64 architecture profile.
        unsafe { aarch64::eval_estrin_f32(values, coefficients) };
        return;
    }

    #[allow(unreachable_code)]
    scalar_estrin_f32(values, coefficients);
}

#[inline]
pub fn eval_horner_f64_in_place<const N: usize>(values: &mut [f64], coefficients: &[f64; N]) {
    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        not(target_os = "none")
    ))]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: both target features were detected immediately above.
        unsafe { x86::eval_horner_f64(values, coefficients) };
        return;
    }

    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        target_os = "none",
        target_feature = "avx2",
        target_feature = "fma"
    ))]
    {
        // SAFETY: the required features are part of the compilation target.
        unsafe { x86::eval_horner_f64(values, coefficients) };
        return;
    }

    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is part of the AArch64 architecture profile.
        unsafe { aarch64::eval_horner_f64(values, coefficients) };
        return;
    }

    #[allow(unreachable_code)]
    scalar_horner_f64(values, coefficients);
}

#[inline]
pub fn eval_estrin_f64_in_place<const N: usize>(values: &mut [f64], coefficients: &[f64; N]) {
    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        not(target_os = "none")
    ))]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: both target features were detected immediately above.
        unsafe { x86::eval_estrin_f64(values, coefficients) };
        return;
    }

    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        target_os = "none",
        target_feature = "avx2",
        target_feature = "fma"
    ))]
    {
        // SAFETY: the required features are part of the compilation target.
        unsafe { x86::eval_estrin_f64(values, coefficients) };
        return;
    }

    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is part of the AArch64 architecture profile.
        unsafe { aarch64::eval_estrin_f64(values, coefficients) };
        return;
    }

    #[allow(unreachable_code)]
    scalar_estrin_f64(values, coefficients);
}

fn scalar_horner_value_f32<const N: usize>(x: f32, coefficients: &[f32; N]) -> f32 {
    let mut result = 0.0_f32;
    for &coefficient in coefficients.iter().rev() {
        result = result.mul_add(x, coefficient);
    }
    result
}

fn scalar_horner_f32<const N: usize>(values: &mut [f32], coefficients: &[f32; N]) {
    for value in values {
        *value = scalar_horner_value_f32(*value, coefficients);
    }
}

fn scalar_estrin_f32<const N: usize>(values: &mut [f32], coefficients: &[f32; N]) {
    for value in values {
        *value = crate::trig::common::eval_estrin_f32_slice(*value, coefficients);
    }
}

fn scalar_horner_value_f64<const N: usize>(x: f64, coefficients: &[f64; N]) -> f64 {
    let mut result = 0.0_f64;
    for &coefficient in coefficients.iter().rev() {
        result = result.mul_add(x, coefficient);
    }
    result
}

fn scalar_horner_f64<const N: usize>(values: &mut [f64], coefficients: &[f64; N]) {
    for value in values {
        *value = scalar_horner_value_f64(*value, coefficients);
    }
}

fn scalar_estrin_f64<const N: usize>(values: &mut [f64], coefficients: &[f64; N]) {
    for value in values {
        *value = crate::trig::common::eval_estrin_f64_slice(*value, coefficients);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const COEFFICIENTS_F32: [f32; 7] = [1.0, -0.5, 0.25, -0.125, 0.0625, -0.03125, 0.015625];
    const COEFFICIENTS_F64: [f64; 7] = [1.0, -0.5, 0.25, -0.125, 0.0625, -0.03125, 0.015625];

    #[test]
    fn f32_batch_matches_scalar_evaluators_including_tail() {
        let input = [-0.95, -0.7, -0.3, -0.1, 0.0, 0.2, 0.4, 0.65, 0.9, 1.1, 1.4];
        let mut horner = input;
        let mut estrin = input;
        eval_horner_f32_in_place(&mut horner, &COEFFICIENTS_F32);
        eval_estrin_f32_in_place(&mut estrin, &COEFFICIENTS_F32);

        for (index, x) in input.into_iter().enumerate() {
            let expected_horner = crate::trig::common::eval_horner_f32_slice(x, &COEFFICIENTS_F32);
            let expected_estrin = crate::trig::common::eval_estrin_f32_slice(x, &COEFFICIENTS_F32);
            assert!((horner[index] - expected_horner).abs() <= 2.0 * f32::EPSILON);
            assert!((estrin[index] - expected_estrin).abs() <= 2.0 * f32::EPSILON);
        }
    }

    #[test]
    fn f64_batch_matches_scalar_evaluators_including_tail() {
        let input = [-0.95, -0.7, -0.3, -0.1, 0.0, 0.2, 0.4, 0.65, 0.9];
        let mut horner = input;
        let mut estrin = input;
        eval_horner_f64_in_place(&mut horner, &COEFFICIENTS_F64);
        eval_estrin_f64_in_place(&mut estrin, &COEFFICIENTS_F64);

        for (index, x) in input.into_iter().enumerate() {
            let expected_horner = crate::trig::common::eval_horner_f64_slice(x, &COEFFICIENTS_F64);
            let expected_estrin = crate::trig::common::eval_estrin_f64_slice(x, &COEFFICIENTS_F64);
            assert!((horner[index] - expected_horner).abs() <= 2.0 * f64::EPSILON);
            assert!((estrin[index] - expected_estrin).abs() <= 2.0 * f64::EPSILON);
        }
    }
}
