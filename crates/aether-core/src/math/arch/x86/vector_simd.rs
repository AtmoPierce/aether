use crate::math::Vector;

#[cfg(target_arch = "x86")]
use core::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

#[target_feature(enable = "avx,fma")]
pub unsafe fn dot4_avx_fma_f64(a: &Vector<f64, 4>, b: &Vector<f64, 4>) -> f64 {
    let av = unsafe { _mm256_loadu_pd(a.data.as_ptr()) };
    let bv = unsafe { _mm256_loadu_pd(b.data.as_ptr()) };
    let p = _mm256_fmadd_pd(av, bv, _mm256_setzero_pd());
    let mut lanes = [0.0_f64; 4];
    unsafe { _mm256_storeu_pd(lanes.as_mut_ptr(), p) };
    lanes[0] + lanes[1] + lanes[2] + lanes[3]
}

#[target_feature(enable = "avx")]
pub unsafe fn dot4_avx_f64(a: &Vector<f64, 4>, b: &Vector<f64, 4>) -> f64 {
    let av = unsafe { _mm256_loadu_pd(a.data.as_ptr()) };
    let bv = unsafe { _mm256_loadu_pd(b.data.as_ptr()) };
    let p = _mm256_mul_pd(av, bv);
    let mut lanes = [0.0_f64; 4];
    unsafe { _mm256_storeu_pd(lanes.as_mut_ptr(), p) };
    lanes[0] + lanes[1] + lanes[2] + lanes[3]
}

#[target_feature(enable = "sse,fma")]
pub unsafe fn dot4_sse_fma_f32(a: &Vector<f32, 4>, b: &Vector<f32, 4>) -> f32 {
    let av = unsafe { _mm_loadu_ps(a.data.as_ptr()) };
    let bv = unsafe { _mm_loadu_ps(b.data.as_ptr()) };
    let p = _mm_fmadd_ps(av, bv, _mm_setzero_ps());
    let mut lanes = [0.0_f32; 4];
    unsafe { _mm_storeu_ps(lanes.as_mut_ptr(), p) };
    lanes[0] + lanes[1] + lanes[2] + lanes[3]
}

#[target_feature(enable = "sse")]
pub unsafe fn dot4_sse_f32(a: &Vector<f32, 4>, b: &Vector<f32, 4>) -> f32 {
    let av = unsafe { _mm_loadu_ps(a.data.as_ptr()) };
    let bv = unsafe { _mm_loadu_ps(b.data.as_ptr()) };
    let p = _mm_mul_ps(av, bv);
    let mut lanes = [0.0_f32; 4];
    unsafe { _mm_storeu_ps(lanes.as_mut_ptr(), p) };
    lanes[0] + lanes[1] + lanes[2] + lanes[3]
}
