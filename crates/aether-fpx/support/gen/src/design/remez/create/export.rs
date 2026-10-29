//! Coefficient export helpers (Rust/C/JSON) for aether-fpx-gen.

use crate::core::fp4::Fp4E2M1;
use crate::core::fp6::Fp6E3M2;
use crate::core::fp8::{Fp8E4M3, Fp8E5M2};
use crate::remez::RemezApprox;
use core::fmt::Debug;
use core::{f128, f16};
use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use super::analysis::solve_interval_from_x_for_name;
use super::coefficients::SegmentPlan;
use super::coefficients::{
    poly_domain, CoeffSet128, CoeffSet16, CoeffSet32, CoeffSet64, FunctionPlan, PolyDomain,
    TrigKernel,
};
use super::config::{load_fpx_remez_toml, FunctionSpec};

mod acos;
mod asin;
mod atan;
mod common;
mod cos;
mod cosh;
mod exp;
mod expm1;
mod log;
mod log1p;
mod sin;
mod sinh;
mod tan;
mod tanh;

/* ------------------------------------------------------------------------------------------------
 * Utilities
 * ------------------------------------------------------------------------------------------------ */

fn out_dir_for(kind: &str, default_rel: &str) -> PathBuf {
    let cfg = load_fpx_remez_toml().expect("failed to load fpx-remez.toml");
    match kind {
        "rs" => cfg.output.rs_out_dir,
        "json" => cfg.output.json_out_dir,
        "c" => cfg.output.c_out_dir,
        "gpu" => cfg.output.gpu_out_dir,
        "plots" => cfg.output.plots_out_dir,
        _ => PathBuf::from(default_rel),
    }
}

fn sanitize_ident(name: &str) -> String {
    let mut out = String::new();
    for c in name.chars() {
        if c.is_ascii_alphanumeric() {
            out.push(c.to_ascii_lowercase());
        } else {
            out.push('_');
        }
    }
    if out.is_empty() {
        "fpx".into()
    } else {
        out
    }
}

fn quantize_coeffs_fp4_bits(coeffs: &[f32]) -> Vec<u8> {
    coeffs.iter().map(|&c| Fp4E2M1::from(c).to_bits()).collect()
}

fn quantize_coeffs_fp6_bits(coeffs: &[f32]) -> Vec<u8> {
    coeffs.iter().map(|&c| Fp6E3M2::from(c).to_bits()).collect()
}

fn quantize_coeffs_fp8_e4m3_bits(coeffs: &[f32]) -> Vec<u8> {
    coeffs.iter().map(|&c| Fp8E4M3::from(c).to_bits()).collect()
}

fn quantize_coeffs_fp8_e5m2_bits(coeffs: &[f32]) -> Vec<u8> {
    coeffs.iter().map(|&c| Fp8E5M2::from(c).to_bits()).collect()
}

fn coeff_decl_f32(symbol: &str, coeffs: &[f32], backend: &str) -> String {
    let qualifier = match backend {
        "cuda" | "rocm" => "__device__ __constant__",
        _ => "const",
    };
    let elems = coeffs
        .iter()
        .map(|c| format!("{:.9e}f", c))
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        "{qualifier} float {symbol}[{}] = {{ {elems} }};",
        coeffs.len()
    )
}

fn coeff_decl_f64(symbol: &str, coeffs: &[f64], backend: &str) -> String {
    let qualifier = match backend {
        "cuda" | "rocm" => "__device__ __constant__",
        _ => "const",
    };
    let elems = coeffs
        .iter()
        .map(|c| format!("{:.17e}", c))
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        "{qualifier} double {symbol}[{}] = {{ {elems} }};",
        coeffs.len()
    )
}

fn coeff_decl_u16(symbol: &str, coeffs: &[u16], backend: &str) -> String {
    let qualifier = match backend {
        "cuda" | "rocm" => "__device__ __constant__",
        _ => "const",
    };
    let elems = coeffs
        .iter()
        .map(|c| format!("0x{c:04x}"))
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        "{qualifier} unsigned short {symbol}[{}] = {{ {elems} }};",
        coeffs.len()
    )
}

fn coeff_decl_u8(symbol: &str, coeffs: &[u8], backend: &str) -> String {
    let qualifier = match backend {
        "cuda" | "rocm" => "__device__ __constant__",
        _ => "const",
    };
    let elems = coeffs
        .iter()
        .map(|c| format!("0x{c:02x}"))
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        "{qualifier} unsigned char {symbol}[{}] = {{ {elems} }};",
        coeffs.len()
    )
}

fn gpu_low_precision_prelude(backend: &str) -> String {
    if backend == "cuda" {
        return r#"using fpx_half = __half;
using fpx_fp4_e2m1 = __nv_fp4_e2m1;
using fpx_fp6_e3m2 = __nv_fp6_e3m2;
using fpx_fp8_e4m3 = __nv_fp8_e4m3;
using fpx_fp8_e5m2 = __nv_fp8_e5m2;

#ifndef FPX_CUDA_USE_FMA
#define FPX_CUDA_USE_FMA 0
#endif

static_assert(sizeof(fpx_half) == sizeof(unsigned short), "expected 16-bit __half");
static_assert(sizeof(fpx_fp4_e2m1) == sizeof(unsigned char), "expected 8-bit __nv_fp4_e2m1 storage");
static_assert(sizeof(fpx_fp6_e3m2) == sizeof(unsigned char), "expected 8-bit __nv_fp6_e3m2 storage");
static_assert(sizeof(fpx_fp8_e4m3) == sizeof(unsigned char), "expected 8-bit __nv_fp8_e4m3 storage");
static_assert(sizeof(fpx_fp8_e5m2) == sizeof(unsigned char), "expected 8-bit __nv_fp8_e5m2 storage");

__device__ __forceinline__ fpx_half fpx_hfma_fallback(fpx_half a, fpx_half b, fpx_half c) {
    return static_cast<fpx_half>(static_cast<float>(a) * static_cast<float>(b) + static_cast<float>(c));
}

#if !defined(__CUDA_ARCH__) || (__CUDA_ARCH__ < 530)
#define __hfma(a, b, c) fpx_hfma_fallback((a), (b), (c))
#endif

typedef union {
    unsigned short u;
    fpx_half h;
} fpx_half_bits;

typedef union {
    unsigned char u;
    fpx_fp4_e2m1 v;
} fpx_fp4_e2m1_bits;

typedef union {
    unsigned char u;
    fpx_fp6_e3m2 v;
} fpx_fp6_e3m2_bits;

typedef union {
    unsigned char u;
    fpx_fp8_e4m3 v;
} fpx_fp8_e4m3_bits;

typedef union {
    unsigned char u;
    fpx_fp8_e5m2 v;
} fpx_fp8_e5m2_bits;

__device__ __forceinline__ float fpx_mul_add_f32(float a, float b, float c) {
#if FPX_CUDA_USE_FMA
    return fmaf(a, b, c);
#else
    return a * b + c;
#endif
}

__device__ __forceinline__ double fpx_mul_add_f64(double a, double b, double c) {
#if FPX_CUDA_USE_FMA
    return fma(a, b, c);
#else
    return a * b + c;
#endif
}

__device__ __forceinline__ fpx_half fpx_f32_to_half(float x) { return static_cast<fpx_half>(x); }
__device__ __forceinline__ float fpx_half_to_f32(fpx_half x) { return static_cast<float>(x); }
__device__ __forceinline__ fpx_half fpx_decode_f16_bits(unsigned short bits) { fpx_half_bits v; v.u = bits; return v.h; }

__device__ __forceinline__ fpx_fp4_e2m1 fpx_f32_to_fp4_e2m1(float x) { return static_cast<fpx_fp4_e2m1>(x); }
__device__ __forceinline__ float fpx_fp4_e2m1_to_f32(fpx_fp4_e2m1 x) { return static_cast<float>(x); }
__device__ __forceinline__ fpx_fp4_e2m1 fpx_decode_fp4_e2m1(unsigned char bits) { fpx_fp4_e2m1_bits v; v.u = bits; return v.v; }

__device__ __forceinline__ fpx_fp6_e3m2 fpx_f32_to_fp6_e3m2(float x) { return static_cast<fpx_fp6_e3m2>(x); }
__device__ __forceinline__ float fpx_fp6_e3m2_to_f32(fpx_fp6_e3m2 x) { return static_cast<float>(x); }
__device__ __forceinline__ fpx_fp6_e3m2 fpx_decode_fp6_e3m2(unsigned char bits) { fpx_fp6_e3m2_bits v; v.u = bits; return v.v; }

__device__ __forceinline__ fpx_fp8_e4m3 fpx_f32_to_fp8_e4m3(float x) { return static_cast<fpx_fp8_e4m3>(x); }
__device__ __forceinline__ float fpx_fp8_e4m3_to_f32(fpx_fp8_e4m3 x) { return static_cast<float>(x); }
__device__ __forceinline__ fpx_fp8_e4m3 fpx_decode_fp8_e4m3(unsigned char bits) { fpx_fp8_e4m3_bits v; v.u = bits; return v.v; }

__device__ __forceinline__ fpx_fp8_e5m2 fpx_f32_to_fp8_e5m2(float x) { return static_cast<fpx_fp8_e5m2>(x); }
__device__ __forceinline__ float fpx_fp8_e5m2_to_f32(fpx_fp8_e5m2 x) { return static_cast<float>(x); }
__device__ __forceinline__ fpx_fp8_e5m2 fpx_decode_fp8_e5m2(unsigned char bits) { fpx_fp8_e5m2_bits v; v.u = bits; return v.v; }

__device__ __forceinline__ fpx_half fpx_add_f16(fpx_half a, fpx_half b) { return __hadd(a, b); }
__device__ __forceinline__ fpx_half fpx_sub_f16(fpx_half a, fpx_half b) { return __hsub(a, b); }
__device__ __forceinline__ fpx_half fpx_mul_f16(fpx_half a, fpx_half b) { return __hmul(a, b); }
__device__ __forceinline__ fpx_half fpx_div_f16(fpx_half a, fpx_half b) { return __hdiv(a, b); }
__device__ __forceinline__ fpx_half fpx_abs_f16(fpx_half x) { return __habs(x); }
__device__ __forceinline__ bool fpx_le_f16(fpx_half a, fpx_half b) { return __hle(a, b); }
__device__ __forceinline__ fpx_half fpx_fma_f16(fpx_half a, fpx_half b, fpx_half c) { return __hfma(a, b, c); }
__device__ __forceinline__ fpx_half fpx_mul_add_f16(fpx_half a, fpx_half b, fpx_half c) {
#if FPX_CUDA_USE_FMA
    return __hfma(a, b, c);
#else
    return __hadd(__hmul(a, b), c);
#endif
}

// CUDA FP4/FP6/FP8 scalar wrappers provide storage and conversion, not scalar
// arithmetic. Each helper widens one operation to f32 and requantizes it.
#define FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC(suffix, scalar_ty, to_f32, from_f32)             \
    __device__ __forceinline__ scalar_ty fpx_add_##suffix(scalar_ty a, scalar_ty b) {           \
        return from_f32(to_f32(a) + to_f32(b));                                                  \
    }                                                                                            \
    __device__ __forceinline__ scalar_ty fpx_sub_##suffix(scalar_ty a, scalar_ty b) {           \
        return from_f32(to_f32(a) - to_f32(b));                                                  \
    }                                                                                            \
    __device__ __forceinline__ scalar_ty fpx_mul_##suffix(scalar_ty a, scalar_ty b) {           \
        return from_f32(to_f32(a) * to_f32(b));                                                  \
    }                                                                                            \
    __device__ __forceinline__ scalar_ty fpx_div_##suffix(scalar_ty a, scalar_ty b) {           \
        return from_f32(to_f32(a) / to_f32(b));                                                  \
    }                                                                                            \
    __device__ __forceinline__ scalar_ty fpx_fma_##suffix(                                      \
        scalar_ty a, scalar_ty b, scalar_ty c) {                                                 \
        return from_f32(fmaf(to_f32(a), to_f32(b), to_f32(c)));                                 \
    }                                                                                            \
    __device__ __forceinline__ scalar_ty fpx_mul_add_##suffix(                                  \
        scalar_ty a, scalar_ty b, scalar_ty c) {                                                 \
        if (FPX_CUDA_USE_FMA) return fpx_fma_##suffix(a, b, c);                                 \
        return from_f32(to_f32(a) * to_f32(b) + to_f32(c));                                     \
    }

FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC(fp4_e2m1, fpx_fp4_e2m1, fpx_fp4_e2m1_to_f32, fpx_f32_to_fp4_e2m1)
FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC(fp6_e3m2, fpx_fp6_e3m2, fpx_fp6_e3m2_to_f32, fpx_f32_to_fp6_e3m2)
FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC(fp8_e4m3, fpx_fp8_e4m3, fpx_fp8_e4m3_to_f32, fpx_f32_to_fp8_e4m3)
FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC(fp8_e5m2, fpx_fp8_e5m2, fpx_fp8_e5m2_to_f32, fpx_f32_to_fp8_e5m2)

#undef FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC
"#.to_string();
    }

    if backend == "rocm" {
        return r#"using fpx_half = __half;
using fpx_fp4_e2m1 = __hip_fp4_e2m1;
using fpx_fp6_e3m2 = __hip_fp6_e3m2;
using fpx_fp8_e4m3 = __hip_fp8_e4m3;
using fpx_fp8_e5m2 = __hip_fp8_e5m2;

#ifndef FPX_ROCM_USE_FMA
#define FPX_ROCM_USE_FMA 0
#endif

static_assert(sizeof(fpx_half) == sizeof(unsigned short), "expected 16-bit __half");
static_assert(sizeof(fpx_fp4_e2m1) == sizeof(unsigned char), "expected 8-bit __hip_fp4_e2m1 storage");
static_assert(sizeof(fpx_fp6_e3m2) == sizeof(unsigned char), "expected 8-bit __hip_fp6_e3m2 storage");
static_assert(sizeof(fpx_fp8_e4m3) == sizeof(unsigned char), "expected 8-bit __hip_fp8_e4m3 storage");
static_assert(sizeof(fpx_fp8_e5m2) == sizeof(unsigned char), "expected 8-bit __hip_fp8_e5m2 storage");

typedef union {
    unsigned short u;
    fpx_half h;
} fpx_half_bits;

typedef union {
    unsigned char u;
    fpx_fp4_e2m1 v;
} fpx_fp4_e2m1_bits;

typedef union {
    unsigned char u;
    fpx_fp6_e3m2 v;
} fpx_fp6_e3m2_bits;

typedef union {
    unsigned char u;
    fpx_fp8_e4m3 v;
} fpx_fp8_e4m3_bits;

typedef union {
    unsigned char u;
    fpx_fp8_e5m2 v;
} fpx_fp8_e5m2_bits;

__device__ __forceinline__ float fpx_mul_add_f32(float a, float b, float c) {
#if FPX_ROCM_USE_FMA
    return fmaf(a, b, c);
#else
    return a * b + c;
#endif
}

__device__ __forceinline__ double fpx_mul_add_f64(double a, double b, double c) {
#if FPX_ROCM_USE_FMA
    return fma(a, b, c);
#else
    return a * b + c;
#endif
}

__device__ __forceinline__ fpx_half fpx_f32_to_half(float x) { return __float2half(x); }
__device__ __forceinline__ float fpx_half_to_f32(fpx_half x) { return __half2float(x); }
__device__ __forceinline__ fpx_half fpx_decode_f16_bits(unsigned short bits) { fpx_half_bits v; v.u = bits; return v.h; }

__device__ __forceinline__ fpx_fp4_e2m1 fpx_f32_to_fp4_e2m1(float x) { return fpx_fp4_e2m1(x); }
__device__ __forceinline__ float fpx_fp4_e2m1_to_f32(fpx_fp4_e2m1 x) { return static_cast<float>(x); }
__device__ __forceinline__ fpx_fp4_e2m1 fpx_decode_fp4_e2m1(unsigned char bits) { fpx_fp4_e2m1_bits v; v.u = bits; return v.v; }

__device__ __forceinline__ fpx_fp6_e3m2 fpx_f32_to_fp6_e3m2(float x) { return fpx_fp6_e3m2(x); }
__device__ __forceinline__ float fpx_fp6_e3m2_to_f32(fpx_fp6_e3m2 x) { return static_cast<float>(x); }
__device__ __forceinline__ fpx_fp6_e3m2 fpx_decode_fp6_e3m2(unsigned char bits) { fpx_fp6_e3m2_bits v; v.u = bits; return v.v; }

__device__ __forceinline__ fpx_fp8_e4m3 fpx_f32_to_fp8_e4m3(float x) { return fpx_fp8_e4m3(x); }
__device__ __forceinline__ float fpx_fp8_e4m3_to_f32(fpx_fp8_e4m3 x) { return static_cast<float>(x); }
__device__ __forceinline__ fpx_fp8_e4m3 fpx_decode_fp8_e4m3(unsigned char bits) { fpx_fp8_e4m3_bits v; v.u = bits; return v.v; }

__device__ __forceinline__ fpx_fp8_e5m2 fpx_f32_to_fp8_e5m2(float x) { return fpx_fp8_e5m2(x); }
__device__ __forceinline__ float fpx_fp8_e5m2_to_f32(fpx_fp8_e5m2 x) { return static_cast<float>(x); }
__device__ __forceinline__ fpx_fp8_e5m2 fpx_decode_fp8_e5m2(unsigned char bits) { fpx_fp8_e5m2_bits v; v.u = bits; return v.v; }

__device__ __forceinline__ fpx_half fpx_add_f16(fpx_half a, fpx_half b) { return __hadd(a, b); }
__device__ __forceinline__ fpx_half fpx_sub_f16(fpx_half a, fpx_half b) { return __hsub(a, b); }
__device__ __forceinline__ fpx_half fpx_mul_f16(fpx_half a, fpx_half b) { return __hmul(a, b); }
__device__ __forceinline__ fpx_half fpx_div_f16(fpx_half a, fpx_half b) { return __hdiv(a, b); }
__device__ __forceinline__ fpx_half fpx_abs_f16(fpx_half x) { return __habs(x); }
__device__ __forceinline__ bool fpx_le_f16(fpx_half a, fpx_half b) { return __hle(a, b); }
__device__ __forceinline__ fpx_half fpx_fma_f16(fpx_half a, fpx_half b, fpx_half c) { return __hfma(a, b, c); }
__device__ __forceinline__ fpx_half fpx_mul_add_f16(fpx_half a, fpx_half b, fpx_half c) {
#if FPX_ROCM_USE_FMA
    return __hfma(a, b, c);
#else
    return __hadd(__hmul(a, b), c);
#endif
}

// HIP FP4/FP6/FP8 scalar wrappers provide storage and conversion, not scalar
// arithmetic. Each helper widens one operation to f32 and requantizes it.
#define FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC(suffix, scalar_ty, to_f32, from_f32)             \
    __device__ __forceinline__ scalar_ty fpx_add_##suffix(scalar_ty a, scalar_ty b) {           \
        return from_f32(to_f32(a) + to_f32(b));                                                  \
    }                                                                                            \
    __device__ __forceinline__ scalar_ty fpx_sub_##suffix(scalar_ty a, scalar_ty b) {           \
        return from_f32(to_f32(a) - to_f32(b));                                                  \
    }                                                                                            \
    __device__ __forceinline__ scalar_ty fpx_mul_##suffix(scalar_ty a, scalar_ty b) {           \
        return from_f32(to_f32(a) * to_f32(b));                                                  \
    }                                                                                            \
    __device__ __forceinline__ scalar_ty fpx_div_##suffix(scalar_ty a, scalar_ty b) {           \
        return from_f32(to_f32(a) / to_f32(b));                                                  \
    }                                                                                            \
    __device__ __forceinline__ scalar_ty fpx_fma_##suffix(                                      \
        scalar_ty a, scalar_ty b, scalar_ty c) {                                                 \
        return from_f32(fmaf(to_f32(a), to_f32(b), to_f32(c)));                                 \
    }                                                                                            \
    __device__ __forceinline__ scalar_ty fpx_mul_add_##suffix(                                  \
        scalar_ty a, scalar_ty b, scalar_ty c) {                                                 \
        if (FPX_ROCM_USE_FMA) return fpx_fma_##suffix(a, b, c);                                 \
        return from_f32(to_f32(a) * to_f32(b) + to_f32(c));                                     \
    }

FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC(fp4_e2m1, fpx_fp4_e2m1, fpx_fp4_e2m1_to_f32, fpx_f32_to_fp4_e2m1)
FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC(fp6_e3m2, fpx_fp6_e3m2, fpx_fp6_e3m2_to_f32, fpx_f32_to_fp6_e3m2)
FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC(fp8_e4m3, fpx_fp8_e4m3, fpx_fp8_e4m3_to_f32, fpx_f32_to_fp8_e4m3)
FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC(fp8_e5m2, fpx_fp8_e5m2, fpx_fp8_e5m2_to_f32, fpx_f32_to_fp8_e5m2)

#undef FPX_DEFINE_EMULATED_SCALAR_ARITHMETIC
"#.to_string();
    }

    let half_ty = match backend {
        "cuda" => "__half",
        "rocm" => "half",
        _ => "float",
    };

    format!(
        r#"using fpx_half = {half_ty};

typedef union {{
    unsigned short u;
    fpx_half h;
}} fpx_half_bits;

__device__ __forceinline__ float fpx_abs_f32(float x) {{ return x < 0.0f ? -x : x; }}
__device__ __forceinline__ int fpx_is_nan_f32(float x) {{ return x != x; }}
__device__ __forceinline__ int fpx_is_inf_f32(float x) {{ return !fpx_is_nan_f32(x) && fpx_is_nan_f32(x - x); }}
__device__ __forceinline__ float fpx_exp2i(int e) {{
    float out = 1.0f;
    if (e >= 0) {{ for (int i = 0; i < e; ++i) out *= 2.0f; }}
    else {{ for (int i = 0; i < -e; ++i) out *= 0.5f; }}
    return out;
}}
__device__ __forceinline__ fpx_half fpx_f32_to_half(float x) {{ return __float2half_rn(x); }}
__device__ __forceinline__ float fpx_half_to_f32(fpx_half x) {{ return __half2float(x); }}
__device__ __forceinline__ float fpx_decode_f16_bits(unsigned short bits) {{ fpx_half_bits v; v.u = bits; return fpx_half_to_f32(v.h); }}
__device__ __forceinline__ unsigned short fpx_encode_f16_bits(float x) {{ fpx_half_bits v; v.h = fpx_f32_to_half(x); return v.u; }}
__device__ __forceinline__ float fpx_quantize_f16(float x) {{ return fpx_decode_f16_bits(fpx_encode_f16_bits(x)); }}

__device__ __forceinline__ float fpx_decode_fp4_e2m1(unsigned char bits) {{
    const float mags[8] = {{0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f, 4.0f, 6.0f}};
    bits &= 0x0F;
    unsigned char sign = (bits >> 3) & 0x1;
    unsigned char idx = bits & 0x07;
    float mag = mags[idx];
    if (mag == 0.0f) return sign ? -0.0f : 0.0f;
    return sign ? -mag : mag;
}}

__device__ __forceinline__ unsigned char fpx_quantize_fp4_e2m1_bits(float v) {{
    if (v == 0.0f || fpx_is_nan_f32(v) || fpx_is_inf_f32(v)) return 0;
    unsigned char sign = v < 0.0f ? 1u : 0u;
    float mag = fpx_abs_f32(v);
    const float mags[8] = {{0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f, 4.0f, 6.0f}};
    unsigned char best = 1;
    float best_err = 3.402823466e38f;
    for (unsigned char idx = 1; idx < 8; ++idx) {{
        float err = fpx_abs_f32(mag - mags[idx]);
        if (err < best_err) {{ best_err = err; best = idx; }}
    }}
    return (unsigned char)((sign << 3) | best);
}}
__device__ __forceinline__ float fpx_quantize_fp4_e2m1(float x) {{ return fpx_decode_fp4_e2m1(fpx_quantize_fp4_e2m1_bits(x)); }}

__device__ __forceinline__ float fpx_decode_fp6_e3m2(unsigned char bits) {{
    bits &= 0x3F;
    unsigned char s = (bits >> 5) & 0x1;
    unsigned char e = (bits >> 2) & 0x07;
    unsigned char m = bits & 0x03;
    float sign = s == 0 ? 1.0f : -1.0f;
    if (e == 0) {{
        if (m == 0) return sign * 0.0f;
        return sign * fpx_exp2i(-2) * ((float)m / 4.0f);
    }}
    return sign * fpx_exp2i((int)e - 3) * (1.0f + (float)m / 4.0f);
}}

__device__ __forceinline__ unsigned char fpx_quantize_fp6_e3m2_bits(float v) {{
    if (fpx_is_nan_f32(v)) return 0;
    unsigned char best = 0;
    float best_err = 3.402823466e38f;
    for (unsigned int bits = 0; bits < 64; ++bits) {{
        float decoded = fpx_decode_fp6_e3m2((unsigned char)bits);
        float err = fpx_abs_f32(decoded - v);
        if (err < best_err || (err == best_err && (((unsigned char)bits & 0x1u) == 0u) && ((best & 0x1u) != 0u))) {{
            best_err = err;
            best = (unsigned char)bits;
        }}
    }}
    return best;
}}
__device__ __forceinline__ float fpx_quantize_fp6_e3m2(float x) {{ return fpx_decode_fp6_e3m2(fpx_quantize_fp6_e3m2_bits(x)); }}

__device__ __forceinline__ float fpx_decode_fp8_e5m2(unsigned char bits) {{
    unsigned char s = (bits >> 7) & 0x1;
    unsigned char e = (bits >> 2) & 0x1F;
    unsigned char m = bits & 0x03;
    float sign = s == 0 ? 1.0f : -1.0f;
    if (e == 0) {{
        if (m == 0) return sign * 0.0f;
        return sign * fpx_exp2i(1 - 15) * ((float)m / 4.0f);
    }}
    if (e == 0x1F) {{
        return m == 0 ? sign * (1.0f / 0.0f) : (0.0f / 0.0f);
    }}
    return sign * fpx_exp2i((int)e - 15) * (1.0f + (float)m / 4.0f);
}}

__device__ __forceinline__ unsigned char fpx_quantize_fp8_e5m2_bits(float v) {{
    if (fpx_is_nan_f32(v)) return 0x7F;
    if (fpx_is_inf_f32(v)) return v > 0.0f ? 0x7C : 0xFC;
    unsigned char best = 0;
    float best_err = 3.402823466e38f;
    for (unsigned int bits = 0; bits <= 255; ++bits) {{
        float decoded = fpx_decode_fp8_e5m2((unsigned char)bits);
        if (fpx_is_nan_f32(decoded) || fpx_is_inf_f32(decoded)) continue;
        float err = fpx_abs_f32(decoded - v);
        if (err < best_err || (err == best_err && (((unsigned char)bits & 0x1u) == 0u) && ((best & 0x1u) != 0u))) {{
            best_err = err;
            best = (unsigned char)bits;
        }}
    }}
    return best;
}}
__device__ __forceinline__ float fpx_quantize_fp8_e5m2(float x) {{ return fpx_decode_fp8_e5m2(fpx_quantize_fp8_e5m2_bits(x)); }}

__device__ __forceinline__ float fpx_decode_fp8_e4m3(unsigned char bits) {{
    unsigned char s = (bits >> 7) & 0x1;
    unsigned char e = (bits >> 3) & 0x0F;
    unsigned char m = bits & 0x07;
    float sign = s == 0 ? 1.0f : -1.0f;
    if (e == 0) {{
        if (m == 0) return sign * 0.0f;
        return sign * fpx_exp2i(1 - 7) * ((float)m / 8.0f);
    }}
    if (e == 0x0F && m == 0x07) return 0.0f / 0.0f;
    return sign * fpx_exp2i((int)e - 7) * (1.0f + (float)m / 8.0f);
}}

__device__ __forceinline__ unsigned char fpx_quantize_fp8_e4m3_bits(float v) {{
    if (fpx_is_nan_f32(v)) return 0x7F;
    if (fpx_is_inf_f32(v)) v = v > 0.0f ? 3.402823466e38f : -3.402823466e38f;
    unsigned char best = 0;
    float best_err = 3.402823466e38f;
    for (unsigned int bits = 0; bits <= 255; ++bits) {{
        float decoded = fpx_decode_fp8_e4m3((unsigned char)bits);
        if (fpx_is_nan_f32(decoded)) continue;
        float err = fpx_abs_f32(decoded - v);
        if (err < best_err || (err == best_err && (((unsigned char)bits & 0x1u) == 0u) && ((best & 0x1u) != 0u))) {{
            best_err = err;
            best = (unsigned char)bits;
        }}
    }}
    return best;
}}
__device__ __forceinline__ float fpx_quantize_fp8_e4m3(float x) {{ return fpx_decode_fp8_e4m3(fpx_quantize_fp8_e4m3_bits(x)); }}

"#
    )
}

fn gpu_eval_source(
    backend: &str,
    func_name: &str,
    segment_tag: &str,
    dom: &str,
    scalar_ty: &str,
    coeff_symbol: &str,
    coeff_len: usize,
) -> String {
    let ident = sanitize_ident(func_name);
    let fn_name = format!("{}_{}_{}", ident, segment_tag, scalar_ty);
    let device_inline = match backend {
        "cuda" | "rocm" => "__device__ __forceinline__",
        _ => "inline",
    };
    let mul_add = match (backend, scalar_ty) {
        ("cuda", "float") | ("rocm", "float") => "fpx_mul_add_f32",
        ("cuda", "double") | ("rocm", "double") => "fpx_mul_add_f64",
        _ => "",
    };

    let body = if func_name == "log" {
        // Log acceptance is evaluated with fused Horner for f32/f64.  Emit an
        // explicit backend FMA so the checked kernel does not depend on the
        // optional global FPX_*_USE_FMA switch.
        let log_mul_add = match (backend, scalar_ty) {
            ("cuda" | "rocm", "float") => "fmaf",
            ("cuda" | "rocm", "double") => "fma",
            _ => "",
        };
        let step = if log_mul_add.is_empty() {
            format!("p = p * z + {coeff_symbol}[i];")
        } else {
            format!("p = {log_mul_add}(p, z, {coeff_symbol}[i]);")
        };
        let reconstructed = if log_mul_add.is_empty() {
            "one + z * p".to_string()
        } else {
            format!("{log_mul_add}(p, z, one)")
        };
        format!(
            "    const {scalar_ty} one = ({scalar_ty})1.0;\n    {scalar_ty} s = (x - one) / (x + one);\n    {scalar_ty} z = s * s;\n    {scalar_ty} p = {coeff_symbol}[{last}];\n    for (int i = {last_minus_one}; i >= 0; --i) {{ {step} }}\n    return ({scalar_ty})2.0 * s * ({reconstructed});",
            last = coeff_len.saturating_sub(1),
            last_minus_one = coeff_len.saturating_sub(2),
        )
    } else if dom == "z" {
        match func_name {
            "sin" => {
                let step = if mul_add.is_empty() {
                    format!("p = p * z + {coeff_symbol}[i];")
                } else {
                    format!("p = {mul_add}(p, z, {coeff_symbol}[i]);")
                };
                format!(
                    "    {scalar_ty} z = x * x;\n    {scalar_ty} p = {coeff_symbol}[{last}];\n    for (int i = {last_minus_one}; i >= 0; --i) {{ {step} }}\n    return x * p;",
                    last = coeff_len.saturating_sub(1),
                    last_minus_one = coeff_len.saturating_sub(2),
                )
            }
            "tan" => {
                let step = if mul_add.is_empty() {
                    format!("p = p * z + {coeff_symbol}[i];")
                } else {
                    format!("p = {mul_add}(p, z, {coeff_symbol}[i]);")
                };
                format!(
                    "    {scalar_ty} z = x * x;\n    {scalar_ty} p = {coeff_symbol}[{last}];\n    for (int i = {last_minus_one}; i >= 0; --i) {{ {step} }}\n    return x + x * z * p;",
                    last = coeff_len.saturating_sub(1),
                    last_minus_one = coeff_len.saturating_sub(2),
                )
            }
            _ => {
                let step = if mul_add.is_empty() {
                    format!("p = p * z + {coeff_symbol}[i];")
                } else {
                    format!("p = {mul_add}(p, z, {coeff_symbol}[i]);")
                };
                format!(
                    "    {scalar_ty} z = x * x;\n    {scalar_ty} p = {coeff_symbol}[{last}];\n    for (int i = {last_minus_one}; i >= 0; --i) {{ {step} }}\n    return p;",
                    last = coeff_len.saturating_sub(1),
                    last_minus_one = coeff_len.saturating_sub(2),
                )
            }
        }
    } else {
        let step = if mul_add.is_empty() {
            format!("p = p * x + {coeff_symbol}[i];")
        } else {
            format!("p = {mul_add}(p, x, {coeff_symbol}[i]);")
        };
        format!(
            "    {scalar_ty} p = {coeff_symbol}[{last}];\n    for (int i = {last_minus_one}; i >= 0; --i) {{ {step} }}\n    return p;",
            last = coeff_len.saturating_sub(1),
            last_minus_one = coeff_len.saturating_sub(2),
        )
    };

    format!("{device_inline} {scalar_ty} {fn_name}({scalar_ty} x) {{\n{body}\n}}")
}

fn gpu_file_ext(backend: &str) -> &'static str {
    match backend {
        "cuda" => "cu",
        "rocm" => "hip",
        _ => "txt",
    }
}

fn gpu_rel_path(backend: &str, func_name: &str) -> String {
    format!(
        "{backend}/{}.{}",
        sanitize_ident(func_name),
        gpu_file_ext(backend)
    )
}

fn gpu_include(backend: &str) -> &'static str {
    match backend {
        "cuda" => "#include <cuda_runtime.h>\n#include <cuda_fp16.h>\n#include <cuda_fp4.h>\n#include <cuda_fp6.h>\n#include <cuda_fp8.h>",
        "rocm" => "#include <hip/hip_runtime.h>\n#include <hip/hip_fp16.h>\n#include <hip/hip_fp4.h>\n#include <hip/hip_fp6.h>\n#include <hip/hip_fp8.h>",
        _ => "",
    }
}

fn gpu_backend_lowp_type(backend: &str, suffix: &str) -> Option<&'static str> {
    match (backend, suffix) {
        ("cuda", "f16") => Some("__half"),
        ("cuda", "fp4_e2m1") => Some("__nv_fp4_e2m1"),
        ("cuda", "fp6_e3m2") => Some("__nv_fp6_e3m2"),
        ("cuda", "fp8_e4m3") => Some("__nv_fp8_e4m3"),
        ("cuda", "fp8_e5m2") => Some("__nv_fp8_e5m2"),
        ("rocm", "f16") => Some("__half"),
        ("rocm", "fp4_e2m1") => Some("__hip_fp4_e2m1"),
        ("rocm", "fp6_e3m2") => Some("__hip_fp6_e3m2"),
        ("rocm", "fp8_e4m3") => Some("__hip_fp8_e4m3"),
        ("rocm", "fp8_e5m2") => Some("__hip_fp8_e5m2"),
        _ => None,
    }
}

fn gpu_lowp_scalar_alias(suffix: &str) -> Option<&'static str> {
    match suffix {
        "f16" => Some("fpx_half"),
        "fp4_e2m1" => Some("fpx_fp4_e2m1"),
        "fp6_e3m2" => Some("fpx_fp6_e3m2"),
        "fp8_e4m3" => Some("fpx_fp8_e4m3"),
        "fp8_e5m2" => Some("fpx_fp8_e5m2"),
        _ => None,
    }
}

fn gpu_lowp_to_f32_fn(suffix: &str) -> Option<&'static str> {
    match suffix {
        "f16" => Some("fpx_half_to_f32"),
        "fp4_e2m1" => Some("fpx_fp4_e2m1_to_f32"),
        "fp6_e3m2" => Some("fpx_fp6_e3m2_to_f32"),
        "fp8_e4m3" => Some("fpx_fp8_e4m3_to_f32"),
        "fp8_e5m2" => Some("fpx_fp8_e5m2_to_f32"),
        _ => None,
    }
}

fn gpu_lowp_from_f32_fn(suffix: &str) -> Option<&'static str> {
    match suffix {
        "f16" => Some("fpx_f32_to_half"),
        "fp4_e2m1" => Some("fpx_f32_to_fp4_e2m1"),
        "fp6_e3m2" => Some("fpx_f32_to_fp6_e3m2"),
        "fp8_e4m3" => Some("fpx_f32_to_fp8_e4m3"),
        "fp8_e5m2" => Some("fpx_f32_to_fp8_e5m2"),
        _ => None,
    }
}

fn gpu_eval_fn_name(func_name: &str, segment_idx: usize, scalar_ty: &str) -> String {
    format!(
        "{}_seg{}_{}",
        sanitize_ident(func_name),
        segment_idx,
        scalar_ty
    )
}

fn gpu_raw_eval_fn_name(func_name: &str, scalar_ty: &str) -> String {
    format!("{}_kernel_eval_{}", sanitize_ident(func_name), scalar_ty)
}

fn gpu_public_eval_fn_name(func_name: &str, scalar_ty: &str) -> String {
    format!("{}_eval_{}", sanitize_ident(func_name), scalar_ty)
}

/// Whether the GPU exporter has a complete public implementation for this
/// function, including full-domain reduction, reconstruction, and special
/// values.  A table-domain polynomial alone is not a public implementation.
fn gpu_has_full_domain_implementation(func_name: &str) -> bool {
    matches!(func_name, "log")
}

fn gpu_uses_abs_dispatch(func_name: &str) -> bool {
    matches!(func_name, "sin" | "cos" | "tan")
}

fn gpu_scalar_literal(value: f64, scalar_ty: &str) -> String {
    if scalar_ty == "double" {
        format!("{value:.17e}")
    } else {
        format!("{:.9e}f", value as f32)
    }
}

fn gpu_math_fn_names(scalar_ty: &str) -> (&'static str, &'static str, &'static str, &'static str) {
    if scalar_ty == "double" {
        ("fabs", "remainder", "copysign", "nan(\"\")")
    } else {
        ("fabsf", "remainderf", "copysignf", "nanf(\"\")")
    }
}

fn render_gpu_raw_dispatch(export: &FunctionRemezExport, scalar_ty: &str) -> String {
    let fn_name = gpu_raw_eval_fn_name(&export.function, scalar_ty);
    let dispatch_var = if gpu_uses_abs_dispatch(&export.function) {
        "ax"
    } else {
        "x"
    };
    let abs_line = if gpu_uses_abs_dispatch(&export.function) {
        let abs_fn = if scalar_ty == "double" {
            "fabs"
        } else {
            "fabsf"
        };
        format!("    {scalar_ty} ax = {abs_fn}(x);\n")
    } else {
        String::new()
    };
    let dispatch = if export.segments.len() == 1 {
        format!(
            "    return {}(x);",
            gpu_eval_fn_name(&export.function, 0, scalar_ty)
        )
    } else {
        let mut lines = Vec::new();
        for (idx, seg) in export.segments.iter().enumerate() {
            let kw = if idx == 0 { "if" } else { "else if" };
            let bound = if scalar_ty == "double" {
                format!("{:.17e}", seg.interval_x.1)
            } else {
                format!("{:.9e}f", seg.interval_x.1 as f32)
            };
            lines.push(format!(
                "    {kw} ({dispatch_var} <= {bound}) {{ return {}(x); }}",
                gpu_eval_fn_name(&export.function, idx, scalar_ty)
            ));
        }
        let last = export.segments.len() - 1;
        lines.push(format!(
            "    else {{ return {}(x); }}",
            gpu_eval_fn_name(&export.function, last, scalar_ty)
        ));
        lines.join("\n")
    };

    format!("__device__ __forceinline__ {scalar_ty} {fn_name}({scalar_ty} x) {{\n{abs_line}{dispatch}\n}}")
}

fn render_gpu_public_dispatch(export: &FunctionRemezExport, scalar_ty: &str) -> String {
    match export.function.as_str() {
        "log" => render_gpu_log_public_dispatch(scalar_ty),
        function => panic!(
            "GPU full-domain FPX dispatch is not implemented for {function}; \
             refusing to substitute a vendor intrinsic"
        ),
    }
}

fn render_gpu_log_public_dispatch(scalar_ty: &str) -> String {
    match scalar_ty {
        "double" => r#"__device__ __forceinline__ double log_eval_double(double x) {
    if (isnan(x)) { return x; }
    if (x == 0.0) { return -INFINITY; }
    if (x < 0.0) { return nan(""); }
    if (isinf(x)) { return INFINITY; }
    if (x == 1.0) { return 0.0; }
    int e = 0;
    double m = frexp(x, &e);
    if (m < 7.07106781186547573e-1) { m *= 2.0; --e; }
    double base = log_kernel_eval_double(m);
    double ef = (double)e;
    return (base + ef * 1.90821492927058770002e-10) + ef * 6.93147180369123816490e-1;
}"#
            .to_string(),
        "float" => r#"__device__ __forceinline__ float log_eval_float(float x) {
    if (isnan(x)) { return x; }
    if (x == 0.0f) { return -INFINITY; }
    if (x < 0.0f) { return nanf(""); }
    if (isinf(x)) { return INFINITY; }
    if (x == 1.0f) { return 0.0f; }
    int e = 0;
    float m = frexpf(x, &e);
    if (m < 7.071067691e-1f) { m *= 2.0f; --e; }
    float base = log_kernel_eval_float(m);
    double ef = (double)e;
    return (float)(((double)base + ef * 1.90821492927058770002e-10) + ef * 6.93147180369123816490e-1);
}"#
            .to_string(),
        _ => panic!("unsupported GPU log scalar type {scalar_ty}"),
    }
}

fn render_gpu_kernel_entry(export: &FunctionRemezExport, scalar_ty: &str) -> String {
    let ident = sanitize_ident(&export.function);
    let kernel_name = format!("{}_kernel_{}", ident, scalar_ty);
    let eval_name = gpu_public_eval_fn_name(&export.function, scalar_ty);
    format!(
        "extern \"C\" __global__ void {kernel_name}(const {scalar_ty}* xs, {scalar_ty}* ys, size_t n) {{\n    size_t i = (size_t)(blockIdx.x * blockDim.x + threadIdx.x);\n    if (i < n) {{ ys[i] = {eval_name}(xs[i]); }}\n}}"
    )
}

fn render_gpu_reduced_kernel_entry(export: &FunctionRemezExport, scalar_ty: &str) -> String {
    let ident = sanitize_ident(&export.function);
    let kernel_name = format!("{}_reduced_kernel_{}", ident, scalar_ty);
    let eval_name = gpu_raw_eval_fn_name(&export.function, scalar_ty);
    format!(
        "extern \"C\" __global__ void {kernel_name}(const {scalar_ty}* xs, {scalar_ty}* ys, size_t n) {{\n    size_t i = (size_t)(blockIdx.x * blockDim.x + threadIdx.x);\n    if (i < n) {{ ys[i] = {eval_name}(xs[i]); }}\n}}"
    )
}

fn render_gpu_reduced_domain_accessors(export: &FunctionRemezExport) -> String {
    let first = export.segments.first().unwrap_or_else(|| {
        panic!(
            "cannot export GPU reduced-domain kernels for {} without any segments",
            export.function
        )
    });
    let last = export.segments.last().expect("first segment exists");
    let ident = sanitize_ident(&export.function);
    let lo = gpu_scalar_literal(first.interval_x.0, "double");
    let hi = gpu_scalar_literal(last.interval_x.1, "double");
    format!(
        "// Reduced-domain kernels accept inputs only in [{lo}, {hi}].\n\
         extern \"C\" double {ident}_reduced_domain_lo() {{ return {lo}; }}\n\
         extern \"C\" double {ident}_reduced_domain_hi() {{ return {hi}; }}"
    )
}

fn gpu_lowp_coeff_symbol(func_name: &str, segment_idx: usize, suffix: &str) -> String {
    format!(
        "{}_seg{}_coeffs_{suffix}",
        sanitize_ident(func_name),
        segment_idx
    )
}

fn render_gpu_lowp_segment_eval(
    func_name: &str,
    segment_idx: usize,
    table_domain: &str,
    suffix: &str,
    coeff_symbol: &str,
    coeff_len: usize,
    coeff_decode: &str,
) -> String {
    let scalar_ty = gpu_lowp_scalar_alias(suffix)
        .unwrap_or_else(|| panic!("unsupported GPU low-precision type {suffix}"));
    let from_f32 = gpu_lowp_from_f32_fn(suffix)
        .unwrap_or_else(|| panic!("missing GPU constructor for {suffix}"));
    let fn_name = format!(
        "{}_seg{}_{}",
        sanitize_ident(func_name),
        segment_idx,
        suffix
    );
    let last = coeff_len.saturating_sub(1);
    let last_minus_one = coeff_len.saturating_sub(2);
    let add = format!("fpx_add_{suffix}");
    let sub = format!("fpx_sub_{suffix}");
    let mul = format!("fpx_mul_{suffix}");
    let div = format!("fpx_div_{suffix}");
    let fma = format!("fpx_fma_{suffix}");
    let mul_add = format!("fpx_mul_add_{suffix}");
    let one = format!("{from_f32}(1.0f)");
    let two = format!("{from_f32}(2.0f)");

    let body = if func_name == "log" {
        let step = format!("p = {fma}(p, z, {coeff_decode}({coeff_symbol}[i]));");
        format!(
            "    const {scalar_ty} one = {one};\n    const {scalar_ty} two = {two};\n    {scalar_ty} u = {sub}(x, one);\n    {scalar_ty} s = {div}(u, {add}(x, one));\n    {scalar_ty} z = {mul}(s, s);\n    {scalar_ty} p = {coeff_decode}({coeff_symbol}[{last}]);\n    for (int i = {last_minus_one}; i >= 0; --i) {{ {step} }}\n    {scalar_ty} reconstructed = {fma}(p, z, one);\n    return {mul}({mul}(two, s), reconstructed);"
        )
    } else if table_domain == "z" {
        if func_name == "sin" {
            let step = format!("p = {mul_add}(p, z, {coeff_decode}({coeff_symbol}[i]));");
            format!(
                "    {scalar_ty} z = {mul}(x, x);\n    {scalar_ty} p = {coeff_decode}({coeff_symbol}[{last}]);\n    for (int i = {last_minus_one}; i >= 0; --i) {{ {step} }}\n    return {mul}(x, p);"
            )
        } else if func_name == "tan" {
            let step = format!("p = {mul_add}(p, z, {coeff_decode}({coeff_symbol}[i]));");
            format!(
                "    {scalar_ty} z = {mul}(x, x);\n    {scalar_ty} p = {coeff_decode}({coeff_symbol}[{last}]);\n    for (int i = {last_minus_one}; i >= 0; --i) {{ {step} }}\n    return {add}(x, {mul}({mul}(x, z), p));"
            )
        } else {
            let step = format!("p = {mul_add}(p, z, {coeff_decode}({coeff_symbol}[i]));");
            format!(
                "    {scalar_ty} z = {mul}(x, x);\n    {scalar_ty} p = {coeff_decode}({coeff_symbol}[{last}]);\n    for (int i = {last_minus_one}; i >= 0; --i) {{ {step} }}\n    return p;"
            )
        }
    } else {
        let step = format!("p = {mul_add}(p, x, {coeff_decode}({coeff_symbol}[i]));");
        format!(
            "    {scalar_ty} p = {coeff_decode}({coeff_symbol}[{last}]);\n    for (int i = {last_minus_one}; i >= 0; --i) {{ {step} }}\n    return p;"
        )
    };

    format!("__device__ __forceinline__ {scalar_ty} {fn_name}({scalar_ty} x) {{\n{body}\n}}")
}

fn render_gpu_lowp_raw_dispatch(export: &FunctionRemezExport, suffix: &str) -> String {
    let scalar_ty = gpu_lowp_scalar_alias(suffix)
        .unwrap_or_else(|| panic!("unsupported GPU low-precision type {suffix}"));
    let to_f32 =
        gpu_lowp_to_f32_fn(suffix).unwrap_or_else(|| panic!("missing GPU conversion for {suffix}"));
    let fn_name = gpu_raw_eval_fn_name(&export.function, suffix);
    let native_half_dispatch = suffix == "f16";
    let dispatch_line = if export.segments.len() > 1 && native_half_dispatch {
        if gpu_uses_abs_dispatch(&export.function) {
            "    fpx_half dispatch_x = fpx_abs_f16(x);\n".to_string()
        } else {
            "    fpx_half dispatch_x = x;\n".to_string()
        }
    } else if export.segments.len() > 1 && gpu_uses_abs_dispatch(&export.function) {
        format!("    float dispatch_x = fabsf({to_f32}(x));\n")
    } else if export.segments.len() > 1 {
        format!("    float dispatch_x = {to_f32}(x);\n")
    } else {
        String::new()
    };
    let dispatch = if export.segments.len() == 1 {
        format!(
            "    return {}_seg0_{}(x);",
            sanitize_ident(&export.function),
            suffix
        )
    } else {
        let mut lines = Vec::new();
        for (idx, seg) in export.segments.iter().enumerate() {
            let kw = if idx == 0 { "if" } else { "else if" };
            let bound = format!("{:.9e}f", seg.interval_x.1 as f32);
            let condition = if native_half_dispatch {
                format!("fpx_le_f16(dispatch_x, fpx_f32_to_half({bound}))")
            } else {
                format!("dispatch_x <= {bound}")
            };
            lines.push(format!(
                "    {kw} ({condition}) {{ return {}_seg{}_{}(x); }}",
                sanitize_ident(&export.function),
                idx,
                suffix
            ));
        }
        let last = export.segments.len() - 1;
        lines.push(format!(
            "    else {{ return {}_seg{}_{}(x); }}",
            sanitize_ident(&export.function),
            last,
            suffix
        ));
        lines.join("\n")
    };
    format!(
        "__device__ __forceinline__ {scalar_ty} {fn_name}({scalar_ty} x) {{\n{dispatch_line}{dispatch}\n}}"
    )
}

fn render_gpu_lowp_dispatch(export: &FunctionRemezExport, suffix: &str) -> String {
    match export.function.as_str() {
        "log" => render_gpu_log_lowp_dispatch(suffix),
        function => panic!(
            "GPU full-domain FPX dispatch is not implemented for {function}; \
             refusing to substitute a vendor intrinsic"
        ),
    }
}

fn render_gpu_log_lowp_dispatch(suffix: &str) -> String {
    let scalar_ty = gpu_lowp_scalar_alias(suffix)
        .unwrap_or_else(|| panic!("unsupported GPU low-precision type {suffix}"));
    let to_f32 =
        gpu_lowp_to_f32_fn(suffix).unwrap_or_else(|| panic!("missing GPU conversion for {suffix}"));
    let from_f32 = gpu_lowp_from_f32_fn(suffix)
        .unwrap_or_else(|| panic!("missing GPU constructor for {suffix}"));
    let add = format!("fpx_add_{suffix}");
    let mul = format!("fpx_mul_{suffix}");
    let fn_name = gpu_public_eval_fn_name("log", suffix);
    let raw_name = gpu_raw_eval_fn_name("log", suffix);
    format!(
        r#"__device__ __forceinline__ {scalar_ty} {fn_name}({scalar_ty} x) {{
    float xf = {to_f32}(x);
    if (isnan(xf)) {{ return x; }}
    if (xf == 0.0f) {{ return {from_f32}(-INFINITY); }}
    if (xf < 0.0f) {{ return {from_f32}(nanf("")); }}
    if (isinf(xf)) {{ return {from_f32}(INFINITY); }}
    if (xf == 1.0f) {{ return {from_f32}(0.0f); }}
    int e = 0;
    float mf = frexpf(xf, &e);
    if (mf < 7.071067691e-1f) {{ mf *= 2.0f; --e; }}
    {scalar_ty} m = {from_f32}(mf);
    {scalar_ty} base = {raw_name}(m);
    return {add}(base, {mul}({from_f32}((float)e), {from_f32}(6.931471806e-1f)));
}}"#
    )
}

fn render_gpu_lowp_kernel_entry(
    export: &FunctionRemezExport,
    backend: &str,
    suffix: &str,
) -> String {
    render_gpu_lowp_kernel_entry_impl(export, backend, suffix, false)
}

fn render_gpu_lowp_reduced_kernel_entry(
    export: &FunctionRemezExport,
    backend: &str,
    suffix: &str,
) -> String {
    render_gpu_lowp_kernel_entry_impl(export, backend, suffix, true)
}

fn render_gpu_lowp_kernel_entry_impl(
    export: &FunctionRemezExport,
    backend: &str,
    suffix: &str,
    reduced: bool,
) -> String {
    let ident = sanitize_ident(&export.function);
    let kernel_name = if reduced {
        format!("{}_reduced_kernel_{}", ident, suffix)
    } else {
        format!("{}_kernel_{}", ident, suffix)
    };
    let eval_name = if reduced {
        gpu_raw_eval_fn_name(&export.function, suffix)
    } else {
        gpu_public_eval_fn_name(&export.function, suffix)
    };
    let scalar_ty = gpu_backend_lowp_type(backend, suffix)
        .unwrap_or_else(|| panic!("unsupported GPU {backend} low-precision type {suffix}"));
    format!(
        "extern \"C\" __global__ void {kernel_name}(const {scalar_ty}* xs, {scalar_ty}* ys, size_t n) {{\n    size_t i = (size_t)(blockIdx.x * blockDim.x + threadIdx.x);\n    if (i < n) {{ ys[i] = {eval_name}(xs[i]); }}\n}}"
    )
}

fn render_gpu_backend_file(export: &FunctionRemezExport, backend: &str) -> String {
    let include = gpu_include(backend);
    let mut out = String::new();
    out.push_str("// @generated by aether-fpx-gen\n");
    out.push_str(&format!("// backend: {backend}\n"));
    out.push_str(&format!("// function: {}\n\n", export.function));
    if !include.is_empty() {
        out.push_str(include);
        out.push('\n');
    }
    out.push_str("#include <math.h>\n#include <stddef.h>\n#include <stdint.h>\n\n");
    out.push_str(&gpu_low_precision_prelude(backend));
    out.push_str("\n");

    for (idx, seg) in export.segments.iter().enumerate() {
        let segment_tag = format!("seg{idx}");
        let backend_export = if backend == "cuda" {
            &seg.cuda
        } else {
            &seg.rocm
        };
        out.push_str(&backend_export.coeff_decl_f64);
        out.push('\n');
        out.push_str(&backend_export.coeff_decl_f32);
        out.push_str("\n\n");
        out.push_str(&gpu_eval_source(
            backend,
            &export.function,
            &segment_tag,
            &export.table_domain,
            "double",
            &backend_export.coeff_symbol_f64,
            seg.coeffs_f64.len(),
        ));
        out.push_str("\n\n");
        out.push_str(&gpu_eval_source(
            backend,
            &export.function,
            &segment_tag,
            &export.table_domain,
            "float",
            &backend_export.coeff_symbol_f32,
            seg.coeffs_f32.len(),
        ));
        out.push_str("\n\n");

        let coeff_f16_symbol = gpu_lowp_coeff_symbol(&export.function, idx, "f16_bits");
        let coeff_fp4_symbol = gpu_lowp_coeff_symbol(&export.function, idx, "fp4_e2m1_bits");
        let coeff_fp6_symbol = gpu_lowp_coeff_symbol(&export.function, idx, "fp6_e3m2_bits");
        let coeff_fp8e4_symbol = gpu_lowp_coeff_symbol(&export.function, idx, "fp8_e4m3_bits");
        let coeff_fp8e5_symbol = gpu_lowp_coeff_symbol(&export.function, idx, "fp8_e5m2_bits");

        out.push_str(&coeff_decl_u16(
            &coeff_f16_symbol,
            &seg.coeffs_f16_bits,
            backend,
        ));
        out.push('\n');
        out.push_str(&coeff_decl_u8(
            &coeff_fp4_symbol,
            &seg.coeffs_fp4_e2m1_bits,
            backend,
        ));
        out.push('\n');
        out.push_str(&coeff_decl_u8(
            &coeff_fp6_symbol,
            &seg.coeffs_fp6_e3m2_bits,
            backend,
        ));
        out.push('\n');
        out.push_str(&coeff_decl_u8(
            &coeff_fp8e4_symbol,
            &seg.coeffs_fp8_e4m3_bits,
            backend,
        ));
        out.push('\n');
        out.push_str(&coeff_decl_u8(
            &coeff_fp8e5_symbol,
            &seg.coeffs_fp8_e5m2_bits,
            backend,
        ));
        out.push_str("\n\n");

        out.push_str(&render_gpu_lowp_segment_eval(
            &export.function,
            idx,
            &export.table_domain,
            "f16",
            &coeff_f16_symbol,
            seg.coeffs_f16_bits.len(),
            "fpx_decode_f16_bits",
        ));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_segment_eval(
            &export.function,
            idx,
            &export.table_domain,
            "fp4_e2m1",
            &coeff_fp4_symbol,
            seg.coeffs_fp4_e2m1_bits.len(),
            "fpx_decode_fp4_e2m1",
        ));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_segment_eval(
            &export.function,
            idx,
            &export.table_domain,
            "fp6_e3m2",
            &coeff_fp6_symbol,
            seg.coeffs_fp6_e3m2_bits.len(),
            "fpx_decode_fp6_e3m2",
        ));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_segment_eval(
            &export.function,
            idx,
            &export.table_domain,
            "fp8_e4m3",
            &coeff_fp8e4_symbol,
            seg.coeffs_fp8_e4m3_bits.len(),
            "fpx_decode_fp8_e4m3",
        ));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_segment_eval(
            &export.function,
            idx,
            &export.table_domain,
            "fp8_e5m2",
            &coeff_fp8e5_symbol,
            seg.coeffs_fp8_e5m2_bits.len(),
            "fpx_decode_fp8_e5m2",
        ));
        out.push_str("\n\n");
    }

    out.push_str(&render_gpu_raw_dispatch(export, "double"));
    out.push_str("\n\n");
    out.push_str(&render_gpu_raw_dispatch(export, "float"));
    out.push_str("\n\n");
    out.push_str(&render_gpu_lowp_raw_dispatch(export, "f16"));
    out.push_str("\n\n");
    out.push_str(&render_gpu_lowp_raw_dispatch(export, "fp4_e2m1"));
    out.push_str("\n\n");
    out.push_str(&render_gpu_lowp_raw_dispatch(export, "fp6_e3m2"));
    out.push_str("\n\n");
    out.push_str(&render_gpu_lowp_raw_dispatch(export, "fp8_e4m3"));
    out.push_str("\n\n");
    out.push_str(&render_gpu_lowp_raw_dispatch(export, "fp8_e5m2"));
    out.push_str("\n\n");
    out.push_str(&render_gpu_reduced_domain_accessors(export));
    out.push_str("\n\n");
    out.push_str(&render_gpu_reduced_kernel_entry(export, "double"));
    out.push_str("\n\n");
    out.push_str(&render_gpu_reduced_kernel_entry(export, "float"));
    out.push_str("\n\n");
    out.push_str(&render_gpu_lowp_reduced_kernel_entry(
        export, backend, "f16",
    ));
    out.push_str("\n\n");
    out.push_str(&render_gpu_lowp_reduced_kernel_entry(
        export, backend, "fp4_e2m1",
    ));
    out.push_str("\n\n");
    out.push_str(&render_gpu_lowp_reduced_kernel_entry(
        export, backend, "fp6_e3m2",
    ));
    out.push_str("\n\n");
    out.push_str(&render_gpu_lowp_reduced_kernel_entry(
        export, backend, "fp8_e4m3",
    ));
    out.push_str("\n\n");
    out.push_str(&render_gpu_lowp_reduced_kernel_entry(
        export, backend, "fp8_e5m2",
    ));

    if gpu_has_full_domain_implementation(&export.function) {
        out.push_str("\n\n");
        out.push_str(&render_gpu_public_dispatch(export, "double"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_public_dispatch(export, "float"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_dispatch(export, "f16"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_dispatch(export, "fp4_e2m1"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_dispatch(export, "fp6_e3m2"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_dispatch(export, "fp8_e4m3"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_dispatch(export, "fp8_e5m2"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_kernel_entry(export, "double"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_kernel_entry(export, "float"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_kernel_entry(export, backend, "f16"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_kernel_entry(export, backend, "fp4_e2m1"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_kernel_entry(export, backend, "fp6_e3m2"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_kernel_entry(export, backend, "fp8_e4m3"));
        out.push_str("\n\n");
        out.push_str(&render_gpu_lowp_kernel_entry(export, backend, "fp8_e5m2"));
    }
    out.push('\n');
    out
}

pub fn write_function_gpu_kernels(export: &FunctionRemezExport) {
    use std::fs;

    let generator = common::generator_for(&export.function).unwrap_or_else(|err| panic!("{err}"));
    let base_out = out_dir_for("gpu", "target/gpu");
    for (backend, source) in [
        ("cuda", generator.cuda_generate(export)),
        ("rocm", generator.rocm_generate(export)),
    ] {
        let rel = gpu_rel_path(backend, &export.function);
        let path = base_out.join(&rel);
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).unwrap();
        }
        fs::write(&path, source).unwrap();
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GpuBackendExport {
    pub backend: String,
    pub source_rel_path: String,
    pub coeff_symbol_f64: String,
    pub coeff_symbol_f32: String,
    pub coeff_decl_f64: String,
    pub coeff_decl_f32: String,
    pub evaluator_f64: String,
    pub evaluator_f32: String,
}

fn build_gpu_backend_export(
    backend: &str,
    func_name: &str,
    segment_tag: &str,
    dom: &str,
    coeffs_f64: &[f64],
    coeffs_f32: &[f32],
) -> GpuBackendExport {
    let ident = sanitize_ident(func_name);
    let coeff_symbol_f64 = format!("{}_{}_coeffs_f64", ident, segment_tag);
    let coeff_symbol_f32 = format!("{}_{}_coeffs_f32", ident, segment_tag);

    GpuBackendExport {
        backend: backend.to_string(),
        source_rel_path: gpu_rel_path(backend, func_name),
        coeff_symbol_f64: coeff_symbol_f64.clone(),
        coeff_symbol_f32: coeff_symbol_f32.clone(),
        coeff_decl_f64: coeff_decl_f64(&coeff_symbol_f64, coeffs_f64, backend),
        coeff_decl_f32: coeff_decl_f32(&coeff_symbol_f32, coeffs_f32, backend),
        evaluator_f64: gpu_eval_source(
            backend,
            func_name,
            segment_tag,
            dom,
            "double",
            &coeff_symbol_f64,
            coeffs_f64.len(),
        ),
        evaluator_f32: gpu_eval_source(
            backend,
            func_name,
            segment_tag,
            dom,
            "float",
            &coeff_symbol_f32,
            coeffs_f32.len(),
        ),
    }
}

/* ------------------------------------------------------------------------------------------------
 * Function-level JSON export (this is what generators.rs should build)
 * ------------------------------------------------------------------------------------------------ */

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FunctionRemezSegmentExport {
    pub interval_x: (f64, f64),
    pub interval_solve: (f64, f64),
    pub max_error_f128_abs: f64,
    pub coeffs_f64: Vec<f64>,
    pub coeffs_f32: Vec<f32>,
    pub coeffs_f16_bits: Vec<u16>,
    pub coeffs_fp4_e2m1_bits: Vec<u8>,
    pub coeffs_fp6_e3m2_bits: Vec<u8>,
    pub coeffs_fp8_e4m3_bits: Vec<u8>,
    pub coeffs_fp8_e5m2_bits: Vec<u8>,
    pub cuda: GpuBackendExport,
    pub rocm: GpuBackendExport,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FunctionRemezExport {
    pub function: String,
    /// "x" for PolyDomain::X, "z" for ZSin/ZCos
    pub table_domain: String,
    pub segments: Vec<FunctionRemezSegmentExport>,
}

impl FunctionRemezExport {
    pub fn new(function: String, table_domain: String) -> Self {
        Self {
            function,
            table_domain,
            segments: Vec::new(),
        }
    }

    pub fn add_segment<const N: usize, const K: usize>(
        &mut self,
        rz: &RemezApprox<f128, N, K>,
        interval_x: (f128, f128),
        interval_solve: (f128, f128),
        max_error: f128,
    ) {
        let coeffs_f64: Vec<f64> = rz.coefficients.iter().map(|&c| c as f64).collect();
        let coeffs_f32: Vec<f32> = rz.coefficients.iter().map(|&c| c as f32).collect();
        let coeffs_f16_bits = rz
            .coefficients
            .iter()
            .map(|&c| (c as f32 as f16).to_bits())
            .collect();
        let coeffs_fp4_e2m1_bits = quantize_coeffs_fp4_bits(&coeffs_f32);
        let coeffs_fp6_e3m2_bits = quantize_coeffs_fp6_bits(&coeffs_f32);
        let coeffs_fp8_e4m3_bits = quantize_coeffs_fp8_e4m3_bits(&coeffs_f32);
        let coeffs_fp8_e5m2_bits = quantize_coeffs_fp8_e5m2_bits(&coeffs_f32);
        let seg_idx = self.segments.len();
        let segment_tag = format!("seg{seg_idx}");
        let cuda = build_gpu_backend_export(
            "cuda",
            &self.function,
            &segment_tag,
            &self.table_domain,
            &coeffs_f64,
            &coeffs_f32,
        );
        let rocm = build_gpu_backend_export(
            "rocm",
            &self.function,
            &segment_tag,
            &self.table_domain,
            &coeffs_f64,
            &coeffs_f32,
        );

        self.segments.push(FunctionRemezSegmentExport {
            interval_x: (interval_x.0 as f64, interval_x.1 as f64),
            interval_solve: (interval_solve.0 as f64, interval_solve.1 as f64),
            max_error_f128_abs: max_error as f64,
            coeffs_f64,
            coeffs_f32,
            coeffs_f16_bits,
            coeffs_fp4_e2m1_bits,
            coeffs_fp6_e3m2_bits,
            coeffs_fp8_e4m3_bits,
            coeffs_fp8_e5m2_bits,
            cuda,
            rocm,
        });
    }

    pub fn add_segment_from_plan(&mut self, seg: &SegmentPlan) {
        let coeffs_f64 = seg.c64.as_slice().to_vec();
        let coeffs_f32 = seg.c32.as_slice().to_vec();
        let coeffs_f16_bits = seg.c16.as_slice().iter().map(|c| c.to_bits()).collect();
        let coeffs_fp4_e2m1_bits = quantize_coeffs_fp4_bits(&coeffs_f32);
        let coeffs_fp6_e3m2_bits = quantize_coeffs_fp6_bits(&coeffs_f32);
        let coeffs_fp8_e4m3_bits = quantize_coeffs_fp8_e4m3_bits(&coeffs_f32);
        let coeffs_fp8_e5m2_bits = quantize_coeffs_fp8_e5m2_bits(&coeffs_f32);
        let seg_idx = self.segments.len();
        let segment_tag = format!("seg{seg_idx}");
        let cuda = build_gpu_backend_export(
            "cuda",
            &self.function,
            &segment_tag,
            &self.table_domain,
            &coeffs_f64,
            &coeffs_f32,
        );
        let rocm = build_gpu_backend_export(
            "rocm",
            &self.function,
            &segment_tag,
            &self.table_domain,
            &coeffs_f64,
            &coeffs_f32,
        );
        let dom = poly_domain(&self.function);
        let interval_solve =
            solve_interval_from_x_for_name(&self.function, dom, seg.a128, seg.b128);

        self.segments.push(FunctionRemezSegmentExport {
            interval_x: (seg.a64, seg.b64),
            interval_solve: (interval_solve.0 as f64, interval_solve.1 as f64),
            max_error_f128_abs: f64::NAN,
            coeffs_f64,
            coeffs_f32,
            coeffs_f16_bits,
            coeffs_fp4_e2m1_bits,
            coeffs_fp6_e3m2_bits,
            coeffs_fp8_e4m3_bits,
            coeffs_fp8_e5m2_bits,
            cuda,
            rocm,
        });
    }
}

/* ------------------------------------------------------------------------------------------------
 * JSON exports (segment-level struct kept if you still use it elsewhere)
 * ------------------------------------------------------------------------------------------------ */

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RemezExport<const N: usize, const K: usize> {
    pub function: String,
    pub interval_x: (f64, f64),
    pub interval_solve: (f64, f64),
    pub table_domain: String,
    pub max_error_f128_abs: f64,
    pub coeffs_f64: Vec<f64>,
    pub coeffs_f32: Vec<f32>,
    pub coeffs_f16_bits: Vec<u16>,
    pub coeffs_fp4_e2m1_bits: Vec<u8>,
    pub coeffs_fp6_e3m2_bits: Vec<u8>,
    pub coeffs_fp8_e4m3_bits: Vec<u8>,
    pub coeffs_fp8_e5m2_bits: Vec<u8>,
    pub cuda: GpuBackendExport,
    pub rocm: GpuBackendExport,
}

impl<const N: usize, const K: usize> RemezExport<N, K> {
    pub fn from_remez<F>(
        rz: &RemezApprox<f128, N, K>,
        func_name: String,
        _f_x: F,
        interval_x: (f128, f128),
        interval_solve: (f128, f128),
        max_error: f128,
        dom: PolyDomain,
    ) -> Self
    where
        F: Fn(f128) -> f128,
    {
        let coeffs_f64: Vec<f64> = rz.coefficients.iter().map(|&c| c as f64).collect();
        let coeffs_f32: Vec<f32> = rz.coefficients.iter().map(|&c| c as f32).collect();
        let coeffs_f16_bits = rz
            .coefficients
            .iter()
            .map(|&c| (c as f32 as f16).to_bits())
            .collect();
        let coeffs_fp4_e2m1_bits = quantize_coeffs_fp4_bits(&coeffs_f32);
        let coeffs_fp6_e3m2_bits = quantize_coeffs_fp6_bits(&coeffs_f32);
        let coeffs_fp8_e4m3_bits = quantize_coeffs_fp8_e4m3_bits(&coeffs_f32);
        let coeffs_fp8_e5m2_bits = quantize_coeffs_fp8_e5m2_bits(&coeffs_f32);

        let table_domain = dom.table_domain_str().to_string();
        let cuda = build_gpu_backend_export(
            "cuda",
            &func_name,
            "seg0",
            &table_domain,
            &coeffs_f64,
            &coeffs_f32,
        );
        let rocm = build_gpu_backend_export(
            "rocm",
            &func_name,
            "seg0",
            &table_domain,
            &coeffs_f64,
            &coeffs_f32,
        );

        Self {
            function: func_name,
            interval_x: (interval_x.0 as f64, interval_x.1 as f64),
            interval_solve: (interval_solve.0 as f64, interval_solve.1 as f64),
            table_domain,
            max_error_f128_abs: max_error as f64,
            coeffs_f64,
            coeffs_f32,
            coeffs_f16_bits,
            coeffs_fp4_e2m1_bits,
            coeffs_fp6_e3m2_bits,
            coeffs_fp8_e4m3_bits,
            coeffs_fp8_e5m2_bits,
            cuda,
            rocm,
        }
    }
}

/* ------------------------------------------------------------------------------------------------
 * Rust codegen
 * ------------------------------------------------------------------------------------------------ */
pub fn write_function_to_rs(func_name: &str, func: &FunctionSpec, plan: &FunctionPlan) {
    common::generator_for(func_name)
        .unwrap_or_else(|err| panic!("{err}"))
        .rust_generate(func_name, func, plan);
}

fn write_function_to_rs_shared(func_name: &str, _func: &FunctionSpec, plan: &FunctionPlan) {
    use std::io::BufWriter;
    use std::{fs, io::Write};

    let base_out = out_dir_for("rs", "target/rs");
    fs::create_dir_all(&base_out).unwrap();
    let path = base_out.join(format!("{}.rs", plan.name));
    let mut w = BufWriter::new(fs::File::create(&path).unwrap());

    let helpers = common::generator_for(func_name)
        .unwrap_or_else(|err| panic!("{err}"))
        .helpers();
    let ident = helpers.ident();
    let upper = helpers.upper();

    let is_sin = func_name == "sin";
    let is_cos = func_name == "cos";
    let is_tan = func_name == "tan";
    let is_exp = func_name == "exp";
    let is_expm1 = func_name == "expm1";
    let is_log = func_name == "log";
    let is_log1p = func_name == "log1p";
    let is_simple = matches!(
        func_name,
        "sinh" | "cosh" | "tanh" | "asin" | "acos" | "atan" | "log" | "log1p" | "exp" | "expm1"
    );
    let ln_2_hi_f128_mask = 0xffff_ffff_ffff_0000_0000_0000_0000_0000u128;
    let ln_2_hi_f128 = f128::from_bits(core::f128::consts::LN_2.to_bits() & ln_2_hi_f128_mask);
    let ln_2_rem_f128 = core::f128::consts::LN_2 - ln_2_hi_f128;
    // Each split component is an independently encoded floating-point number.
    // Preserve the remainder's exponent while truncating its significand; masking
    // the original bit positions would turn the middle component into a subnormal.
    let ln_2_mid_f128 = f128::from_bits(ln_2_rem_f128.to_bits() & ln_2_hi_f128_mask);
    let ln_2_lo_f128 = core::f128::consts::LN_2 - ln_2_hi_f128 - ln_2_mid_f128;

    if !(is_sin || is_cos || is_tan || is_simple) {
        panic!("unsupported function {}", func_name);
    }

    // ---------------------------------------------------------------------------------------------
    // Module header
    // ---------------------------------------------------------------------------------------------
    writeln!(w, "// @generated by aether-fpx-gen").unwrap();

    if is_tan {
        writeln!(
            w,
            "use crate::policy::policy::{{
                apply_policy_tan_f64,
                apply_policy_tan_f32
            }};"
        )
        .unwrap();

        writeln!(w, "use core::f64::consts::FRAC_PI_2 as FRAC_PI_2_F64;").unwrap();
        writeln!(w, "use core::f32::consts::FRAC_PI_2 as FRAC_PI_2_F32;").unwrap();
        writeln!(w, "use core::f128::consts::FRAC_PI_2 as FRAC_PI_2_F128;").unwrap();

        writeln!(
            w,
            "use crate::trig::payne_hanek::{{rem_pio2_tan_f64, rem_pio2_tan_f32, rem_pio2_tan_f128, rem_pio2_tan_f16}};"
        )
        .unwrap();
    } else if is_sin || is_cos {
        writeln!(w, "use core::f64::consts::FRAC_PI_4 as FRAC_PI_4_F64;").unwrap();
        writeln!(w, "use core::f32::consts::FRAC_PI_4 as FRAC_PI_4_F32;").unwrap();
        writeln!(w, "use core::f128::consts::FRAC_PI_4 as FRAC_PI_4_F128;").unwrap();
        writeln!(
            w,
            "use crate::trig::payne_hanek::{{rem_pio2_f64, rem_pio2_f32, rem_pio2_f128, rem_pio2_f16}};"
        )
        .unwrap();
    }

    writeln!(w, "use core::f16;").unwrap();
    writeln!(w, "use core::f128;").unwrap();
    writeln!(w, "#[cfg(feature = \"mx\")] use crate::ocp::fp4::Fp4E2M1;").unwrap();
    writeln!(w, "#[cfg(feature = \"mx\")] use crate::ocp::fp6::Fp6E3M2;").unwrap();
    writeln!(
        w,
        "#[cfg(feature = \"mx\")] use crate::ocp::fp8::{{Fp8E4M3, Fp8E5M2}};"
    )
    .unwrap();

    if is_exp || is_expm1 || is_log || is_log1p {
        writeln!(w, "const LN_2_F64: f64 = core::f64::consts::LN_2;").unwrap();
        writeln!(w, "const LN_2_F32: f32 = core::f32::consts::LN_2;").unwrap();
        writeln!(w, "const LN_2_F128: f128 = core::f128::consts::LN_2;").unwrap();
        writeln!(w, "const LN_2_HI_F64: f64 = 6.93147180369123816490e-1_f64;").unwrap();
        writeln!(
            w,
            "const LN_2_LO_F64: f64 = 1.90821492927058770002e-10_f64;"
        )
        .unwrap();
        writeln!(
            w,
            "const LN_2_HI_F128: f128 = f128::from_bits({:#034x}u128);",
            ln_2_hi_f128.to_bits()
        )
        .unwrap();
        writeln!(
            w,
            "const LN_2_MID_F128: f128 = f128::from_bits({:#034x}u128);",
            ln_2_mid_f128.to_bits()
        )
        .unwrap();
        writeln!(
            w,
            "const LN_2_LO_F128: f128 = f128::from_bits({:#034x}u128);",
            ln_2_lo_f128.to_bits()
        )
        .unwrap();
        writeln!(w).unwrap();
    }

    if is_exp || is_expm1 {
        writeln!(w, "const LOG2_E_F64: f64 = core::f64::consts::LOG2_E;").unwrap();
        writeln!(w, "const LOG2_E_F32: f32 = core::f32::consts::LOG2_E;").unwrap();
        writeln!(w, "const LOG2_E_F128: f128 = core::f128::consts::LOG2_E;").unwrap();
        writeln!(w, "const EXP_MAX_F64: f64 = 1024.0_f64 * LN_2_F64;").unwrap();
        writeln!(w, "const EXP_MIN_F64: f64 = -1075.0_f64 * LN_2_F64;").unwrap();
        writeln!(w, "const EXP_MAX_F32: f32 = 128.0_f32 * LN_2_F32;").unwrap();
        writeln!(w, "const EXP_MIN_F32: f32 = -150.0_f32 * LN_2_F32;").unwrap();
        writeln!(w, "const EXP_MAX_F128: f128 = 16384.0_f128 * LN_2_F128;").unwrap();
        writeln!(w, "const EXP_MIN_F128: f128 = -16495.0_f128 * LN_2_F128;").unwrap();
        writeln!(w).unwrap();
    }

    if is_log || is_log1p {
        writeln!(
            w,
            "const SQRT_HALF_F64: f64 = core::f64::consts::FRAC_1_SQRT_2;"
        )
        .unwrap();
        writeln!(
            w,
            "const SQRT_HALF_F32: f32 = core::f32::consts::FRAC_1_SQRT_2;"
        )
        .unwrap();
        writeln!(
            w,
            "const SQRT_HALF_F128: f128 = core::f128::consts::FRAC_1_SQRT_2;"
        )
        .unwrap();
        writeln!(w, "const SQRT_2_F64: f64 = core::f64::consts::SQRT_2;").unwrap();
        writeln!(w, "const SQRT_2_F32: f32 = core::f32::consts::SQRT_2;").unwrap();
        writeln!(w, "const SQRT_2_F128: f128 = core::f128::consts::SQRT_2;").unwrap();
        writeln!(w).unwrap();
    }

    writeln!(w).unwrap();

    // ---------------------------------------------------------------------------------------------
    // Interval constants
    // ---------------------------------------------------------------------------------------------
    if is_sin || is_cos {
        writeln!(w, "pub const {upper}_INTERVAL_MIN_F64: f64 = 0.0_f64;").unwrap();
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MAX_F64: f64 = FRAC_PI_4_F64;"
        )
        .unwrap();
        writeln!(w, "pub const {upper}_INTERVAL_MIN_F32: f32 = 0.0_f32;").unwrap();
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MAX_F32: f32 = FRAC_PI_4_F32;"
        )
        .unwrap();
        writeln!(w, "pub const {upper}_INTERVAL_MIN_F128: f128 = 0.0_f128;").unwrap();
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MAX_F128: f128 = FRAC_PI_4_F128;"
        )
        .unwrap();
    } else if is_tan {
        writeln!(w, "pub const {upper}_INTERVAL_MIN_F64: f64 = 0.0_f64;").unwrap();
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MAX_F64: f64 = 0.5_f64 * FRAC_PI_2_F64;"
        )
        .unwrap();
        writeln!(w, "pub const {upper}_INTERVAL_MIN_F32: f32 = 0.0_f32;").unwrap();
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MAX_F32: f32 = 0.5_f32 * FRAC_PI_2_F32;"
        )
        .unwrap();
        writeln!(w, "pub const {upper}_INTERVAL_MIN_F128: f128 = 0.0_f128;").unwrap();
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MAX_F128: f128 = 0.5_f128 * FRAC_PI_2_F128;"
        )
        .unwrap();
    } else {
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MIN_F64: f64 = {:17e}_f64;",
            plan.interval_min_f64
        )
        .unwrap();
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MAX_F64: f64 = {:17e}_f64;",
            plan.interval_max_f64
        )
        .unwrap();
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MIN_F32: f32 = {:9e}_f32;",
            plan.interval_min_f32
        )
        .unwrap();
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MAX_F32: f32 = {:9e}_f32;",
            plan.interval_max_f32
        )
        .unwrap();
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MIN_F128: f128 = f128::from_bits({:#034x}u128);",
            plan.interval_min_f128.to_bits()
        )
        .unwrap();
        writeln!(
            w,
            "pub const {upper}_INTERVAL_MAX_F128: f128 = f128::from_bits({:#034x}u128);",
            plan.interval_max_f128.to_bits()
        )
        .unwrap();
    }

    writeln!(w).unwrap();

    // ---------------------------------------------------------------------------------------------
    // Coefficients
    // ---------------------------------------------------------------------------------------------
    if plan.segments.len() == 1 {
        let s = &plan.segments[0];

        write!(
            w,
            "pub const COEFFS128: [f128; {}] = ",
            s.c128.as_slice().len()
        )
        .unwrap();
        write_coeffset128_inline(&mut w, &s.c128).unwrap();
        writeln!(w, ";").unwrap();

        write!(
            w,
            "pub const COEFFS64: [f64; {}] = ",
            s.c64.as_slice().len()
        )
        .unwrap();
        write_coeffset64_inline(&mut w, &s.c64).unwrap();
        writeln!(w, ";").unwrap();

        write!(
            w,
            "pub const COEFFS32: [f32; {}] = ",
            s.c32.as_slice().len()
        )
        .unwrap();
        write_coeffset32_inline(&mut w, &s.c32).unwrap();
        writeln!(w, ";").unwrap();

        if !(is_log || is_log1p) {
            write!(
                w,
                "pub const COEFFS16: [f16; {}] = ",
                s.c16.as_slice().len()
            )
            .unwrap();
            write_coeffset16_inline(&mut w, &s.c16).unwrap();
            writeln!(w, ";").unwrap();
        }

        writeln!(w).unwrap();
    } else {
        for (idx, s) in plan.segments.iter().enumerate() {
            write!(
                w,
                "pub const COEFFS128_SEG{idx}: [f128; {}] = ",
                s.c128.as_slice().len()
            )
            .unwrap();
            write_coeffset128_inline(&mut w, &s.c128).unwrap();
            writeln!(w, ";").unwrap();

            write!(
                w,
                "pub const COEFFS64_SEG{idx}: [f64; {}] = ",
                s.c64.as_slice().len()
            )
            .unwrap();
            write_coeffset64_inline(&mut w, &s.c64).unwrap();
            writeln!(w, ";").unwrap();

            write!(
                w,
                "pub const COEFFS32_SEG{idx}: [f32; {}] = ",
                s.c32.as_slice().len()
            )
            .unwrap();
            write_coeffset32_inline(&mut w, &s.c32).unwrap();
            writeln!(w, ";").unwrap();

            if !(is_log || is_log1p) {
                write!(
                    w,
                    "pub const COEFFS16_SEG{idx}: [f16; {}] = ",
                    s.c16.as_slice().len()
                )
                .unwrap();
                write_coeffset16_inline(&mut w, &s.c16).unwrap();
                writeln!(w, ";").unwrap();
            }
        }

        write!(
            w,
            "pub const {upper}_SEG_BOUNDS_F64: [f64; {}] = [",
            plan.segments.len()
        )
        .unwrap();
        for (idx, s) in plan.segments.iter().enumerate() {
            if idx != 0 {
                write!(w, ", ").unwrap();
            }
            write!(w, "{:17e}_f64", s.b64).unwrap();
        }
        writeln!(w, "];").unwrap();

        write!(
            w,
            "pub const {upper}_SEG_BOUNDS_F32: [f32; {}] = [",
            plan.segments.len()
        )
        .unwrap();
        for (idx, s) in plan.segments.iter().enumerate() {
            if idx != 0 {
                write!(w, ", ").unwrap();
            }
            write!(w, "{:9e}_f32", s.b32).unwrap();
        }
        writeln!(w, "];").unwrap();

        write!(
            w,
            "pub const {upper}_SEG_BOUNDS_F128: [f128; {}] = [",
            plan.segments.len()
        )
        .unwrap();
        for (idx, s) in plan.segments.iter().enumerate() {
            if idx != 0 {
                write!(w, ", ").unwrap();
            }
            write!(w, "f128::from_bits({:#034x}u128)", s.b128.to_bits()).unwrap();
        }
        writeln!(w, "];").unwrap();

        writeln!(w).unwrap();
    }

    write_mx_coefficients(&mut w, plan).unwrap();

    // ---------------------------------------------------------------------------------------------
    // Cross-kernel imports for sin/cos
    // ---------------------------------------------------------------------------------------------
    write_mx_trig_imports(&mut w, is_sin, is_cos).unwrap();

    if is_sin || is_cos || is_tan {
        writeln!(
            w,
            "use crate::common::{{eval_estrin_f32_slice, eval_horner_f64_slice, eval_horner_f32_slice, eval_horner_f16_slice, eval_horner_f128_slice}};"
        ).unwrap();
    } else {
        writeln!(
            w,
            "use crate::common::{{eval_horner_f64_slice, eval_horner_f32_slice, eval_horner_f16_slice, eval_horner_f128_slice}};"
        ).unwrap();
    }
    writeln!(w).unwrap();

    write_mx_horner_helpers(&mut w).unwrap();

    let generator = common::generator_for(func_name).unwrap_or_else(|err| panic!("{err}"));
    generator
        .write_kernel_helpers(&mut w, plan, &upper)
        .unwrap();

    write_mx_kernel_helpers(&mut w, &upper, func_name, plan, is_sin, is_cos, is_tan).unwrap();

    generator.write_helpers(&mut w).unwrap();

    let public_ctx = common::PublicFnCtx {
        ident: &ident,
        upper: &upper,
    };

    // ---------------------------------------------------------------------------------------------
    // Public f64
    // ---------------------------------------------------------------------------------------------
    writeln!(w, "#[inline]").unwrap();
    writeln!(w, "pub fn {ident}(x: f64) -> f64 {{").unwrap();
    common::write_scalar_public_body(&mut w, func_name, &public_ctx, common::ScalarKind::F64)
        .unwrap();

    writeln!(w, "}}").unwrap();

    // ---------------------------------------------------------------------------------------------
    // Public f128
    // ---------------------------------------------------------------------------------------------
    writeln!(w).unwrap();
    writeln!(w, "#[inline]").unwrap();
    writeln!(w, "pub fn {ident}l(x: f128) -> f128 {{").unwrap();
    common::write_scalar_public_body(&mut w, func_name, &public_ctx, common::ScalarKind::F128)
        .unwrap();

    writeln!(w, "}}").unwrap();

    // ---------------------------------------------------------------------------------------------
    // Public f32
    // ---------------------------------------------------------------------------------------------
    writeln!(w).unwrap();
    writeln!(w, "#[inline]").unwrap();
    writeln!(w, "pub fn {ident}f(x: f32) -> f32 {{").unwrap();
    common::write_scalar_public_body(&mut w, func_name, &public_ctx, common::ScalarKind::F32)
        .unwrap();

    writeln!(w, "}}").unwrap();

    // ---------------------------------------------------------------------------------------------
    // Public f16
    // ---------------------------------------------------------------------------------------------
    writeln!(w).unwrap();
    writeln!(w, "#[inline]").unwrap();
    writeln!(w, "pub fn {ident}s(x: f16) -> f16 {{").unwrap();
    common::write_scalar_public_body(&mut w, func_name, &public_ctx, common::ScalarKind::F16)
        .unwrap();

    writeln!(w, "}}").unwrap();

    write_mx_public_functions(&mut w, func_name, &ident, &upper).unwrap();

    println!("Wrote Rust function module {}", path.display());
}

fn build_straight_estrin_terms(
    input_var: &str,
    coeff_symbol: &str,
    coeff_len: usize,
    prefix: &str,
) -> (Vec<String>, String) {
    let mut terms = (0..coeff_len)
        .map(|idx| format!("{coeff_symbol}[{idx}]"))
        .collect::<Vec<_>>();
    let mut lines = Vec::new();
    let mut power_expr = input_var.to_string();
    let mut level = 0usize;

    while terms.len() > 1 {
        let mut next_terms = Vec::new();
        let mut term_idx = 0usize;
        let mut pair_idx = 0usize;

        while term_idx + 1 < terms.len() {
            let var = format!("{prefix}_{level}_{pair_idx}");
            lines.push(format!(
                "let {var} = {} + {} * {};",
                terms[term_idx],
                terms[term_idx + 1],
                power_expr
            ));
            next_terms.push(var);
            term_idx += 2;
            pair_idx += 1;
        }

        if term_idx < terms.len() {
            next_terms.push(terms[term_idx].clone());
        }

        terms = next_terms;
        if terms.len() > 1 {
            let next_power = format!("{prefix}_pow_{}", level + 1);
            lines.push(format!("let {next_power} = {power_expr} * {power_expr};"));
            power_expr = next_power;
        }

        level += 1;
    }

    (
        lines,
        terms
            .pop()
            .expect("Estrin evaluation requires coefficients"),
    )
}

fn write_straight_estrin_z_f64<W: std::io::Write>(
    w: &mut W,
    indent: &str,
    input_var: &str,
    coeff_symbol: &str,
    coeff_len: usize,
    kernel: TrigKernel,
) -> Result<(), std::io::Error> {
    writeln!(w, "{indent}let z = {input_var} * {input_var};")?;
    let (lines, expr) = build_straight_estrin_terms("z", coeff_symbol, coeff_len, "estrin");
    for line in lines {
        writeln!(w, "{indent}{line}")?;
    }

    let out = match kernel {
        TrigKernel::Sin => format!("{input_var} * {expr}"),
        TrigKernel::Cos => expr,
        TrigKernel::Tan => format!("{input_var} + {input_var} * z * {expr}"),
    };
    writeln!(w, "{indent}{out}")?;
    Ok(())
}

fn write_sin_cos_kernel_helpers<W: std::io::Write>(
    w: &mut W,
    plan: &FunctionPlan,
    upper: &str,
    is_sin: bool,
    is_cos: bool,
) -> Result<(), std::io::Error> {
    if is_sin {
        if plan.segments.len() == 1 {
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_f64(x: f64) -> f64 {{ let z = x * x; x * eval_horner_f64_slice(z, &COEFFS64) }}"
            )?;
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_estrin_f64(x: f64) -> f64 {{"
            )?;
            write_straight_estrin_z_f64(
                w,
                "    ",
                "x",
                "COEFFS64",
                plan.segments[0].c64.as_slice().len(),
                TrigKernel::Sin,
            )?;
            writeln!(w, "}}")?;
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_f32(x: f32) -> f32 {{ let z = x * x; x * eval_horner_f32_slice(z, &COEFFS32) }}"
            )?;
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_estrin_f32(x: f32) -> f32 {{ let z = x * x; x * eval_estrin_f32_slice(z, &COEFFS32) }}"
            )?;
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_f16(x: f16) -> f16 {{ let z = x * x; x * eval_horner_f16_slice(z, &COEFFS16) }}"
            )?;
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_f128(x: f128) -> f128 {{ let z = x * x; x * eval_horner_f128_slice(z, &COEFFS128) }}"
            )?;
        } else {
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_f64(x: f64) -> f64 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{ let z = x * x; x * eval_horner_f64_slice(z, &COEFFS64_SEG{idx}) }}")?;
            }
            let last = plan.segments.len() - 1;
            writeln!(
                w,
                "    else {{ let z = x * x; x * eval_horner_f64_slice(z, &COEFFS64_SEG{last}) }}"
            )?;
            writeln!(w, "}}")?;

            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_estrin_f64(x: f64) -> f64 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{")?;
                write_straight_estrin_z_f64(
                    w,
                    "        ",
                    "x",
                    &format!("COEFFS64_SEG{idx}"),
                    plan.segments[idx].c64.as_slice().len(),
                    TrigKernel::Sin,
                )?;
                writeln!(w, "    }}")?;
            }
            writeln!(w, "    else {{")?;
            write_straight_estrin_z_f64(
                w,
                "        ",
                "x",
                &format!("COEFFS64_SEG{last}"),
                plan.segments[last].c64.as_slice().len(),
                TrigKernel::Sin,
            )?;
            writeln!(w, "    }}")?;
            writeln!(w, "}}")?;

            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_f32(x: f32) -> f32 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = x * x; x * eval_horner_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
            }
            writeln!(
                w,
                "    else {{ let z = x * x; x * eval_horner_f32_slice(z, &COEFFS32_SEG{last}) }}"
            )?;
            writeln!(w, "}}")?;

            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_estrin_f32(x: f32) -> f32 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = x * x; x * eval_estrin_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
            }
            writeln!(
                w,
                "    else {{ let z = x * x; x * eval_estrin_f32_slice(z, &COEFFS32_SEG{last}) }}"
            )?;
            writeln!(w, "}}")?;

            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_f16(x: f16) -> f16 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= ({upper}_SEG_BOUNDS_F32[{idx}] as f16) {{ let z = x * x; x * eval_horner_f16_slice(z, &COEFFS16_SEG{idx}) }}")?;
            }
            let last = plan.segments.len() - 1;
            writeln!(
                w,
                "    else {{ let z = x * x; x * eval_horner_f16_slice(z, &COEFFS16_SEG{last}) }}"
            )?;
            writeln!(w, "}}")?;

            writeln!(
                w,
                "#[inline(always)] pub(crate) fn sin_kernel_f128(x: f128) -> f128 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F128[{idx}] {{ let z = x * x; x * eval_horner_f128_slice(z, &COEFFS128_SEG{idx}) }}")?;
            }
            let last = plan.segments.len() - 1;
            writeln!(
                w,
                "    else {{ let z = x * x; x * eval_horner_f128_slice(z, &COEFFS128_SEG{last}) }}"
            )?;
            writeln!(w, "}}")?;
        }
    } else if is_cos {
        if plan.segments.len() == 1 {
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_f64(x: f64) -> f64 {{ let z = x * x; eval_horner_f64_slice(z, &COEFFS64) }}"
            )?;
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_estrin_f64(x: f64) -> f64 {{"
            )?;
            write_straight_estrin_z_f64(
                w,
                "    ",
                "x",
                "COEFFS64",
                plan.segments[0].c64.as_slice().len(),
                TrigKernel::Cos,
            )?;
            writeln!(w, "}}")?;
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_f32(x: f32) -> f32 {{ let z = x * x; eval_horner_f32_slice(z, &COEFFS32) }}"
            )?;
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_estrin_f32(x: f32) -> f32 {{ let z = x * x; eval_estrin_f32_slice(z, &COEFFS32) }}"
            )?;
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_f16(x: f16) -> f16 {{ let z = x * x; eval_horner_f16_slice(z, &COEFFS16) }}"
            )?;
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_f128(x: f128) -> f128 {{ let z = x * x; eval_horner_f128_slice(z, &COEFFS128) }}"
            )?;
        } else {
            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_f64(x: f64) -> f64 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{ let z = x * x; eval_horner_f64_slice(z, &COEFFS64_SEG{idx}) }}")?;
            }
            let last = plan.segments.len() - 1;
            writeln!(
                w,
                "    else {{ let z = x * x; eval_horner_f64_slice(z, &COEFFS64_SEG{last}) }}"
            )?;
            writeln!(w, "}}")?;

            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_estrin_f64(x: f64) -> f64 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{")?;
                write_straight_estrin_z_f64(
                    w,
                    "        ",
                    "x",
                    &format!("COEFFS64_SEG{idx}"),
                    plan.segments[idx].c64.as_slice().len(),
                    TrigKernel::Cos,
                )?;
                writeln!(w, "    }}")?;
            }
            writeln!(w, "    else {{")?;
            write_straight_estrin_z_f64(
                w,
                "        ",
                "x",
                &format!("COEFFS64_SEG{last}"),
                plan.segments[last].c64.as_slice().len(),
                TrigKernel::Cos,
            )?;
            writeln!(w, "    }}")?;
            writeln!(w, "}}")?;

            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_f32(x: f32) -> f32 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = x * x; eval_horner_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
            }
            writeln!(
                w,
                "    else {{ let z = x * x; eval_horner_f32_slice(z, &COEFFS32_SEG{last}) }}"
            )?;
            writeln!(w, "}}")?;

            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_estrin_f32(x: f32) -> f32 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = x * x; eval_estrin_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
            }
            writeln!(
                w,
                "    else {{ let z = x * x; eval_estrin_f32_slice(z, &COEFFS32_SEG{last}) }}"
            )?;
            writeln!(w, "}}")?;

            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_f16(x: f16) -> f16 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= ({upper}_SEG_BOUNDS_F32[{idx}] as f16) {{ let z = x * x; eval_horner_f16_slice(z, &COEFFS16_SEG{idx}) }}")?;
            }
            let last = plan.segments.len() - 1;
            writeln!(
                w,
                "    else {{ let z = x * x; eval_horner_f16_slice(z, &COEFFS16_SEG{last}) }}"
            )?;
            writeln!(w, "}}")?;

            writeln!(
                w,
                "#[inline(always)] pub(crate) fn cos_kernel_f128(x: f128) -> f128 {{"
            )?;
            writeln!(w, "    let ax = x.abs();")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F128[{idx}] {{ let z = x * x; eval_horner_f128_slice(z, &COEFFS128_SEG{idx}) }}")?;
            }
            let last = plan.segments.len() - 1;
            writeln!(
                w,
                "    else {{ let z = x * x; eval_horner_f128_slice(z, &COEFFS128_SEG{last}) }}"
            )?;
            writeln!(w, "}}")?;
        }
    }

    writeln!(w, "#[inline(always)] fn k_sin_f64(x: f64) -> f64 {{ crate::trig::sin::sin_kernel_estrin_f64(x) }}")?;
    writeln!(w, "#[inline(always)] fn k_cos_f64(x: f64) -> f64 {{ crate::trig::cos::cos_kernel_estrin_f64(x) }}")?;
    writeln!(
        w,
        "#[inline(always)] fn k_sin_f32(x: f32) -> f32 {{ crate::trig::sin::sin_kernel_f32(x) }}"
    )?;
    writeln!(
        w,
        "#[inline(always)] fn k_cos_f32(x: f32) -> f32 {{ crate::trig::cos::cos_kernel_f32(x) }}"
    )?;
    writeln!(
        w,
        "#[inline(always)] fn k_sin_f16(x: f16) -> f16 {{ crate::trig::sin::sin_kernel_f16(x) }}"
    )?;
    writeln!(
        w,
        "#[inline(always)] fn k_cos_f16(x: f16) -> f16 {{ crate::trig::cos::cos_kernel_f16(x) }}"
    )?;
    writeln!(w, "#[inline(always)] fn k_sin_f128(x: f128) -> f128 {{ crate::trig::sin::sin_kernel_f128(x) }}")?;
    writeln!(w, "#[inline(always)] fn k_cos_f128(x: f128) -> f128 {{ crate::trig::cos::cos_kernel_f128(x) }}")?;
    writeln!(w)?;

    Ok(())
}

fn write_tan_kernel_helpers<W: std::io::Write>(
    w: &mut W,
    plan: &FunctionPlan,
    upper: &str,
) -> std::io::Result<()> {
    if plan.segments.len() == 1 {
        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_f64(x: f64) -> f64 {{ let ax = x.abs(); let z = ax * ax; let t = ax + ax * z * eval_horner_f64_slice(z, &COEFFS64); if x.is_sign_negative() {{ -t }} else {{ t }} }}"
        )?;
        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_estrin_f64(x: f64) -> f64 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        writeln!(w, "    let t = {{")?;
        write_straight_estrin_z_f64(
            w,
            "        ",
            "ax",
            "COEFFS64",
            plan.segments[0].c64.as_slice().len(),
            TrigKernel::Tan,
        )?;
        writeln!(w, "    }};")?;
        writeln!(w, "    if x.is_sign_negative() {{ -t }} else {{ t }}")?;
        writeln!(w, "}}")?;
        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_f32(x: f32) -> f32 {{ let ax = x.abs(); let z = ax * ax; let t = ax + ax * z * eval_horner_f32_slice(z, &COEFFS32); if x.is_sign_negative() {{ -t }} else {{ t }} }}"
        )?;
        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_estrin_f32(x: f32) -> f32 {{ let ax = x.abs(); let z = ax * ax; let t = ax + ax * z * eval_estrin_f32_slice(z, &COEFFS32); if x.is_sign_negative() {{ -t }} else {{ t }} }}"
        )?;
        writeln!(
            w,
            "#[inline(always)] fn k_tan_f64(x: f64) -> f64 {{ tan_kernel_estrin_f64(x) }}"
        )?;
        writeln!(
            w,
            "#[inline(always)] fn k_tan_f32(x: f32) -> f32 {{ tan_kernel_f32(x) }}"
        )?;
        writeln!(w, "#[inline(always)] fn k_tan_f16(x: f16) -> f16 {{ let z = x * x; x + x * z * eval_horner_f16_slice(z, &COEFFS16) }}")?;
        writeln!(w, "#[inline(always)] fn k_tan_f128(x: f128) -> f128 {{ let z = x * x; x + x * z * eval_horner_f128_slice(z, &COEFFS128) }}")?;
        writeln!(w)?;
    } else {
        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_f64(x: f64) -> f64 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        writeln!(w, "    let t = if ax <= {upper}_SEG_BOUNDS_F64[0] {{ let z = ax * ax; ax + ax * z * eval_horner_f64_slice(z, &COEFFS64_SEG0) }}")?;
        for idx in 1..plan.segments.len() {
            writeln!(w, "    else if ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{ let z = ax * ax; ax + ax * z * eval_horner_f64_slice(z, &COEFFS64_SEG{idx}) }}")?;
        }
        let last = plan.segments.len() - 1;
        writeln!(w, "    else {{ let z = ax * ax; ax + ax * z * eval_horner_f64_slice(z, &COEFFS64_SEG{last}) }};")?;
        writeln!(w, "    if x.is_sign_negative() {{ -t }} else {{ t }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_estrin_f64(x: f64) -> f64 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        writeln!(w, "    let t = if ax <= {upper}_SEG_BOUNDS_F64[0] {{")?;
        write_straight_estrin_z_f64(
            w,
            "        ",
            "ax",
            "COEFFS64_SEG0",
            plan.segments[0].c64.as_slice().len(),
            TrigKernel::Tan,
        )?;
        writeln!(w, "    }}")?;
        for idx in 1..plan.segments.len() {
            writeln!(w, "    else if ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{")?;
            write_straight_estrin_z_f64(
                w,
                "        ",
                "ax",
                &format!("COEFFS64_SEG{idx}"),
                plan.segments[idx].c64.as_slice().len(),
                TrigKernel::Tan,
            )?;
            writeln!(w, "    }}")?;
        }
        writeln!(w, "    else {{")?;
        write_straight_estrin_z_f64(
            w,
            "        ",
            "ax",
            &format!("COEFFS64_SEG{last}"),
            plan.segments[last].c64.as_slice().len(),
            TrigKernel::Tan,
        )?;
        writeln!(w, "    }};")?;
        writeln!(w, "    if x.is_sign_negative() {{ -t }} else {{ t }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_f32(x: f32) -> f32 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        writeln!(w, "    let t = if ax <= {upper}_SEG_BOUNDS_F32[0] {{ let z = ax * ax; ax + ax * z * eval_horner_f32_slice(z, &COEFFS32_SEG0) }}")?;
        for idx in 1..plan.segments.len() {
            writeln!(w, "    else if ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = ax * ax; ax + ax * z * eval_horner_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
        }
        writeln!(w, "    else {{ let z = ax * ax; ax + ax * z * eval_horner_f32_slice(z, &COEFFS32_SEG{last}) }};")?;
        writeln!(w, "    if x.is_sign_negative() {{ -t }} else {{ t }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_estrin_f32(x: f32) -> f32 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        writeln!(w, "    let t = if ax <= {upper}_SEG_BOUNDS_F32[0] {{ let z = ax * ax; ax + ax * z * eval_estrin_f32_slice(z, &COEFFS32_SEG0) }}")?;
        for idx in 1..plan.segments.len() {
            writeln!(w, "    else if ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = ax * ax; ax + ax * z * eval_estrin_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
        }
        writeln!(w, "    else {{ let z = ax * ax; ax + ax * z * eval_estrin_f32_slice(z, &COEFFS32_SEG{last}) }};")?;
        writeln!(w, "    if x.is_sign_negative() {{ -t }} else {{ t }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;

        writeln!(
            w,
            "#[inline(always)] fn k_tan_f64(x: f64) -> f64 {{ tan_kernel_estrin_f64(x) }}"
        )?;
        writeln!(
            w,
            "#[inline(always)] fn k_tan_f32(x: f32) -> f32 {{ tan_kernel_f32(x) }}"
        )?;

        writeln!(w, "#[inline(always)] fn k_tan_f16(x: f16) -> f16 {{")?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(
            w,
            "    {kw} ax <= ({upper}_SEG_BOUNDS_F32[{idx}] as f16) {{ let z = x * x; x + x * z * eval_horner_f16_slice(z, &COEFFS16_SEG{idx}) }}"
        )?;
        }
        writeln!(w, "    else {{ let z = x * x; x + x * z * eval_horner_f16_slice(z, &COEFFS16_SEG{last}) }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;

        writeln!(w, "#[inline(always)] fn k_tan_f128(x: f128) -> f128 {{")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(
                w,
                "    {kw} x <= {upper}_SEG_BOUNDS_F128[{idx}] {{ let z = x * x; x + x * z * eval_horner_f128_slice(z, &COEFFS128_SEG{idx}) }}"
            )?;
        }
        writeln!(w, "    else {{ let z = x * x; x + x * z * eval_horner_f128_slice(z, &COEFFS128_SEG{last}) }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;
    }

    writeln!(
        w,
        "#[inline(always)] fn reconstruct_tan_f64(q: u8, t: f64) -> f64 {{"
    )?;
    writeln!(w, "    match q & 1 {{")?;
    writeln!(w, "        0 => t,")?;
    writeln!(w, "        1 => -1.0_f64 / t,")?;
    writeln!(w, "        _ => unreachable!(),")?;
    writeln!(w, "    }}")?;
    writeln!(w, "}}")?;
    writeln!(w)?;

    writeln!(
        w,
        "#[inline(always)] fn reconstruct_tan_f32(q: u8, t: f32) -> f32 {{"
    )?;
    writeln!(w, "    match q & 1 {{")?;
    writeln!(w, "        0 => t,")?;
    writeln!(w, "        1 => -1.0_f32 / t,")?;
    writeln!(w, "        _ => unreachable!(),")?;
    writeln!(w, "    }}")?;
    writeln!(w, "}}")?;
    writeln!(w)?;

    writeln!(
        w,
        "#[inline(always)] fn reconstruct_tan_f16(q: u8, t: f16) -> f16 {{"
    )?;
    writeln!(w, "    let tf = t as f32;")?;
    writeln!(w, "    let y = match q & 1 {{")?;
    writeln!(w, "        0 => tf,")?;
    writeln!(w, "        1 => -1.0_f32 / tf,")?;
    writeln!(w, "        _ => unreachable!(),")?;
    writeln!(w, "    }};")?;
    writeln!(w, "    y as f16")?;
    writeln!(w, "}}")?;
    writeln!(w)?;

    writeln!(
        w,
        "#[inline(always)] fn reconstruct_tan_f128(q: u8, t: f128) -> f128 {{"
    )?;
    writeln!(w, "    match q & 1 {{")?;
    writeln!(w, "        0 => t,")?;
    writeln!(w, "        1 => -1.0_f128 / t,")?;
    writeln!(w, "        _ => unreachable!(),")?;
    writeln!(w, "    }}")?;
    writeln!(w, "}}")?;
    writeln!(w)?;

    Ok(())
}

fn write_simple_kernel_helpers<W: std::io::Write>(
    w: &mut W,
    plan: &FunctionPlan,
    upper: &str,
    is_log: bool,
    prefer_f32_for_f16: bool,
) -> std::io::Result<()> {
    if plan.segments.len() == 1 {
        if is_log {
            writeln!(w, "#[inline(always)] fn k_simple_f64(x: f64) -> f64 {{ let u = x - 1.0_f64; let s = u / (x + 1.0_f64); let z = s * s; let q = eval_horner_f64_slice(z, &COEFFS64); 2.0_f64 * s * q.mul_add(z, 1.0_f64) }}")?;
            writeln!(w, "#[inline(always)] fn k_simple_f32(x: f32) -> f32 {{ let u = x - 1.0_f32; let s = u / (x + 1.0_f32); let z = s * s; let q = eval_horner_f32_slice(z, &COEFFS32); 2.0_f32 * s * q.mul_add(z, 1.0_f32) }}")?;
        } else {
            writeln!(
                w,
                "#[inline(always)] fn k_simple_f64(x: f64) -> f64 {{ eval_horner_f64_slice(x, &COEFFS64) }}"
            )?;
            writeln!(
                w,
                "#[inline(always)] fn k_simple_f32(x: f32) -> f32 {{ eval_horner_f32_slice(x, &COEFFS32) }}"
            )?;
        }
        if prefer_f32_for_f16 {
            if is_log {
                writeln!(w, "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{ let xf = x as f32; let u = xf - 1.0_f32; let s = u / (xf + 1.0_f32); let z = s * s; let q = eval_horner_f32_slice(z, &COEFFS32); (2.0_f32 * s * q.mul_add(z, 1.0_f32)) as f16 }}")?;
            } else {
                writeln!(
                    w,
                    "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{ k_simple_f32(x as f32) as f16 }}"
                )?;
            }
        } else {
            if is_log {
                writeln!(w, "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{ let u = x - 1.0_f16; let s = u / (x + 1.0_f16); let z = s * s; let q = eval_horner_f16_slice(z, &COEFFS16); 2.0_f16 * s * (1.0_f16 + z * q) }}")?;
            } else {
                writeln!(
                    w,
                    "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{ eval_horner_f16_slice(x, &COEFFS16) }}"
                )?;
            }
        }
        if is_log {
            writeln!(w, "#[inline(always)] fn k_simple_f128(x: f128) -> f128 {{ let u = x - 1.0_f128; let s = u / (x + 1.0_f128); let z = s * s; let q = eval_horner_f128_slice(z, &COEFFS128); 2.0_f128 * s * (1.0_f128 + z * q) }}")?;
        } else {
            writeln!(
                w,
                "#[inline(always)] fn k_simple_f128(x: f128) -> f128 {{ eval_horner_f128_slice(x, &COEFFS128) }}"
            )?;
        }
        writeln!(w)?;
    } else {
        writeln!(w, "#[inline(always)] fn k_simple_f64(x: f64) -> f64 {{")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            if is_log {
                writeln!(
                    w,
                    "    {kw} x <= {upper}_SEG_BOUNDS_F64[{idx}] {{ let u = x - 1.0_f64; let s = u / (x + 1.0_f64); let z = s * s; let q = eval_horner_f64_slice(z, &COEFFS64_SEG{idx}); 2.0_f64 * s * q.mul_add(z, 1.0_f64) }}"
                )?;
            } else {
                writeln!(
                    w,
                    "    {kw} x <= {upper}_SEG_BOUNDS_F64[{idx}] {{ eval_horner_f64_slice(x, &COEFFS64_SEG{idx}) }}"
                )?;
            }
        }
        let last = plan.segments.len() - 1;
        if is_log {
            writeln!(w, "    else {{ let u = x - 1.0_f64; let s = u / (x + 1.0_f64); let z = s * s; let q = eval_horner_f64_slice(z, &COEFFS64_SEG{last}); 2.0_f64 * s * q.mul_add(z, 1.0_f64) }}")?;
        } else {
            writeln!(
                w,
                "    else {{ eval_horner_f64_slice(x, &COEFFS64_SEG{last}) }}"
            )?;
        }
        writeln!(w, "}}")?;

        writeln!(w, "#[inline(always)] fn k_simple_f32(x: f32) -> f32 {{")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            if is_log {
                writeln!(
                    w,
                    "    {kw} x <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let u = x - 1.0_f32; let s = u / (x + 1.0_f32); let z = s * s; let q = eval_horner_f32_slice(z, &COEFFS32_SEG{idx}); 2.0_f32 * s * q.mul_add(z, 1.0_f32) }}"
                )?;
            } else {
                writeln!(
                    w,
                    "    {kw} x <= {upper}_SEG_BOUNDS_F32[{idx}] {{ eval_horner_f32_slice(x, &COEFFS32_SEG{idx}) }}"
                )?;
            }
        }
        let last = plan.segments.len() - 1;
        if is_log {
            writeln!(w, "    else {{ let u = x - 1.0_f32; let s = u / (x + 1.0_f32); let z = s * s; let q = eval_horner_f32_slice(z, &COEFFS32_SEG{last}); 2.0_f32 * s * q.mul_add(z, 1.0_f32) }}")?;
        } else {
            writeln!(
                w,
                "    else {{ eval_horner_f32_slice(x, &COEFFS32_SEG{last}) }}"
            )?;
        }
        writeln!(w, "}}")?;

        if prefer_f32_for_f16 {
            if is_log {
                writeln!(w, "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{ k_simple_f32(x as f32) as f16 }}")?;
            } else {
                writeln!(w, "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{ k_simple_f32(x as f32) as f16 }}")?;
            }
        } else {
            writeln!(w, "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                if is_log {
                    writeln!(
                        w,
                        "    {kw} x <= ({upper}_SEG_BOUNDS_F32[{idx}] as f16) {{ let u = x - 1.0_f16; let s = u / (x + 1.0_f16); let z = s * s; let q = eval_horner_f16_slice(z, &COEFFS16_SEG{idx}); 2.0_f16 * s * (1.0_f16 + z * q) }}"
                    )?;
                } else {
                    writeln!(
                        w,
                        "    {kw} x <= ({upper}_SEG_BOUNDS_F32[{idx}] as f16) {{ eval_horner_f16_slice(x, &COEFFS16_SEG{idx}) }}"
                    )?;
                }
            }
            let last = plan.segments.len() - 1;
            if is_log {
                writeln!(w, "    else {{ let u = x - 1.0_f16; let s = u / (x + 1.0_f16); let z = s * s; let q = eval_horner_f16_slice(z, &COEFFS16_SEG{last}); 2.0_f16 * s * (1.0_f16 + z * q) }}")?;
            } else {
                writeln!(
                    w,
                    "    else {{ eval_horner_f16_slice(x, &COEFFS16_SEG{last}) }}"
                )?;
            }
            writeln!(w, "}}")?;
        }

        writeln!(w, "#[inline(always)] fn k_simple_f128(x: f128) -> f128 {{")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            if is_log {
                writeln!(
                    w,
                    "    {kw} x <= {upper}_SEG_BOUNDS_F128[{idx}] {{ let u = x - 1.0_f128; let s = u / (x + 1.0_f128); let z = s * s; let q = eval_horner_f128_slice(z, &COEFFS128_SEG{idx}); 2.0_f128 * s * (1.0_f128 + z * q) }}"
                )?;
            } else {
                writeln!(
                    w,
                    "    {kw} x <= {upper}_SEG_BOUNDS_F128[{idx}] {{ eval_horner_f128_slice(x, &COEFFS128_SEG{idx}) }}"
                )?;
            }
        }
        let last = plan.segments.len() - 1;
        if is_log {
            writeln!(w, "    else {{ let u = x - 1.0_f128; let s = u / (x + 1.0_f128); let z = s * s; let q = eval_horner_f128_slice(z, &COEFFS128_SEG{last}); 2.0_f128 * s * (1.0_f128 + z * q) }}")?;
        } else {
            writeln!(
                w,
                "    else {{ eval_horner_f128_slice(x, &COEFFS128_SEG{last}) }}"
            )?;
        }
        writeln!(w, "}}")?;
        writeln!(w)?;
    }

    Ok(())
}

struct MxSpec {
    ty: &'static str,
    suffix: &'static str,
    coeff_base: &'static str,
}

fn mx_specs() -> [MxSpec; 4] {
    [
        MxSpec {
            ty: "Fp4E2M1",
            suffix: "fp4_e2m1",
            coeff_base: "COEFFS_FP4_E2M1",
        },
        MxSpec {
            ty: "Fp6E3M2",
            suffix: "fp6_e3m2",
            coeff_base: "COEFFS_FP6_E3M2",
        },
        MxSpec {
            ty: "Fp8E4M3",
            suffix: "fp8_e4m3",
            coeff_base: "COEFFS_FP8_E4M3",
        },
        MxSpec {
            ty: "Fp8E5M2",
            suffix: "fp8_e5m2",
            coeff_base: "COEFFS_FP8_E5M2",
        },
    ]
}

fn write_mx_coefficients<W: Write>(w: &mut W, plan: &FunctionPlan) -> io::Result<()> {
    for spec in mx_specs() {
        if plan.segments.len() == 1 {
            let s = &plan.segments[0];
            write!(
                w,
                "#[cfg(feature = \"mx\")] pub const {}: [{}; {}] = ",
                spec.coeff_base,
                spec.ty,
                s.c32.as_slice().len()
            )?;
            match spec.suffix {
                "fp4_e2m1" => write_coeffset_fp4_inline(w, &s.c32)?,
                "fp6_e3m2" => write_coeffset_fp6_inline(w, &s.c32)?,
                "fp8_e4m3" => write_coeffset_fp8e4_inline(w, &s.c32)?,
                "fp8_e5m2" => write_coeffset_fp8e5_inline(w, &s.c32)?,
                _ => unreachable!(),
            }
            writeln!(w, ";")?;
        } else {
            for (idx, s) in plan.segments.iter().enumerate() {
                write!(
                    w,
                    "#[cfg(feature = \"mx\")] pub const {}_SEG{idx}: [{}; {}] = ",
                    spec.coeff_base,
                    spec.ty,
                    s.c32.as_slice().len()
                )?;
                match spec.suffix {
                    "fp4_e2m1" => write_coeffset_fp4_inline(w, &s.c32)?,
                    "fp6_e3m2" => write_coeffset_fp6_inline(w, &s.c32)?,
                    "fp8_e4m3" => write_coeffset_fp8e4_inline(w, &s.c32)?,
                    "fp8_e5m2" => write_coeffset_fp8e5_inline(w, &s.c32)?,
                    _ => unreachable!(),
                }
                writeln!(w, ";")?;
            }
        }
    }
    writeln!(w)?;
    Ok(())
}

fn write_mx_trig_imports<W: Write>(w: &mut W, is_sin: bool, is_cos: bool) -> io::Result<()> {
    let _ = (w, is_sin, is_cos);
    Ok(())
}

fn write_mx_horner_helpers<W: Write>(w: &mut W) -> io::Result<()> {
    writeln!(w)?;
    Ok(())
}

fn write_mx_kernel_helpers<W: Write>(
    w: &mut W,
    upper: &str,
    _func_name: &str,
    plan: &FunctionPlan,
    is_sin: bool,
    is_cos: bool,
    is_tan: bool,
) -> io::Result<()> {
    for spec in mx_specs() {
        writeln!(w, "#[cfg(feature = \"mx\")]")?;
        writeln!(
            w,
            "#[inline(always)] fn eval_horner_{}_slice(x: {}, cs: &[{}]) -> {} {{",
            spec.suffix, spec.ty, spec.ty, spec.ty
        )?;
        writeln!(w, "    let mut p = {}::from(0.0_f32);", spec.ty)?;
        writeln!(w, "    for &c in cs.iter().rev() {{ p = p * x + c; }}")?;
        writeln!(w, "    p")?;
        writeln!(w, "}}")?;
        writeln!(w)?;

        if is_sin || is_cos {
            if is_sin {
                writeln!(w, "#[cfg(feature = \"mx\")]")?;
                writeln!(
                    w,
                    "#[inline(always)] pub(crate) fn sin_kernel_{}(x: {}) -> {} {{",
                    spec.suffix, spec.ty, spec.ty
                )?;
                if plan.segments.len() == 1 {
                    writeln!(w, "    let z = x * x;")?;
                    writeln!(
                        w,
                        "    x * eval_horner_{}_slice(z, &{})",
                        spec.suffix, spec.coeff_base
                    )?;
                } else {
                    writeln!(w, "    let axf = f32::from(x).abs();")?;
                    for (idx, _) in plan.segments.iter().enumerate() {
                        let kw = if idx == 0 { "if" } else { "else if" };
                        writeln!(w, "    {kw} axf <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = x * x; x * eval_horner_{}_slice(z, &{}_SEG{idx}) }}", spec.suffix, spec.coeff_base)?;
                    }
                    let last = plan.segments.len() - 1;
                    writeln!(
                        w,
                        "    else {{ let z = x * x; x * eval_horner_{}_slice(z, &{}_SEG{last}) }}",
                        spec.suffix, spec.coeff_base
                    )?;
                }
                writeln!(w, "}}")?;
            } else {
                writeln!(w, "#[cfg(feature = \"mx\")]")?;
                writeln!(
                    w,
                    "#[inline(always)] pub(crate) fn cos_kernel_{}(x: {}) -> {} {{",
                    spec.suffix, spec.ty, spec.ty
                )?;
                if plan.segments.len() == 1 {
                    writeln!(w, "    let z = x * x;")?;
                    writeln!(
                        w,
                        "    eval_horner_{}_slice(z, &{})",
                        spec.suffix, spec.coeff_base
                    )?;
                } else {
                    writeln!(w, "    let axf = f32::from(x).abs();")?;
                    for (idx, _) in plan.segments.iter().enumerate() {
                        let kw = if idx == 0 { "if" } else { "else if" };
                        writeln!(w, "    {kw} axf <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = x * x; eval_horner_{}_slice(z, &{}_SEG{idx}) }}", spec.suffix, spec.coeff_base)?;
                    }
                    let last = plan.segments.len() - 1;
                    writeln!(
                        w,
                        "    else {{ let z = x * x; eval_horner_{}_slice(z, &{}_SEG{last}) }}",
                        spec.suffix, spec.coeff_base
                    )?;
                }
                writeln!(w, "}}")?;
            }

            writeln!(w, "#[cfg(feature = \"mx\")]")?;
            writeln!(w, "#[inline(always)] fn k_sin_{}(x: {}) -> {} {{ crate::trig::sin::sin_kernel_{}(x) }}", spec.suffix, spec.ty, spec.ty, spec.suffix)?;
            writeln!(w, "#[cfg(feature = \"mx\")]")?;
            writeln!(w, "#[inline(always)] fn k_cos_{}(x: {}) -> {} {{ crate::trig::cos::cos_kernel_{}(x) }}", spec.suffix, spec.ty, spec.ty, spec.suffix)?;
            writeln!(w)?;
            continue;
        }

        let helper_name = if is_tan { "k_tan" } else { "k_simple" };
        writeln!(w, "#[cfg(feature = \"mx\")]")?;
        writeln!(
            w,
            "#[inline(always)] fn {}_{}(x: {}) -> {} {{",
            helper_name, spec.suffix, spec.ty, spec.ty
        )?;
        if plan.segments.len() == 1 {
            if is_tan {
                writeln!(w, "    let z = x * x;")?;
                writeln!(
                    w,
                    "    x + x * z * eval_horner_{}_slice(z, &{})",
                    spec.suffix, spec.coeff_base
                )?;
            } else {
                writeln!(
                    w,
                    "    eval_horner_{}_slice(x, &{})",
                    spec.suffix, spec.coeff_base
                )?;
            }
        } else {
            writeln!(w, "    let xf = f32::from(x);")?;
            for (idx, _) in plan.segments.iter().enumerate() {
                let kw = if idx == 0 { "if" } else { "else if" };
                if is_tan {
                    writeln!(w, "    {kw} xf <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = x * x; x + x * z * eval_horner_{}_slice(z, &{}_SEG{idx}) }}", spec.suffix, spec.coeff_base)?;
                } else {
                    writeln!(w, "    {kw} xf <= {upper}_SEG_BOUNDS_F32[{idx}] {{ eval_horner_{}_slice(x, &{}_SEG{idx}) }}", spec.suffix, spec.coeff_base)?;
                }
            }
            let last = plan.segments.len() - 1;
            if is_tan {
                writeln!(w, "    else {{ let z = x * x; x + x * z * eval_horner_{}_slice(z, &{}_SEG{last}) }}", spec.suffix, spec.coeff_base)?;
            } else {
                writeln!(
                    w,
                    "    else {{ eval_horner_{}_slice(x, &{}_SEG{last}) }}",
                    spec.suffix, spec.coeff_base
                )?;
            }
        }
        writeln!(w, "}}")?;
        writeln!(w)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_mx_public_functions<W: Write>(
    w: &mut W,
    func_name: &str,
    ident: &str,
    upper: &str,
) -> io::Result<()> {
    let ctx = common::PublicFnCtx { ident, upper };

    for spec in mx_specs() {
        let fn_name = format!("{}_{}", ident, spec.suffix);
        writeln!(w)?;
        writeln!(w, "#[cfg(feature = \"mx\")]")?;
        writeln!(w, "#[inline]")?;
        writeln!(w, "pub fn {fn_name}(x: {}) -> {} {{", spec.ty, spec.ty)?;
        common::write_mx_public_body(w, func_name, &ctx, &fn_name, &spec)?;
        writeln!(w, "}}")?;
    }
    Ok(())
}

/* ------------------------------------------------------------------------------------------------
 * C codegen
 * ------------------------------------------------------------------------------------------------ */

pub fn write_coeffs_to_c<const N: usize, const K: usize, F>(
    rz: &RemezApprox<f128, N, K>,
    func_name: &str,
    f_x: F,
    interval_x: (f128, f128),
    max_error: f128,
    file_prefix: &str,
) where
    F: Fn(f128) -> f128,
{
    common::generator_for(func_name)
        .unwrap_or_else(|err| panic!("{err}"))
        .c_generate(rz, func_name, f_x, interval_x, max_error, file_prefix);
}

fn write_coeffs_to_c_shared<const N: usize, const K: usize, F>(
    rz: &RemezApprox<f128, N, K>,
    func_name: &str,
    f_x: F,
    interval_x: (f128, f128),
    max_error: f128,
    file_prefix: &str,
) where
    F: Fn(f128) -> f128,
{
    use std::io::BufWriter;
    use std::{fs, io::Write};

    let dom = poly_domain(func_name);
    let base_out = out_dir_for("c", "target/c");
    fs::create_dir_all(&base_out).unwrap();
    let path = base_out.join(format!("{file_prefix}.c"));
    let mut w = BufWriter::new(fs::File::create(&path).unwrap());

    let export = RemezExport::<N, K>::from_remez(
        rz,
        func_name.to_string(),
        f_x,
        interval_x,
        interval_x,
        max_error,
        dom,
    );

    writeln!(w, "/* Autogenerated by aether-fpx-gen (remez) for {func_name} */").unwrap();
    writeln!(
        w,
        "/* interval_f64: [{:.17e}, {:.17e}] */",
        export.interval_x.0, export.interval_x.1
    )
    .unwrap();
    writeln!(
        w,
        "/* interval_f32: [{:.9e}, {:.9e}] */",
        export.interval_x.0 as f32, export.interval_x.1 as f32
    )
    .unwrap();
    writeln!(w, "/* table_domain: {} */", export.table_domain).unwrap();
    writeln!(w, "#include <math.h>").unwrap();
    writeln!(w, "#include \"{file_prefix}.h\"").unwrap();
    writeln!(w).unwrap();

    writeln!(
        w,
        "static inline double mirror_f64(double x) {{ return x < 0.0 ? -x : x; }}"
    )
    .unwrap();
    writeln!(
        w,
        "static inline float  mirror_f32(float  x) {{ return x < 0.0f ? -x : x; }}"
    )
    .unwrap();
    writeln!(w).unwrap();
}

/* ------------------------------------------------------------------------------------------------
 * Function JSON writer
 * ------------------------------------------------------------------------------------------------ */

pub fn write_function_coeffs_to_json(export: &FunctionRemezExport) {
    use std::{fs, io::BufWriter};

    let base_out = out_dir_for("json", "target/coeffs");
    fs::create_dir_all(&base_out).unwrap();

    let path = base_out.join(format!("{}.json", export.function));
    let file = fs::File::create(&path).unwrap();
    let writer = BufWriter::new(file);

    serde_json::to_writer_pretty(writer, export).unwrap();
}

/* ------------------------------------------------------------------------------------------------
 * Coefficient writers
 * ------------------------------------------------------------------------------------------------ */

use std::io::{self, Write};

pub fn write_coeffset64_inline<W: Write>(w: &mut W, c: &CoeffSet64) -> io::Result<()> {
    match c {
        CoeffSet64::N2(a) => write_array64(w, "CoeffSet64::N2", a),
        CoeffSet64::N3(a) => write_array64(w, "CoeffSet64::N3", a),
        CoeffSet64::N4(a) => write_array64(w, "CoeffSet64::N4", a),
        CoeffSet64::N5(a) => write_array64(w, "CoeffSet64::N5", a),
        CoeffSet64::N6(a) => write_array64(w, "CoeffSet64::N6", a),
        CoeffSet64::N7(a) => write_array64(w, "CoeffSet64::N7", a),
        CoeffSet64::N8(a) => write_array64(w, "CoeffSet64::N8", a),
        CoeffSet64::N9(a) => write_array64(w, "CoeffSet64::N9", a),
        CoeffSet64::N10(a) => write_array64(w, "CoeffSet64::N10", a),
        CoeffSet64::N11(a) => write_array64(w, "CoeffSet64::N11", a),
        CoeffSet64::N12(a) => write_array64(w, "CoeffSet64::N12", a),
        CoeffSet64::N13(a) => write_array64(w, "CoeffSet64::N13", a),
        CoeffSet64::N14(a) => write_array64(w, "CoeffSet64::N14", a),
        CoeffSet64::N15(a) => write_array64(w, "CoeffSet64::N15", a),
        CoeffSet64::N16(a) => write_array64(w, "CoeffSet64::N16", a),
        CoeffSet64::N17(a) => write_array64(w, "CoeffSet64::N17", a),
        CoeffSet64::N18(a) => write_array64(w, "CoeffSet64::N18", a),
        CoeffSet64::N19(a) => write_array64(w, "CoeffSet64::N19", a),
        CoeffSet64::N20(a) => write_array64(w, "CoeffSet64::N20", a),
        CoeffSet64::N21(a) => write_array64(w, "CoeffSet64::N21", a),
        CoeffSet64::N22(a) => write_array64(w, "CoeffSet64::N22", a),
        CoeffSet64::N23(a) => write_array64(w, "CoeffSet64::N23", a),
        CoeffSet64::N24(a) => write_array64(w, "CoeffSet64::N24", a),
        CoeffSet64::N25(a) => write_array64(w, "CoeffSet64::N25", a),
    }
}

pub fn write_coeffset128_inline<W: Write>(w: &mut W, c: &CoeffSet128) -> io::Result<()> {
    match c {
        CoeffSet128::N2(a) => write_array128(w, "CoeffSet128::N2", a),
        CoeffSet128::N3(a) => write_array128(w, "CoeffSet128::N3", a),
        CoeffSet128::N4(a) => write_array128(w, "CoeffSet128::N4", a),
        CoeffSet128::N5(a) => write_array128(w, "CoeffSet128::N5", a),
        CoeffSet128::N6(a) => write_array128(w, "CoeffSet128::N6", a),
        CoeffSet128::N7(a) => write_array128(w, "CoeffSet128::N7", a),
        CoeffSet128::N8(a) => write_array128(w, "CoeffSet128::N8", a),
        CoeffSet128::N9(a) => write_array128(w, "CoeffSet128::N9", a),
        CoeffSet128::N10(a) => write_array128(w, "CoeffSet128::N10", a),
        CoeffSet128::N11(a) => write_array128(w, "CoeffSet128::N11", a),
        CoeffSet128::N12(a) => write_array128(w, "CoeffSet128::N12", a),
        CoeffSet128::N13(a) => write_array128(w, "CoeffSet128::N13", a),
        CoeffSet128::N14(a) => write_array128(w, "CoeffSet128::N14", a),
        CoeffSet128::N15(a) => write_array128(w, "CoeffSet128::N15", a),
        CoeffSet128::N16(a) => write_array128(w, "CoeffSet128::N16", a),
        CoeffSet128::N17(a) => write_array128(w, "CoeffSet128::N17", a),
        CoeffSet128::N18(a) => write_array128(w, "CoeffSet128::N18", a),
        CoeffSet128::N19(a) => write_array128(w, "CoeffSet128::N19", a),
        CoeffSet128::N20(a) => write_array128(w, "CoeffSet128::N20", a),
        CoeffSet128::N21(a) => write_array128(w, "CoeffSet128::N21", a),
        CoeffSet128::N22(a) => write_array128(w, "CoeffSet128::N22", a),
        CoeffSet128::N23(a) => write_array128(w, "CoeffSet128::N23", a),
        CoeffSet128::N24(a) => write_array128(w, "CoeffSet128::N24", a),
        CoeffSet128::N25(a) => write_array128(w, "CoeffSet128::N25", a),
    }
}

fn write_array64<W: Write, const N: usize>(w: &mut W, _tag: &str, a: &[f64; N]) -> io::Result<()> {
    write!(w, "[")?;
    for (i, v) in a.iter().enumerate() {
        if i != 0 {
            write!(w, ", ")?;
        }
        if v.is_nan() {
            write!(w, "f64::NAN")?;
        } else if v.is_infinite() {
            if v.is_sign_positive() {
                write!(w, "f64::INFINITY")?;
            } else {
                write!(w, "f64::NEG_INFINITY")?;
            }
        } else {
            write!(w, "{:.17e}_f64", v)?;
        }
    }
    write!(w, "]")
}

fn write_array128<W: Write, const N: usize>(
    w: &mut W,
    _tag: &str,
    a: &[f128; N],
) -> io::Result<()> {
    write!(w, "[")?;
    for (i, v) in a.iter().enumerate() {
        if i != 0 {
            write!(w, ", ")?;
        }
        write!(w, "f128::from_bits(0x{:032x}u128)", v.to_bits())?;
    }
    write!(w, "]")
}

pub fn write_coeffset32_inline<W: Write>(w: &mut W, c: &CoeffSet32) -> io::Result<()> {
    match c {
        CoeffSet32::N2(a) => write_array32(w, "CoeffSet32::N2", a),
        CoeffSet32::N3(a) => write_array32(w, "CoeffSet32::N3", a),
        CoeffSet32::N4(a) => write_array32(w, "CoeffSet32::N4", a),
        CoeffSet32::N5(a) => write_array32(w, "CoeffSet32::N5", a),
        CoeffSet32::N6(a) => write_array32(w, "CoeffSet32::N6", a),
        CoeffSet32::N7(a) => write_array32(w, "CoeffSet32::N7", a),
        CoeffSet32::N8(a) => write_array32(w, "CoeffSet32::N8", a),
        CoeffSet32::N9(a) => write_array32(w, "CoeffSet32::N9", a),
        CoeffSet32::N10(a) => write_array32(w, "CoeffSet32::N10", a),
        CoeffSet32::N11(a) => write_array32(w, "CoeffSet32::N11", a),
        CoeffSet32::N12(a) => write_array32(w, "CoeffSet32::N12", a),
        CoeffSet32::N13(a) => write_array32(w, "CoeffSet32::N13", a),
        CoeffSet32::N14(a) => write_array32(w, "CoeffSet32::N14", a),
        CoeffSet32::N15(a) => write_array32(w, "CoeffSet32::N15", a),
        CoeffSet32::N16(a) => write_array32(w, "CoeffSet32::N16", a),
        CoeffSet32::N17(a) => write_array32(w, "CoeffSet32::N17", a),
        CoeffSet32::N18(a) => write_array32(w, "CoeffSet32::N18", a),
        CoeffSet32::N19(a) => write_array32(w, "CoeffSet32::N19", a),
        CoeffSet32::N20(a) => write_array32(w, "CoeffSet32::N20", a),
        CoeffSet32::N21(a) => write_array32(w, "CoeffSet32::N21", a),
        CoeffSet32::N22(a) => write_array32(w, "CoeffSet32::N22", a),
        CoeffSet32::N23(a) => write_array32(w, "CoeffSet32::N23", a),
        CoeffSet32::N24(a) => write_array32(w, "CoeffSet32::N24", a),
        CoeffSet32::N25(a) => write_array32(w, "CoeffSet32::N25", a),
    }
}

fn write_array32<W: Write, const N: usize>(w: &mut W, _tag: &str, a: &[f32; N]) -> io::Result<()> {
    write!(w, "[")?;
    for (i, v) in a.iter().enumerate() {
        if i != 0 {
            write!(w, ", ")?;
        }
        if v.is_nan() {
            write!(w, "f32::NAN")?;
        } else if v.is_infinite() {
            if v.is_sign_positive() {
                write!(w, "f32::INFINITY")?;
            } else {
                write!(w, "f32::NEG_INFINITY")?;
            }
        } else {
            write!(w, "{:.9e}_f32", v)?;
        }
    }
    write!(w, "]")
}

pub fn write_coeffset16_inline<W: Write>(w: &mut W, c: &CoeffSet16) -> io::Result<()> {
    match c {
        CoeffSet16::N2(a) => write_array16(w, "CoeffSet16::N2", a),
        CoeffSet16::N3(a) => write_array16(w, "CoeffSet16::N3", a),
        CoeffSet16::N4(a) => write_array16(w, "CoeffSet16::N4", a),
        CoeffSet16::N5(a) => write_array16(w, "CoeffSet16::N5", a),
        CoeffSet16::N6(a) => write_array16(w, "CoeffSet16::N6", a),
        CoeffSet16::N7(a) => write_array16(w, "CoeffSet16::N7", a),
        CoeffSet16::N8(a) => write_array16(w, "CoeffSet16::N8", a),
        CoeffSet16::N9(a) => write_array16(w, "CoeffSet16::N9", a),
        CoeffSet16::N10(a) => write_array16(w, "CoeffSet16::N10", a),
        CoeffSet16::N11(a) => write_array16(w, "CoeffSet16::N11", a),
        CoeffSet16::N12(a) => write_array16(w, "CoeffSet16::N12", a),
        CoeffSet16::N13(a) => write_array16(w, "CoeffSet16::N13", a),
        CoeffSet16::N14(a) => write_array16(w, "CoeffSet16::N14", a),
        CoeffSet16::N15(a) => write_array16(w, "CoeffSet16::N15", a),
        CoeffSet16::N16(a) => write_array16(w, "CoeffSet16::N16", a),
        CoeffSet16::N17(a) => write_array16(w, "CoeffSet16::N17", a),
        CoeffSet16::N18(a) => write_array16(w, "CoeffSet16::N18", a),
        CoeffSet16::N19(a) => write_array16(w, "CoeffSet16::N19", a),
        CoeffSet16::N20(a) => write_array16(w, "CoeffSet16::N20", a),
        CoeffSet16::N21(a) => write_array16(w, "CoeffSet16::N21", a),
        CoeffSet16::N22(a) => write_array16(w, "CoeffSet16::N22", a),
        CoeffSet16::N23(a) => write_array16(w, "CoeffSet16::N23", a),
        CoeffSet16::N24(a) => write_array16(w, "CoeffSet16::N24", a),
        CoeffSet16::N25(a) => write_array16(w, "CoeffSet16::N25", a),
    }
}

pub fn write_coeffset_fp4_inline<W: Write>(w: &mut W, c: &CoeffSet32) -> io::Result<()> {
    write_array_fp4(w, c.as_slice())
}

pub fn write_coeffset_fp6_inline<W: Write>(w: &mut W, c: &CoeffSet32) -> io::Result<()> {
    write_array_fp6(w, c.as_slice())
}

pub fn write_coeffset_fp8e4_inline<W: Write>(w: &mut W, c: &CoeffSet32) -> io::Result<()> {
    write_array_fp8e4(w, c.as_slice())
}

pub fn write_coeffset_fp8e5_inline<W: Write>(w: &mut W, c: &CoeffSet32) -> io::Result<()> {
    write_array_fp8e5(w, c.as_slice())
}

fn write_array16<W: Write, const N: usize>(w: &mut W, _tag: &str, a: &[f16; N]) -> io::Result<()> {
    write!(w, "[")?;
    for (i, v) in a.iter().enumerate() {
        if i != 0 {
            write!(w, ", ")?;
        }
        let vf32 = (*v) as f32;
        if vf32.is_nan() {
            write!(w, "f16::NAN")?;
        } else if vf32.is_infinite() {
            if vf32.is_sign_positive() {
                write!(w, "f16::INFINITY")?;
            } else {
                write!(w, "f16::NEG_INFINITY")?;
            }
        } else {
            write!(w, "{:.5e}_f16", vf32)?;
        }
    }
    write!(w, "]")
}

fn write_array_fp4<W: Write>(w: &mut W, a: &[f32]) -> io::Result<()> {
    write!(w, "[")?;
    for (i, v) in a.iter().enumerate() {
        if i != 0 {
            write!(w, ", ")?;
        }
        write!(
            w,
            "Fp4E2M1::from_bits(0x{:02x})",
            Fp4E2M1::from(*v).to_bits()
        )?;
    }
    write!(w, "]")
}

fn write_array_fp6<W: Write>(w: &mut W, a: &[f32]) -> io::Result<()> {
    write!(w, "[")?;
    for (i, v) in a.iter().enumerate() {
        if i != 0 {
            write!(w, ", ")?;
        }
        write!(
            w,
            "Fp6E3M2::from_bits(0x{:02x})",
            Fp6E3M2::from(*v).to_bits()
        )?;
    }
    write!(w, "]")
}

fn write_array_fp8e4<W: Write>(w: &mut W, a: &[f32]) -> io::Result<()> {
    write!(w, "[")?;
    for (i, v) in a.iter().enumerate() {
        if i != 0 {
            write!(w, ", ")?;
        }
        write!(
            w,
            "Fp8E4M3::from_bits(0x{:02x})",
            Fp8E4M3::from(*v).to_bits()
        )?;
    }
    write!(w, "]")
}

fn write_array_fp8e5<W: Write>(w: &mut W, a: &[f32]) -> io::Result<()> {
    write!(w, "[")?;
    for (i, v) in a.iter().enumerate() {
        if i != 0 {
            write!(w, ", ")?;
        }
        write!(
            w,
            "Fp8E5M2::from_bits(0x{:02x})",
            Fp8E5M2::from(*v).to_bits()
        )?;
    }
    write!(w, "]")
}

#[cfg(test)]
mod gpu_export_tests {
    use super::{
        build_gpu_backend_export, gpu_eval_source, gpu_has_full_domain_implementation,
        render_gpu_backend_file, render_gpu_lowp_dispatch, render_gpu_public_dispatch,
        FunctionRemezExport, FunctionRemezSegmentExport,
    };

    const GPU_FUNCTIONS: [&str; 13] = [
        "sin", "cos", "tan", "sinh", "cosh", "tanh", "asin", "acos", "atan", "log", "log1p", "exp",
        "expm1",
    ];
    const GPU_FUNCTIONS_WITHOUT_FULL_DOMAIN_DISPATCH: [&str; 12] = [
        "sin", "cos", "tan", "sinh", "cosh", "tanh", "asin", "acos", "atan", "log1p", "exp",
        "expm1",
    ];
    const GPU_PRECISIONS: [&str; 7] = [
        "double", "float", "f16", "fp4_e2m1", "fp6_e3m2", "fp8_e4m3", "fp8_e5m2",
    ];

    fn test_export(function: &str) -> FunctionRemezExport {
        let table_domain = if matches!(function, "sin" | "cos" | "tan") {
            "z"
        } else {
            "x"
        };
        let mut export = FunctionRemezExport::new(function.into(), table_domain.into());
        for (idx, interval_x) in [(-0.25, 0.125), (0.125, 0.75)].into_iter().enumerate() {
            let coeffs_f64 = vec![1.0, 0.5];
            let coeffs_f32 = vec![1.0, 0.5];
            let segment_tag = format!("seg{idx}");
            export.segments.push(FunctionRemezSegmentExport {
                interval_x,
                interval_solve: interval_x,
                max_error_f128_abs: 0.0,
                coeffs_f64: coeffs_f64.clone(),
                coeffs_f32: coeffs_f32.clone(),
                coeffs_f16_bits: vec![0x3c00, 0x3800],
                coeffs_fp4_e2m1_bits: vec![0x02, 0x01],
                coeffs_fp6_e3m2_bits: vec![0x0c, 0x08],
                coeffs_fp8_e4m3_bits: vec![0x38, 0x30],
                coeffs_fp8_e5m2_bits: vec![0x3c, 0x38],
                cuda: build_gpu_backend_export(
                    "cuda",
                    function,
                    &segment_tag,
                    table_domain,
                    &coeffs_f64,
                    &coeffs_f32,
                ),
                rocm: build_gpu_backend_export(
                    "rocm",
                    function,
                    &segment_tag,
                    table_domain,
                    &coeffs_f64,
                    &coeffs_f32,
                ),
            });
        }
        export
    }

    fn global_kernel_body<'a>(source: &'a str, kernel_name: &str) -> &'a str {
        let signature = format!("extern \"C\" __global__ void {kernel_name}(");
        device_function_body(source, &signature)
    }

    fn device_function_body<'a>(source: &'a str, signature: &str) -> &'a str {
        let start = source
            .find(&signature)
            .unwrap_or_else(|| panic!("missing {signature} in generated source"));
        let rest = &source[start..];
        let end = rest
            .find("\n}")
            .unwrap_or_else(|| panic!("unterminated {signature} in generated source"));
        &rest[..end + 2]
    }

    fn lowp_types(backend: &str, suffix: &str) -> (&'static str, &'static str) {
        let alias = match suffix {
            "f16" => "fpx_half",
            "fp4_e2m1" => "fpx_fp4_e2m1",
            "fp6_e3m2" => "fpx_fp6_e3m2",
            "fp8_e4m3" => "fpx_fp8_e4m3",
            "fp8_e5m2" => "fpx_fp8_e5m2",
            _ => panic!("not a low-precision suffix: {suffix}"),
        };
        let backend_ty = match (backend, suffix) {
            ("cuda", "f16") | ("rocm", "f16") => "__half",
            ("cuda", "fp4_e2m1") => "__nv_fp4_e2m1",
            ("cuda", "fp6_e3m2") => "__nv_fp6_e3m2",
            ("cuda", "fp8_e4m3") => "__nv_fp8_e4m3",
            ("cuda", "fp8_e5m2") => "__nv_fp8_e5m2",
            ("rocm", "fp4_e2m1") => "__hip_fp4_e2m1",
            ("rocm", "fp6_e3m2") => "__hip_fp6_e3m2",
            ("rocm", "fp8_e4m3") => "__hip_fp8_e4m3",
            ("rocm", "fp8_e5m2") => "__hip_fp8_e5m2",
            _ => panic!("unsupported backend/type pair: {backend} {suffix}"),
        };
        (alias, backend_ty)
    }

    #[test]
    fn only_log_is_advertised_as_a_full_domain_gpu_implementation() {
        assert!(gpu_has_full_domain_implementation("log"));
        for function in GPU_FUNCTIONS_WITHOUT_FULL_DOMAIN_DISPATCH {
            assert!(
                !gpu_has_full_domain_implementation(function),
                "{function} must remain disabled until its generated approximation has full-domain dispatch"
            );
        }
    }

    #[test]
    fn every_gpu_function_exports_explicit_reduced_domain_kernels_and_accessors() {
        for backend in ["cuda", "rocm"] {
            for function in GPU_FUNCTIONS {
                let export = test_export(function);
                let source = render_gpu_backend_file(&export, backend);

                assert!(
                    source.contains(&format!(
                        "extern \"C\" double {function}_reduced_domain_lo() {{ return -2.50000000000000000e-1; }}"
                    )),
                    "missing first-segment lower bound for {backend} {function}:\n{source}"
                );
                assert!(
                    source.contains(&format!(
                        "extern \"C\" double {function}_reduced_domain_hi() {{ return 7.50000000000000000e-1; }}"
                    )),
                    "missing last-segment upper bound for {backend} {function}:\n{source}"
                );

                for precision in GPU_PRECISIONS {
                    let kernel_name = format!("{function}_reduced_kernel_{precision}");
                    let body = global_kernel_body(&source, &kernel_name);
                    let eval_name = format!("{function}_kernel_eval_{precision}");
                    assert!(
                        body.contains(&format!("ys[i] = {eval_name}(xs[i]);")),
                        "{backend} {kernel_name} does not call its generated raw evaluator:\n{body}"
                    );

                    if !matches!(precision, "double" | "float") {
                        let (alias, backend_ty) = lowp_types(backend, precision);
                        assert!(
                            source.contains(&format!(
                                "__device__ __forceinline__ {alias} {eval_name}({alias} x)"
                            )),
                            "{backend} {eval_name} does not use the tested native scalar type"
                        );
                        assert!(
                            body.contains(&format!(
                                "void {kernel_name}(const {backend_ty}* xs, {backend_ty}* ys"
                            )),
                            "{backend} {kernel_name} does not expose native buffers:\n{body}"
                        );
                        assert!(!body.contains("_to_f32(xs[i])"), "{body}");
                        assert!(!body.contains("_from_f32("), "{body}");
                    }
                }
            }
        }
    }

    #[test]
    fn full_domain_kernel_symbols_are_emitted_only_for_implemented_functions() {
        for backend in ["cuda", "rocm"] {
            for function in GPU_FUNCTIONS {
                let source = render_gpu_backend_file(&test_export(function), backend);
                for precision in GPU_PRECISIONS {
                    let full_signature =
                        format!("extern \"C\" __global__ void {function}_kernel_{precision}(");
                    assert_eq!(
                        source.contains(&full_signature),
                        function == "log",
                        "incorrect full-domain symbol availability for {backend} {function} {precision}"
                    );
                    if function == "log" {
                        let body =
                            global_kernel_body(&source, &format!("{function}_kernel_{precision}"));
                        assert!(
                            body.contains(&format!("ys[i] = {function}_eval_{precision}(xs[i]);")),
                            "{backend} full-domain {precision} call is not type-direct:\n{body}"
                        );
                        if !matches!(precision, "double" | "float") {
                            let (alias, backend_ty) = lowp_types(backend, precision);
                            assert!(
                                source.contains(&format!(
                                    "__device__ __forceinline__ {alias} log_eval_{precision}({alias} x)"
                                )),
                                "{backend} full-domain log {precision} evaluator is not native"
                            );
                            assert!(
                                body.contains(&format!(
                                    "void log_kernel_{precision}(const {backend_ty}* xs, {backend_ty}* ys"
                                )),
                                "{backend} full-domain log {precision} buffers are not native:\n{body}"
                            );
                        }
                    }
                }
                assert!(!source.contains(&format!("return {function}(x);")));
                assert!(!source.contains(&format!("return {function}f(x);")));
            }
        }
    }

    #[test]
    fn public_gpu_log_uses_range_reduced_generated_kernels() {
        let export = FunctionRemezExport::new("log".into(), "x".into());

        for scalar_ty in ["double", "float"] {
            let source = render_gpu_public_dispatch(&export, scalar_ty);
            assert!(
                source.contains("frexp"),
                "missing range reduction:\n{source}"
            );
            assert!(
                source.contains("log_kernel_eval_"),
                "missing generated kernel call:\n{source}"
            );
            assert!(!source.contains("return log("));
            assert!(!source.contains("return logf("));
        }

        for suffix in ["f16", "fp4_e2m1", "fp6_e3m2", "fp8_e4m3", "fp8_e5m2"] {
            let source = render_gpu_lowp_dispatch(&export, suffix);
            let (alias, _) = lowp_types("cuda", suffix);
            assert!(
                source.contains("frexpf"),
                "missing range reduction:\n{source}"
            );
            assert!(
                source.contains(&format!("{alias} m = fpx_f32_to_")),
                "missing native range-reduced value:\n{source}"
            );
            assert!(
                source.contains(&format!("base = log_kernel_eval_{suffix}(m);")),
                "missing generated low-precision kernel call:\n{source}"
            );
            assert!(!source.contains("return logf("));
        }
    }

    #[test]
    fn generated_f16_horner_uses_native_half_arithmetic() {
        for backend in ["cuda", "rocm"] {
            let source = render_gpu_backend_file(&test_export("sin"), backend);
            let signature = "__device__ __forceinline__ fpx_half sin_seg0_f16(fpx_half x)";
            let body = device_function_body(&source, signature);

            assert!(body.contains("fpx_mul_f16(x, x)"), "{backend}:\n{body}");
            assert!(body.contains("fpx_mul_add_f16("), "{backend}:\n{body}");
            assert!(!body.contains("float"), "{backend}:\n{body}");
            assert!(!body.contains("_to_f32"), "{backend}:\n{body}");
            assert!(!body.contains("_from_f32"), "{backend}:\n{body}");

            for helper in [
                "fpx_half fpx_add_f16(fpx_half a, fpx_half b) { return __hadd(a, b); }",
                "fpx_half fpx_mul_f16(fpx_half a, fpx_half b) { return __hmul(a, b); }",
                "fpx_half fpx_fma_f16(fpx_half a, fpx_half b, fpx_half c) { return __hfma(a, b, c); }",
            ] {
                assert!(source.contains(helper), "{backend} missing native helper: {helper}");
            }

            let dispatch = device_function_body(
                &source,
                "__device__ __forceinline__ fpx_half sin_kernel_eval_f16(fpx_half x)",
            );
            assert!(
                dispatch.contains("fpx_abs_f16(x)"),
                "{backend}:\n{dispatch}"
            );
            assert!(dispatch.contains("fpx_le_f16("), "{backend}:\n{dispatch}");
            assert!(!dispatch.contains("float"), "{backend}:\n{dispatch}");
            assert!(!dispatch.contains("_to_f32"), "{backend}:\n{dispatch}");
        }
    }

    #[test]
    fn generated_f16_horner_is_native_for_x_domain_and_log_paths() {
        for backend in ["cuda", "rocm"] {
            for (function, expected_step) in [("exp", "fpx_mul_add_f16("), ("log", "fpx_fma_f16(")]
            {
                let source = render_gpu_backend_file(&test_export(function), backend);
                let signature =
                    format!("__device__ __forceinline__ fpx_half {function}_seg0_f16(fpx_half x)");
                let body = device_function_body(&source, &signature);
                assert!(
                    body.contains(expected_step),
                    "{backend} {function}:\n{body}"
                );
                assert!(!body.contains("float"), "{backend} {function}:\n{body}");
                assert!(!body.contains("_to_f32"), "{backend} {function}:\n{body}");
            }
        }
    }

    #[test]
    fn sub_f16_formats_use_native_boundaries_and_explicit_scalar_emulation() {
        for backend in ["cuda", "rocm"] {
            let source = render_gpu_backend_file(&test_export("sin"), backend);
            assert!(source.contains("scalar wrappers provide storage and conversion, not scalar"));
            assert!(!source.contains("fpx_quantize_"), "{backend}");
            assert!(!source.contains("fpx_encode_f16_bits"), "{backend}");

            for suffix in ["fp4_e2m1", "fp6_e3m2", "fp8_e4m3", "fp8_e5m2"] {
                let (alias, _) = lowp_types(backend, suffix);
                let signature =
                    format!("__device__ __forceinline__ {alias} sin_seg0_{suffix}({alias} x)");
                let body = device_function_body(&source, &signature);
                assert!(
                    body.contains(&format!("fpx_mul_add_{suffix}(")),
                    "{backend} {suffix}:\n{body}"
                );
                assert!(!body.contains("float"), "{backend} {suffix}:\n{body}");
            }
        }
    }

    #[test]
    fn gpu_log_segment_horner_is_explicitly_fused() {
        let f64_source = gpu_eval_source("cuda", "log", "seg0", "x", "double", "coeffs", 4);
        let f32_source = gpu_eval_source("rocm", "log", "seg0", "x", "float", "coeffs", 4);

        assert!(f64_source.contains("fma(p, z"), "{f64_source}");
        assert!(f32_source.contains("fmaf(p, z"), "{f32_source}");
        assert!(!f64_source.contains("fpx_mul_add"));
        assert!(!f32_source.contains("fpx_mul_add"));
    }
}
