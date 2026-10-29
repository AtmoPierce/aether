use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(Log1pGenerator, "log1p");

pub(super) fn write_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    super::common::write_simple_kernel_helpers(w, plan, upper, false, true)
}

pub(super) fn write_helpers<W: Write>(w: &mut W) -> io::Result<()> {
    macro_rules! line {
		($($arg:tt)*) => {
			writeln!(w, $($arg)*)?;
		};
	}

    line!("#[inline(always)] fn reduce_log_f64(x: f64) -> (f64, i32) {{");
    line!("    let bits = x.to_bits();");
    line!("    let exp_bits = ((bits >> 52) & 0x7ff) as i32;");
    line!(
        "    if exp_bits == 0 {{ let (m, e) = reduce_log_f64(x * 4.503599627370496e15_f64); return (m, e - 52); }}"
    );
    line!("    let mut e = exp_bits - 1023;");
    line!("    let mant_bits = (bits & ((1_u64 << 52) - 1)) | (1023_u64 << 52);");
    line!("    let mut m = f64::from_bits(mant_bits);");
    line!("    if m < SQRT_HALF_F64 {{ m *= 2.0_f64; e -= 1; }}");
    line!("    if m >= SQRT_2_F64 {{ m *= 0.5_f64; e += 1; }}");
    line!("    (m, e)");
    line!("}}");
    line!("#[inline(always)] fn reduce_log_f32(x: f32) -> (f32, i32) {{");
    line!("    let bits = x.to_bits();");
    line!("    let exp_bits = ((bits >> 23) & 0xff) as i32;");
    line!(
        "    if exp_bits == 0 {{ let (m, e) = reduce_log_f32(x * 8.388608e6_f32); return (m, e - 23); }}"
    );
    line!("    let mut e = exp_bits - 127;");
    line!("    let mant_bits = (bits & ((1_u32 << 23) - 1)) | (127_u32 << 23);");
    line!("    let mut m = f32::from_bits(mant_bits);");
    line!("    if m < SQRT_HALF_F32 {{ m *= 2.0_f32; e -= 1; }}");
    line!("    if m >= SQRT_2_F32 {{ m *= 0.5_f32; e += 1; }}");
    line!("    (m, e)");
    line!("}}");
    line!("#[inline(always)] fn reduce_log_f128(x: f128) -> (f128, i32) {{");
    line!("    let bits = x.to_bits();");
    line!("    let exp_bits = ((bits >> 112) & 0x7fff_u128) as i32;");
    line!(
        "    if exp_bits == 0 {{ let (m, e) = reduce_log_f128(x * 5.192296858534827628530496329220096e33_f128); return (m, e - 112); }}"
    );
    line!("    let mut e = exp_bits - 16383;");
    line!("    let mant_bits = (bits & ((1_u128 << 112) - 1)) | ((16383_u128) << 112);");
    line!("    let mut m = f128::from_bits(mant_bits);");
    line!("    if m < SQRT_HALF_F128 {{ m *= 2.0_f128; e -= 1; }}");
    line!("    if m >= SQRT_2_F128 {{ m *= 0.5_f128; e += 1; }}");
    line!("    (m, e)");
    line!("}}");
    line!(
        "#[inline(always)] fn log_reduced_f64(x: f64) -> f64 {{ let (m, e) = reduce_log_f64(x); let u = m - 1.0_f64; let base = if u >= crate::trig::log1p::LOG1P_INTERVAL_MIN_F64 && u <= crate::trig::log1p::LOG1P_INTERVAL_MAX_F64 {{ crate::trig::log1p::log1p(u) }} else {{ k_simple_f64(m) }}; let ef = e as f64; base + ef * LN_2_HI_F64 + ef * LN_2_LO_F64 }}"
    );
    line!(
        "#[inline(always)] fn log_reduced_f32(x: f32) -> f32 {{ let (m, e) = reduce_log_f32(x); let u = m - 1.0_f32; let base = if u >= crate::trig::log1p::LOG1P_INTERVAL_MIN_F32 && u <= crate::trig::log1p::LOG1P_INTERVAL_MAX_F32 {{ crate::trig::log1p::log1pf(u) }} else {{ k_simple_f32(m) }}; let ef = e as f64; (base as f64 + ef * LN_2_HI_F64 + ef * LN_2_LO_F64) as f32 }}"
    );
    line!(
        "#[inline(always)] fn log_reduced_f128(x: f128) -> f128 {{ let (m, e) = reduce_log_f128(x); let u = m - 1.0_f128; let base = if u >= crate::trig::log1p::LOG1P_INTERVAL_MIN_F128 && u <= crate::trig::log1p::LOG1P_INTERVAL_MAX_F128 {{ crate::trig::log1p::log1pl(u) }} else {{ k_simple_f128(m) }}; let ef = e as f128; base + ef * LN_2_HI_F128 + ef * LN_2_MID_F128 + ef * LN_2_LO_F128 }}"
    );
    line!(
        "#[inline(always)] fn log1p_large_f64(x: f64) -> f64 {{ let y = 1.0_f64 + x; let c = x - (y - 1.0_f64); crate::trig::log::log(y) + c / y }}"
    );
    line!(
        "#[inline(always)] fn log1p_large_f32(x: f32) -> f32 {{ let y = 1.0_f32 + x; let c = x - (y - 1.0_f32); crate::trig::log::logf(y) + c / y }}"
    );
    line!(
        "#[inline(always)] fn log1p_large_f128(x: f128) -> f128 {{ let y = 1.0_f128 + x; let c = x - (y - 1.0_f128); crate::trig::log::logl(y) + c / y }}"
    );
    line!(
        "#[inline(always)] fn log1p_small_f64(x: f64) -> f64 {{ let t = x / (2.0_f64 + x); let z = t * t; let mut term = t; let mut sum = t; let mut denominator = 3_u32; while denominator <= 63 {{ term *= z; sum += term / denominator as f64; denominator += 2; }} 2.0_f64 * sum }}"
    );
    line!(
        "#[inline(always)] fn log1p_small_f128(x: f128) -> f128 {{ let t = x / (2.0_f128 + x); let z = t * t; let mut term = t; let mut sum = t; let mut denominator = 3_u32; while denominator <= 127 {{ term *= z; sum += term / denominator as f128; denominator += 2; }} 2.0_f128 * sum }}"
    );
    writeln!(w)?;
    Ok(())
}

fn write_approx_scalar_body<W: Write>(
    w: &mut W,
    ctx: &PublicFnCtx<'_>,
    kind: ScalarKind,
) -> io::Result<()> {
    let upper = ctx.upper;

    macro_rules! line {
		($($arg:tt)*) => {
			writeln!(w, $($arg)*)?;
		};
	}

    match kind {
        ScalarKind::F64 => {
            line!("    if x.is_nan() {{ return f64::NAN; }}");
            line!("    if x == -1.0_f64 {{ return f64::NEG_INFINITY; }}");
            line!("    if x < -1.0_f64 {{ return f64::NAN; }}");
            line!("    if x == f64::INFINITY {{ return f64::INFINITY; }}");
            line!("    if x == 0.0_f64 {{ return x; }}");
            line!("    if x.abs() <= f64::EPSILON {{ return x; }}");
            line!("    if x.abs() <= 0.5_f64 {{ return log1p_small_f64(x); }}");
            line!("    log1p_large_f64(x)");
        }
        ScalarKind::F128 => {
            line!("    if x.is_nan() {{ return f128::NAN; }}");
            line!("    if x == -1.0_f128 {{ return f128::NEG_INFINITY; }}");
            line!("    if x < -1.0_f128 {{ return f128::NAN; }}");
            line!("    if x == f128::INFINITY {{ return f128::INFINITY; }}");
            line!("    if x == 0.0_f128 {{ return x; }}");
            line!("    if x.abs() <= f128::EPSILON {{ return x; }}");
            line!("    if x.abs() <= 0.125_f128 {{ return log1p_small_f128(x); }}");
            line!("    log1p_large_f128(x)");
        }
        ScalarKind::F32 => {
            line!("    if x.is_nan() {{ return f32::NAN; }}");
            line!("    if x == -1.0_f32 {{ return f32::NEG_INFINITY; }}");
            line!("    if x < -1.0_f32 {{ return f32::NAN; }}");
            line!("    if x == f32::INFINITY {{ return f32::INFINITY; }}");
            line!("    if x == 0.0_f32 {{ return x; }}");
            line!("    if x.abs() <= f32::EPSILON {{ return x; }}");
            line!(
                "    if x.abs() <= f32::EPSILON.sqrt() {{ let xx = x * x; return x - 0.5_f32 * xx; }}"
            );
            line!(
                "    if x >= {upper}_INTERVAL_MIN_F32 && x <= {upper}_INTERVAL_MAX_F32 {{ return k_simple_f32(x); }}"
            );
            line!("    log1p_large_f32(x)");
        }
        ScalarKind::F16 => {
            line!("    let x = x as f32;");
            line!("    if x.is_nan() {{ return f16::NAN; }}");
            line!("    if x == -1.0_f32 {{ return f16::NEG_INFINITY; }}");
            line!("    if x < -1.0_f32 {{ return f16::NAN; }}");
            line!("    if x == f32::INFINITY {{ return f16::INFINITY; }}");
            line!("    if x == 0.0_f32 {{ return x as f16; }}");
            line!("    if x.abs() <= f32::EPSILON {{ return x as f16; }}");
            line!(
                "    if x.abs() <= f32::EPSILON.sqrt() {{ let xx = x * x; return (x - 0.5_f32 * xx) as f16; }}"
            );
            line!(
                "    if x >= {upper}_INTERVAL_MIN_F32 && x <= {upper}_INTERVAL_MAX_F32 {{ return k_simple_f16(x as f16); }}"
            );
            line!("    log1p_large_f32(x) as f16");
        }
    }

    Ok(())
}

fn write_approx_mx_body<W: Write>(
    w: &mut W,
    ctx: &PublicFnCtx<'_>,
    _fn_name: &str,
    spec: &super::MxSpec,
) -> io::Result<()> {
    let upper = ctx.upper;

    macro_rules! line {
		($($arg:tt)*) => {
			writeln!(w, $($arg)*)?;
		};
	}

    line!("    let x32 = f32::from(x);");
    line!(
        "    if x32.is_nan() {{ return {}::from(f32::NAN); }}",
        spec.ty
    );
    line!(
        "    if x32 == -1.0_f32 {{ return {}::from(f32::NEG_INFINITY); }}",
        spec.ty
    );
    line!(
        "    if x32 < -1.0_f32 {{ return {}::from(f32::NAN); }}",
        spec.ty
    );
    line!(
        "    if x32 == f32::INFINITY {{ return {}::from(f32::INFINITY); }}",
        spec.ty
    );
    line!("    if x32 == 0.0_f32 {{ return x; }}");
    line!("    if x32.abs() <= f32::EPSILON {{ return x; }}");
    line!(
        "    if x32.abs() <= f32::EPSILON.sqrt() {{ let xx = x32 * x32; return {}::from(x32 - 0.5_f32 * xx); }}",
        spec.ty
    );
    line!(
        "    if x32.abs() <= {upper}_INTERVAL_MAX_F32 {{ return k_simple_{}(x); }}",
        spec.suffix
    );
    line!(
        "    crate::trig::log::log_{}(x + {}::from(1.0_f32))",
        spec.suffix,
        spec.ty
    );

    Ok(())
}

pub(super) fn write_scalar_body<W: Write>(
    w: &mut W,
    ctx: &PublicFnCtx<'_>,
    kind: ScalarKind,
) -> io::Result<()> {
    write_approx_scalar_body(w, ctx, kind)
}

pub(super) fn write_mx_body<W: Write>(
    w: &mut W,
    ctx: &PublicFnCtx<'_>,
    fn_name: &str,
    spec: &super::MxSpec,
) -> io::Result<()> {
    write_approx_mx_body(w, ctx, fn_name, spec)
}
