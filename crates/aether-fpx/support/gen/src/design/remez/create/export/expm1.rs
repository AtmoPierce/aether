use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(Expm1Generator, "expm1");

pub(super) fn write_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    super::common::write_simple_kernel_helpers(w, plan, upper, false, false)
}

pub(super) fn write_helpers<W: Write>(w: &mut W) -> io::Result<()> {
    macro_rules! line {
		($($arg:tt)*) => {
			writeln!(w, $($arg)*)?;
		};
	}

    line!(
        "#[inline(always)] fn scale_pow2_f64(k: i32) -> f64 {{ if (-1022..=1023).contains(&k) {{ f64::from_bits(((k + 1023) as u64) << 52) }} else if (-1074..=-1023).contains(&k) {{ f64::from_bits(1_u64 << (k + 1074)) }} else if k > 1023 {{ f64::INFINITY }} else {{ 0.0_f64 }} }}"
    );
    line!(
        "#[inline(always)] fn scale_pow2_f32(k: i32) -> f32 {{ if (-126..=127).contains(&k) {{ f32::from_bits(((k + 127) as u32) << 23) }} else if (-149..=-127).contains(&k) {{ f32::from_bits(1_u32 << (k + 149)) }} else if k > 127 {{ f32::INFINITY }} else {{ 0.0_f32 }} }}"
    );
    line!(
        "#[inline(always)] fn scale_pow2_f128(k: i32) -> f128 {{ if (-16382..=16383).contains(&k) {{ f128::from_bits(((k + 16383) as u128) << 112) }} else if (-16494..=-16383).contains(&k) {{ f128::from_bits(1_u128 << (k + 16494)) }} else if k > 16383 {{ f128::INFINITY }} else {{ 0.0_f128 }} }}"
    );
    line!(
        "#[inline(always)] fn expm1_small_f64(x: f64) -> f64 {{ let mut sum = x; let mut term = x; let mut n = 2_u32; while n <= 32 {{ term *= x / n as f64; sum += term; n += 1; }} sum }}"
    );
    line!(
        "#[inline(always)] fn expm1_small_f32(x: f32) -> f32 {{ x * (1.0_f32 + x * (0.5_f32 + x * (1.6666667e-1_f32 + x * (4.1666668e-2_f32 + x * (8.3333338e-3_f32 + x * 1.3888889e-3_f32))))) }}"
    );
    line!(
        "#[inline(always)] fn expm1_small_f128(x: f128) -> f128 {{ let mut sum = x; let mut term = x; let mut n = 2_u32; while n <= 64 {{ term *= x / n as f128; sum += term; n += 1; }} sum }}"
    );
    line!(
        "#[inline(always)] fn expm1_reduced_f64(x: f64) -> f64 {{ crate::trig::exp::exp(x) - 1.0_f64 }}"
    );
    line!(
        "#[inline(always)] fn expm1_reduced_f32(x: f32) -> f32 {{ let k = (x * LOG2_E_F32).round() as i32; let kd = k as f64; let r = ((x as f64 - kd * LN_2_HI_F64) - kd * LN_2_LO_F64) as f32; let s = scale_pow2_f32(k); s * k_simple_f32(r) + (s - 1.0_f32) }}"
    );
    line!(
        "#[inline(always)] fn expm1_reduced_f128(x: f128) -> f128 {{ let k = (x * LOG2_E_F128).round() as i32; let kd = k as f128; let r = ((x - kd * LN_2_HI_F128) - kd * LN_2_MID_F128) - kd * LN_2_LO_F128; let s = scale_pow2_f128(k); s * expm1_small_f128(r) + (s - 1.0_f128) }}"
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
            line!("    if x == f64::INFINITY {{ return f64::INFINITY; }}");
            line!("    if x == f64::NEG_INFINITY {{ return -1.0_f64; }}");
            line!("    if x == 0.0_f64 {{ return x; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= f64::EPSILON {{ return x; }}");
            line!("    if ax <= 0.5_f64 {{ return expm1_small_f64(x); }}");
            line!("    expm1_reduced_f64(x)");
        }
        ScalarKind::F128 => {
            line!("    if x.is_nan() {{ return f128::NAN; }}");
            line!("    if x == f128::INFINITY {{ return f128::INFINITY; }}");
            line!("    if x == f128::NEG_INFINITY {{ return -1.0_f128; }}");
            line!("    if x == 0.0_f128 {{ return x; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= f128::EPSILON {{ return x; }}");
            line!("    if ax <= 0.5_f128 {{ return expm1_small_f128(x); }}");
            line!("    expm1_reduced_f128(x)");
        }
        ScalarKind::F32 => {
            line!("    if x.is_nan() {{ return f32::NAN; }}");
            line!("    if x == f32::INFINITY {{ return f32::INFINITY; }}");
            line!("    if x == f32::NEG_INFINITY {{ return -1.0_f32; }}");
            line!("    if x == 0.0_f32 {{ return x; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= f32::EPSILON {{ return x; }}");
            line!("    if ax <= 0.03125_f32 {{ return expm1_small_f32(x); }}");
            line!("    if ax <= {upper}_INTERVAL_MAX_F32 {{ return k_simple_f32(x); }}");
            line!("    expm1_reduced_f32(x)");
        }
        ScalarKind::F16 => {
            line!("    let x32 = x as f32;");
            line!("    if x32.is_nan() {{ return f16::NAN; }}");
            line!("    if x32 == f32::INFINITY {{ return f16::INFINITY; }}");
            line!("    if x32 == f32::NEG_INFINITY {{ return -1.0_f16; }}");
            line!("    expm1f(x32) as f16");
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
        "    if x32 == f32::INFINITY {{ return {}::from(f32::INFINITY); }}",
        spec.ty
    );
    line!(
        "    if x32 == f32::NEG_INFINITY {{ return {}::from(-1.0_f32); }}",
        spec.ty
    );
    line!("    if x32 == 0.0_f32 {{ return x; }}");
    line!("    if x32.abs() <= f32::EPSILON {{ return x; }}");
    line!(
        "    if x32 >= {upper}_INTERVAL_MIN_F32 && x32 <= {upper}_INTERVAL_MAX_F32 {{ return k_simple_{}(x); }}",
        spec.suffix
    );
    line!("    {}::from(expm1_reduced_f32(x32))", spec.ty);

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
