use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(ExpGenerator, "exp");

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

    line!("#[inline(always)] fn scale_pow2_f64(k: i32) -> f64 {{ if (-1022..=1023).contains(&k) {{ f64::from_bits(((k + 1023) as u64) << 52) }} else if (-1074..=-1023).contains(&k) {{ f64::from_bits(1_u64 << (k + 1074)) }} else if k > 1023 {{ f64::INFINITY }} else {{ 0.0_f64 }} }}");
    line!("#[inline(always)] fn scale_pow2_f32(k: i32) -> f32 {{ if (-126..=127).contains(&k) {{ f32::from_bits(((k + 127) as u32) << 23) }} else if (-149..=-127).contains(&k) {{ f32::from_bits(1_u32 << (k + 149)) }} else if k > 127 {{ f32::INFINITY }} else {{ 0.0_f32 }} }}");
    line!("#[inline(always)] fn scale_pow2_f128(k: i32) -> f128 {{ if (-16382..=16383).contains(&k) {{ f128::from_bits(((k + 16383) as u128) << 112) }} else if (-16494..=-16383).contains(&k) {{ f128::from_bits(1_u128 << (k + 16494)) }} else if k > 16383 {{ f128::INFINITY }} else {{ 0.0_f128 }} }}");
    line!("#[inline(always)] fn exp_series_f128(x: f128) -> f128 {{ let mut sum = 1.0_f128; let mut term = 1.0_f128; let mut n = 1_u32; while n <= 48 {{ term *= x / n as f128; sum += term; n += 1; }} sum }}");
    line!("#[inline(always)] fn exp_reduced_f64(x: f64) -> f64 {{ let k = (x * LOG2_E_F64).round() as i32; let kd = k as f64; let r = (x - kd * LN_2_HI_F64) - kd * LN_2_LO_F64; k_simple_f64(r) * scale_pow2_f64(k) }}");
    line!("#[inline(always)] fn exp_reduced_f32(x: f32) -> f32 {{ let k = ((x as f64) * LOG2_E_F64).round() as i32; let kd = k as f64; let r = ((x as f64 - kd * LN_2_HI_F64) - kd * LN_2_LO_F64) as f32; k_simple_f32(r) * scale_pow2_f32(k) }}");
    line!("#[inline(always)] fn exp_reduced_f128(x: f128) -> f128 {{ let k = (x * LOG2_E_F128).round() as i32; let kd = k as f128; let r = ((x - kd * LN_2_HI_F128) - kd * LN_2_MID_F128) - kd * LN_2_LO_F128; let y = exp_series_f128(r); if k == 16384 {{ (y * scale_pow2_f128(16383)) * 2.0_f128 }} else {{ y * scale_pow2_f128(k) }} }}");
    writeln!(w)?;
    Ok(())
}

fn write_approx_scalar_body<W: Write>(
    w: &mut W,
    _ctx: &PublicFnCtx<'_>,
    kind: ScalarKind,
) -> io::Result<()> {
    macro_rules! line {
		($($arg:tt)*) => {
			writeln!(w, $($arg)*)?;
		};
	}

    match kind {
        ScalarKind::F64 => {
            line!("    if x.is_nan() {{ return f64::NAN; }}");
            line!("    if x == f64::INFINITY {{ return f64::INFINITY; }}");
            line!("    if x == f64::NEG_INFINITY {{ return 0.0_f64; }}");
            line!("    if x > EXP_MAX_F64 {{ return f64::INFINITY; }}");
            line!("    if x < EXP_MIN_F64 {{ return 0.0_f64; }}");
            line!("    exp_reduced_f64(x)");
        }
        ScalarKind::F128 => {
            line!("    if x.is_nan() {{ return f128::NAN; }}");
            line!("    if x == f128::INFINITY {{ return f128::INFINITY; }}");
            line!("    if x == f128::NEG_INFINITY {{ return 0.0_f128; }}");
            line!("    if x > EXP_MAX_F128 {{ return f128::INFINITY; }}");
            line!("    if x < EXP_MIN_F128 {{ return 0.0_f128; }}");
            line!("    exp_reduced_f128(x)");
        }
        ScalarKind::F32 => {
            line!("    if x.is_nan() {{ return f32::NAN; }}");
            line!("    if x == f32::INFINITY {{ return f32::INFINITY; }}");
            line!("    if x == f32::NEG_INFINITY {{ return 0.0_f32; }}");
            line!("    if x > EXP_MAX_F32 {{ return f32::INFINITY; }}");
            line!("    if x < EXP_MIN_F32 {{ return 0.0_f32; }}");
            line!("    exp_reduced_f32(x)");
        }
        ScalarKind::F16 => {
            line!("    let x32 = x as f32;");
            line!("    if x32.is_nan() {{ return f16::NAN; }}");
            line!("    if x32 == f32::INFINITY {{ return f16::INFINITY; }}");
            line!("    if x32 == f32::NEG_INFINITY {{ return 0.0_f16; }}");
            line!("    let k = (x32 * LOG2_E_F32).round() as i32;");
            line!("    let r = (x32 - (k as f32) * LN_2_F32) as f16;");
            line!("    k_simple_f16(r) * (scale_pow2_f32(k) as f16)");
        }
    }

    Ok(())
}

fn write_approx_mx_body<W: Write>(
    w: &mut W,
    _ctx: &PublicFnCtx<'_>,
    _fn_name: &str,
    spec: &super::MxSpec,
) -> io::Result<()> {
    macro_rules! line {
		($($arg:tt)*) => {
			writeln!(w, $($arg)*)?;
		};
	}

    line!("    let x = f32::from(x);");
    line!(
        "    if !x.is_finite() {{ return {}::from(f32::NAN); }}",
        spec.ty
    );
    line!(
        "    if x > EXP_MAX_F32 {{ return {}::from(f32::INFINITY); }}",
        spec.ty
    );
    line!(
        "    if x < EXP_MIN_F32 {{ return {}::from(0.0_f32); }}",
        spec.ty
    );
    line!("    let k = (x * LOG2_E_F32).round() as i32;");
    line!("    let r = {}::from(x - (k as f32) * LN_2_F32);", spec.ty);
    line!(
        "    k_simple_{}(r) * {}::from(scale_pow2_f32(k))",
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
