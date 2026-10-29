use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(AcosGenerator, "acos");

pub(super) fn write_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    super::common::write_simple_kernel_helpers(w, plan, upper, false, false)
}

fn write_approx_scalar_body<W: Write>(
    w: &mut W,
    ctx: &PublicFnCtx<'_>,
    kind: ScalarKind,
) -> io::Result<()> {
    let ident = ctx.ident;
    let upper = ctx.upper;

    macro_rules! line {
		($($arg:tt)*) => {
			writeln!(w, $($arg)*)?;
		};
	}

    match kind {
        ScalarKind::F64 => {
            line!("    if !x.is_finite() {{ return f64::NAN; }}");
            line!("    if x.abs() > 1.0_f64 {{ return f64::NAN; }}");
            line!("    if x > 0.5_f64 {{ return 2.0_f64 * crate::trig::asin::asin(((1.0_f64 - x) * 0.5_f64).sqrt()); }}");
            line!("    if x < -0.5_f64 {{ return core::f64::consts::PI - 2.0_f64 * crate::trig::asin::asin(((1.0_f64 + x) * 0.5_f64).sqrt()); }}");
            line!("    core::f64::consts::FRAC_PI_2 - crate::trig::asin::asin(x)");
        }
        ScalarKind::F128 => {
            line!("    if !x.is_finite() {{ return f128::NAN; }}");
            line!("    if x.abs() > 1.0_f128 {{ return f128::NAN; }}");
            line!("    if x > 0.5_f128 {{ return 2.0_f128 * crate::trig::asin::asinl(((1.0_f128 - x) * 0.5_f128).sqrt()); }}");
            line!("    if x < -0.5_f128 {{ return core::f128::consts::PI - 2.0_f128 * crate::trig::asin::asinl(((1.0_f128 + x) * 0.5_f128).sqrt()); }}");
            line!("    core::f128::consts::FRAC_PI_2 - crate::trig::asin::asinl(x)");
        }
        ScalarKind::F32 => {
            line!("    if !x.is_finite() {{ return f32::NAN; }}");
            line!("    if x.abs() > 1.0_f32 {{ return f32::NAN; }}");
            line!("    if x.is_sign_negative() {{ return core::f32::consts::PI - {ident}f(-x); }}");
            line!("    if x <= {upper}_INTERVAL_MAX_F32 {{ return k_simple_f32(x); }}");
            line!("    let t = ((1.0_f32 - x) * 0.5_f32).sqrt();");
            line!("    2.0_f32 * crate::trig::asin::asinf(t)");
        }
        ScalarKind::F16 => {
            line!("    let x32 = x as f32;");
            line!("    if !x32.is_finite() || x32.abs() > 1.0_f32 {{ return f16::NAN; }}");
            line!("    if x32.is_sign_negative() {{ return (core::f32::consts::PI as f16) - acoss(-x); }}");
            line!("    if x <= ({upper}_INTERVAL_MAX_F32 as f16) {{ return k_simple_f16(x); }}");
            line!("    let t = (((1.0_f32 - x32) * 0.5_f32).sqrt()) as f16;");
            line!("    2.0_f16 * crate::trig::asin::asins(t)");
        }
    }

    Ok(())
}

fn write_approx_mx_body<W: Write>(
    w: &mut W,
    ctx: &PublicFnCtx<'_>,
    fn_name: &str,
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
        "    if !x32.is_finite() {{ return {}::from(f32::NAN); }}",
        spec.ty
    );
    line!(
        "    if x32.abs() > 1.0_f32 {{ return {}::from(f32::NAN); }}",
        spec.ty
    );
    line!("    if x32.is_sign_negative() {{ return {}::from(core::f32::consts::PI) - {fn_name}({}::from(-x32)); }}", spec.ty, spec.ty);
    line!(
        "    if x32 <= {upper}_INTERVAL_MAX_F32 {{ return k_simple_{}(x); }}",
        spec.suffix
    );
    line!(
        "    let t = {}::from(((1.0_f32 - x32) * 0.5_f32).sqrt());",
        spec.ty
    );
    line!(
        "    {}::from(2.0_f32) * crate::trig::asin::asin_{}(t)",
        spec.ty,
        spec.suffix
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
