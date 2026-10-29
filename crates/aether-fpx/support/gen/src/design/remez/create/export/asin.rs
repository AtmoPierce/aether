use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(AsinGenerator, "asin");

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
            line!("    let ax = x.abs();");
            line!("    if ax > 1.0_f64 {{ return f64::NAN; }}");
            line!("    if x == 0.0_f64 {{ return x; }}");
            line!("    if ax == 1.0_f64 {{ return if x.is_sign_negative() {{ -core::f64::consts::FRAC_PI_2 }} else {{ core::f64::consts::FRAC_PI_2 }}; }}");
            line!("    crate::trig::atan::atan(x / ((1.0_f64 - x) * (1.0_f64 + x)).sqrt())");
        }
        ScalarKind::F128 => {
            line!("    if !x.is_finite() {{ return f128::NAN; }}");
            line!("    let ax = x.abs();");
            line!("    if ax > 1.0_f128 {{ return f128::NAN; }}");
            line!("    if x == 0.0_f128 {{ return x; }}");
            line!("    if ax == 1.0_f128 {{ return if x.is_sign_negative() {{ -core::f128::consts::FRAC_PI_2 }} else {{ core::f128::consts::FRAC_PI_2 }}; }}");
            line!("    let y = if ax <= 0.75_f128 {{");
            line!("        crate::trig::atan::atanl(ax / (-ax).mul_add(ax, 1.0_f128).sqrt())");
            line!("    }} else {{");
            line!("        let t = ((1.0_f128 - ax) * 0.5_f128).sqrt();");
            line!("        let inner = crate::trig::atan::atanl(t / (-t).mul_add(t, 1.0_f128).sqrt());");
            line!("        core::f128::consts::FRAC_PI_2 - 2.0_f128 * inner");
            line!("    }};");
            line!("    if x.is_sign_negative() {{ -y }} else {{ y }}");
        }
        ScalarKind::F32 => {
            line!("    if !x.is_finite() {{ return f32::NAN; }}");
            line!("    let ax = x.abs();");
            line!("    if ax > 1.0_f32 {{ return f32::NAN; }}");
            line!("    if x == 0.0_f32 {{ return x; }}");
            line!("    let y = if ax <= {upper}_INTERVAL_MAX_F32 {{ k_simple_f32(ax) }} else {{ let t = ((1.0_f32 - ax) * 0.5_f32).sqrt(); core::f32::consts::FRAC_PI_2 - 2.0_f32 * {ident}f(t) }};");
            line!("    if x.is_sign_negative() {{ -y }} else {{ y }}");
        }
        ScalarKind::F16 => {
            line!("    let x32 = x as f32;");
            line!("    if !x32.is_finite() || x32.abs() > 1.0_f32 {{ return f16::NAN; }}");
            line!("    if x32 == 0.0_f32 {{ return x; }}");
            line!("    let ax = x.abs();");
            line!("    let y = if ax <= ({upper}_INTERVAL_MAX_F32 as f16) {{ k_simple_f16(ax) }} else {{ let t = (((1.0_f32 - x32.abs()) * 0.5_f32).sqrt()) as f16; (core::f32::consts::FRAC_PI_2 as f16) - 2.0_f16 * asins(t) }};");
            line!("    if x32.is_sign_negative() {{ -y }} else {{ y }}");
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
    line!("    let ax32 = x32.abs();");
    line!(
        "    if ax32 > 1.0_f32 {{ return {}::from(f32::NAN); }}",
        spec.ty
    );
    line!("    if x32 == 0.0_f32 {{ return x; }}");
    line!("    let y = if ax32 <= {upper}_INTERVAL_MAX_F32 {{");
    line!("        k_simple_{}({}::from(ax32))", spec.suffix, spec.ty);
    line!("    }} else {{");
    line!(
        "        let t = {}::from(((1.0_f32 - ax32) * 0.5_f32).sqrt());",
        spec.ty
    );
    line!(
        "        {}::from(core::f32::consts::FRAC_PI_2) - {}::from(2.0_f32) * {fn_name}(t)",
        spec.ty,
        spec.ty
    );
    line!("    }};");
    line!("    if x32.is_sign_negative() {{ -y }} else {{ y }}");

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
