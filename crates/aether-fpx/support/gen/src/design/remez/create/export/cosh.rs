use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(CoshGenerator, "cosh");

pub(super) fn write_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    // The configured identity plan is a placeholder, so emit the FPX small
    // argument kernel explicitly. Larger arguments are reconstructed through
    // FPX expm1 rather than the platform cosh intrinsic.
    writeln!(w, "#[inline(always)] fn cosh_small_f64(x: f64) -> f64 {{ let z = x * x; 1.0_f64 + z * (5.0e-1_f64 + z * (4.16666666666666666667e-2_f64 + z * (1.38888888888888888889e-3_f64 + z * (2.48015873015873015873e-5_f64 + z * (2.75573192239858906526e-7_f64 + z * 2.08767569878680989792e-9_f64))))) }}")?;
    writeln!(w, "#[inline(always)] fn cosh_small_f32(x: f32) -> f32 {{ let z = x * x; 1.0_f32 + z * (5.0e-1_f32 + z * (4.1666668e-2_f32 + z * (1.3888889e-3_f32 + z * (2.4801588e-5_f32 + z * (2.755732e-7_f32 + z * 2.087676e-9_f32))))) }}")?;
    writeln!(w, "#[inline(always)] fn cosh_small_f128(x: f128) -> f128 {{ let z = x * x; 1.0_f128 + z * (5.0e-1_f128 + z * (4.1666666666666666666666666666666667e-2_f128 + z * (1.3888888888888888888888888888888889e-3_f128 + z * (2.4801587301587301587301587301587302e-5_f128 + z * (2.7557319223985890652557319223985891e-7_f128 + z * 2.0876756987868098979210090321201432e-9_f128))))) }}")?;
    writeln!(w, "#[inline(always)] fn cosh_series_f128(x: f128) -> f128 {{ let z = x * x; let mut sum = 1.0_f128; let mut term = 1.0_f128; let mut n = 1_u32; while n <= 36 {{ let two_n = 2_u32 * n; term *= z / (((two_n - 1_u32) as f128) * (two_n as f128)); sum += term; n += 1; }} sum }}")?;
    let _ = (plan, upper);
    Ok(())
}

fn write_approx_scalar_body<W: Write>(
    w: &mut W,
    ctx: &PublicFnCtx<'_>,
    kind: ScalarKind,
) -> io::Result<()> {
    macro_rules! line {
		($($arg:tt)*) => {
			writeln!(w, $($arg)*)?;
		};
	}

    let _ = ctx;
    match kind {
        ScalarKind::F64 => {
            line!("    if x.is_nan() {{ return f64::NAN; }}");
            line!("    if x.is_infinite() {{ return f64::INFINITY; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= 0.125_f64 {{ return cosh_small_f64(x); }}");
            line!("    let e = crate::trig::expm1::expm1(ax);");
            line!("    1.0_f64 + e * e / (2.0_f64 * (e + 1.0_f64))");
        }
        ScalarKind::F128 => {
            line!("    if x.is_nan() {{ return f128::NAN; }}");
            line!("    if x.is_infinite() {{ return f128::INFINITY; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= 0.00390625_f128 {{ return cosh_small_f128(x); }}");
            line!("    if ax <= COSH_INTERVAL_MAX_F128 {{ return cosh_series_f128(x); }}");
            line!("    let e = crate::trig::expm1::expm1l(ax);");
            line!("    1.0_f128 + e * e / (2.0_f128 * (e + 1.0_f128))");
        }
        ScalarKind::F32 => {
            line!("    if x.is_nan() {{ return f32::NAN; }}");
            line!("    if x.is_infinite() {{ return f32::INFINITY; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= 0.5_f32 {{ return cosh_small_f32(x); }}");
            line!("    let e = crate::trig::expm1::expm1f(ax);");
            line!("    1.0_f32 + e * e / (2.0_f32 * (e + 1.0_f32))");
        }
        ScalarKind::F16 => {
            line!("    let x32 = x as f32;");
            line!("    if x32.is_nan() {{ return f16::NAN; }}");
            line!("    if x32.is_infinite() {{ return f16::INFINITY; }}");
            line!("    coshf(x32) as f16");
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
    writeln!(w, "    let x32 = f32::from(x);")?;
    writeln!(
        w,
        "    if !x32.is_finite() {{ return {}::from(f32::NAN); }}",
        spec.ty
    )?;
    writeln!(w, "    {}::from(coshf(x32))", spec.ty)?;
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
