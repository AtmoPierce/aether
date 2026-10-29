use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(SinhGenerator, "sinh");

pub(super) fn write_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    // The configured identity plan is a placeholder, so emit the FPX small
    // argument kernel explicitly. Larger arguments are reconstructed through
    // FPX expm1 rather than the platform sinh intrinsic.
    writeln!(
        w,
        "#[inline(always)] fn sinh_small_f64(x: f64) -> f64 {{ let z = x * x; x * (1.0_f64 + z * (1.66666666666666666667e-1_f64 + z * (8.33333333333333333333e-3_f64 + z * (1.98412698412698412698e-4_f64 + z * (2.75573192239858906526e-6_f64 + z * 2.50521083854417187751e-8_f64))))) }}"
    )?;
    writeln!(
        w,
        "#[inline(always)] fn sinh_small_f32(x: f32) -> f32 {{ let z = x * x; x * (1.0_f32 + z * (1.6666667e-1_f32 + z * (8.333333e-3_f32 + z * (1.984127e-4_f32 + z * (2.755732e-6_f32 + z * 2.505211e-8_f32))))) }}"
    )?;
    writeln!(
        w,
        "#[inline(always)] fn sinh_small_f128(x: f128) -> f128 {{ let z = x * x; x * (1.0_f128 + z * (1.6666666666666666666666666666666667e-1_f128 + z * (8.3333333333333333333333333333333333e-3_f128 + z * (1.9841269841269841269841269841269841e-4_f128 + z * (2.7557319223985890652557319223985891e-6_f128 + z * 2.5052108385441718775052108385441719e-8_f128))))) }}"
    )?;
    writeln!(
        w,
        "#[inline(always)] fn sinh_series_f128(x: f128) -> f128 {{ let z = x * x; let mut sum = x; let mut term = x; let mut n = 1_u32; while n <= 36 {{ let two_n = 2_u32 * n; term *= z / ((two_n as f128) * ((two_n + 1_u32) as f128)); sum += term; n += 1; }} sum }}"
    )?;
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
            line!("    if x == f64::INFINITY {{ return f64::INFINITY; }}");
            line!("    if x == f64::NEG_INFINITY {{ return f64::NEG_INFINITY; }}");
            line!("    if x == 0.0_f64 {{ return x; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= 0.125_f64 {{ return sinh_small_f64(x); }}");
            line!("    let e = crate::trig::expm1::expm1(ax);");
            line!("    let y = e * (e + 2.0_f64) / (2.0_f64 * (e + 1.0_f64));");
            line!("    if x.is_sign_negative() {{ -y }} else {{ y }}");
        }
        ScalarKind::F128 => {
            line!("    if x.is_nan() {{ return f128::NAN; }}");
            line!("    if x == f128::INFINITY {{ return f128::INFINITY; }}");
            line!("    if x == f128::NEG_INFINITY {{ return f128::NEG_INFINITY; }}");
            line!("    if x == 0.0_f128 {{ return x; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= 0.00390625_f128 {{ return sinh_small_f128(x); }}");
            line!("    if ax <= SINH_INTERVAL_MAX_F128 {{ return sinh_series_f128(x); }}");
            line!("    let e = crate::trig::expm1::expm1l(ax);");
            line!("    let y = e * (e + 2.0_f128) / (2.0_f128 * (e + 1.0_f128));");
            line!("    if x.is_sign_negative() {{ -y }} else {{ y }}");
        }
        ScalarKind::F32 => {
            line!("    if x.is_nan() {{ return f32::NAN; }}");
            line!("    if x == f32::INFINITY {{ return f32::INFINITY; }}");
            line!("    if x == f32::NEG_INFINITY {{ return f32::NEG_INFINITY; }}");
            line!("    if x == 0.0_f32 {{ return x; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= 0.5_f32 {{ return sinh_small_f32(x); }}");
            line!("    let e = crate::trig::expm1::expm1f(ax);");
            line!("    let y = e * (e + 2.0_f32) / (2.0_f32 * (e + 1.0_f32));");
            line!("    if x.is_sign_negative() {{ -y }} else {{ y }}");
        }
        ScalarKind::F16 => {
            line!("    let x32 = x as f32;");
            line!("    if x32.is_nan() {{ return f16::NAN; }}");
            line!("    if x32.is_infinite() {{ return x; }}");
            line!("    sinhf(x32) as f16");
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

    line!("    let x32 = f32::from(x);");
    line!(
        "    if !x32.is_finite() {{ return {}::from(f32::NAN); }}",
        spec.ty
    );
    line!("    {}::from(sinhf(x32))", spec.ty);

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
