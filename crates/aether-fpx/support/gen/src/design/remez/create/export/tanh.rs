use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(TanhGenerator, "tanh");

pub(super) fn write_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    // The configured identity plan is a placeholder, so emit the FPX small
    // argument kernel explicitly. Larger arguments are reconstructed through
    // FPX expm1 rather than the platform tanh intrinsic.
    writeln!(w, "#[inline(always)] fn tanh_small_f64(x: f64) -> f64 {{ let z = x * x; x * (1.0_f64 + z * (-3.33333333333333333333e-1_f64 + z * (1.33333333333333333333e-1_f64 + z * (-5.39682539682539682540e-2_f64 + z * (2.18694885361552028219e-2_f64 + z * (-8.86323552990219656886e-3_f64 + z * 3.59212803657248101693e-3_f64)))))) }}")?;
    writeln!(w, "#[inline(always)] fn tanh_small_f32(x: f32) -> f32 {{ let z = x * x; x * (1.0_f32 + z * (-3.3333334e-1_f32 + z * (1.3333334e-1_f32 + z * (-5.3968254e-2_f32 + z * (2.1869488e-2_f32 + z * (-8.863236e-3_f32 + z * 3.592128e-3_f32)))))) }}")?;
    writeln!(w, "#[inline(always)] fn tanh_small_f128(x: f128) -> f128 {{ let z = x * x; x * (1.0_f128 + z * (-3.3333333333333333333333333333333333e-1_f128 + z * (1.3333333333333333333333333333333333e-1_f128 + z * (-5.3968253968253968253968253968253968e-2_f128 + z * (2.1869488536155202821869488536155203e-2_f128 + z * (-8.8632355299021965688632355299021966e-3_f128 + z * 3.5921280365724810169254613699058144e-3_f128)))))) }}")?;
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
            line!("    if x == f64::INFINITY {{ return 1.0_f64; }}");
            line!("    if x == f64::NEG_INFINITY {{ return -1.0_f64; }}");
            line!("    if x == 0.0_f64 {{ return x; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= 0.0625_f64 {{ return tanh_small_f64(x); }}");
            line!("    if ax >= 19.0_f64 {{ return if x.is_sign_negative() {{ -1.0_f64 }} else {{ 1.0_f64 }}; }}");
            line!("    let e = crate::trig::expm1::expm1(2.0_f64 * ax);");
            line!("    let y = e / (e + 2.0_f64);");
            line!("    if x.is_sign_negative() {{ -y }} else {{ y }}");
        }
        ScalarKind::F128 => {
            line!("    if x.is_nan() {{ return f128::NAN; }}");
            line!("    if x == f128::INFINITY {{ return 1.0_f128; }}");
            line!("    if x == f128::NEG_INFINITY {{ return -1.0_f128; }}");
            line!("    if x == 0.0_f128 {{ return x; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= 0.00390625_f128 {{ return tanh_small_f128(x); }}");
            line!("    if ax >= 40.0_f128 {{ return if x.is_sign_negative() {{ -1.0_f128 }} else {{ 1.0_f128 }}; }}");
            line!("    let e = crate::trig::expm1::expm1l(2.0_f128 * ax);");
            line!("    let y = e / (e + 2.0_f128);");
            line!("    if x.is_sign_negative() {{ -y }} else {{ y }}");
        }
        ScalarKind::F32 => {
            line!("    if x.is_nan() {{ return f32::NAN; }}");
            line!("    if x == f32::INFINITY {{ return 1.0_f32; }}");
            line!("    if x == f32::NEG_INFINITY {{ return -1.0_f32; }}");
            line!("    if x == 0.0_f32 {{ return x; }}");
            line!("    let ax = x.abs();");
            line!("    if ax <= 0.25_f32 {{ return tanh_small_f32(x); }}");
            line!("    if ax >= 9.0_f32 {{ return if x.is_sign_negative() {{ -1.0_f32 }} else {{ 1.0_f32 }}; }}");
            line!("    let e = crate::trig::expm1::expm1f(2.0_f32 * ax);");
            line!("    let y = e / (e + 2.0_f32);");
            line!("    if x.is_sign_negative() {{ -y }} else {{ y }}");
        }
        ScalarKind::F16 => {
            line!("    let x32 = x as f32;");
            line!("    if x32.is_nan() {{ return f16::NAN; }}");
            line!("    if x32 == f32::INFINITY {{ return 1.0_f16; }}");
            line!("    if x32 == f32::NEG_INFINITY {{ return -1.0_f16; }}");
            line!("    tanhf(x32) as f16");
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
    line!("    {}::from(tanhf(x32))", spec.ty);

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
