use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(AtanGenerator, "atan");

pub(super) fn write_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    super::common::write_simple_kernel_helpers(w, plan, upper, false, false)?;

    writeln!(
        w,
        r#"
#[inline(always)]
fn atan_positive_f32(mut x: f32) -> f32 {{
    let reciprocal = x > 1.0_f32;
    if reciprocal {{ x = 1.0_f32 / x; }}

    let mut scale = 1.0_f32;
    while x > 0.125_f32 {{
        x /= 1.0_f32 + (1.0_f32 + x * x).sqrt();
        scale *= 2.0_f32;
    }}

    let z = -(x * x);
    let mut term = x;
    let mut sum = x;
    let mut denominator = 3_u32;
    while denominator <= 19 {{
        term *= z;
        sum += term / denominator as f32;
        denominator += 2;
    }}
    let y = scale * sum;
    if reciprocal {{ core::f32::consts::FRAC_PI_2 - y }} else {{ y }}
}}

#[inline(always)]
fn atan_positive_f64(mut x: f64) -> f64 {{
    let reciprocal = x > 1.0_f64;
    if reciprocal {{ x = 1.0_f64 / x; }}

    let mut scale = 1.0_f64;
    while x > 0.125_f64 {{
        x /= 1.0_f64 + (1.0_f64 + x * x).sqrt();
        scale *= 2.0_f64;
    }}

    let z = -(x * x);
    let mut term = x;
    let mut sum = x;
    let mut denominator = 3_u32;
    while denominator <= 35 {{
        term *= z;
        sum += term / denominator as f64;
        denominator += 2;
    }}
    let y = scale * sum;
    if reciprocal {{ core::f64::consts::FRAC_PI_2 - y }} else {{ y }}
}}

#[inline(always)]
fn atan_positive_f128(mut x: f128) -> f128 {{
    let reciprocal = x > 1.0_f128;
    if reciprocal {{ x = 1.0_f128 / x; }}

    let mut scale = 1.0_f128;
    while x > 0.125_f128 {{
        x /= 1.0_f128 + (1.0_f128 + x * x).sqrt();
        scale *= 2.0_f128;
    }}

    let z = -(x * x);
    let mut denominator = 71_u32;
    let mut poly = 1.0_f128 / denominator as f128;
    while denominator > 3 {{
        denominator -= 2;
        poly = poly.mul_add(z, 1.0_f128 / denominator as f128);
    }}
    let sum = poly.mul_add(z, 1.0_f128) * x;
    let y = scale * sum;
    if reciprocal {{ core::f128::consts::FRAC_PI_2 - y }} else {{ y }}
}}
"#
    )
}

fn write_approx_scalar_body<W: Write>(
    w: &mut W,
    ctx: &PublicFnCtx<'_>,
    kind: ScalarKind,
) -> io::Result<()> {
    let _ = ctx;

    macro_rules! line {
		($($arg:tt)*) => {
			writeln!(w, $($arg)*)?;
		};
	}

    match kind {
        ScalarKind::F64 => {
            line!("    if x.is_nan() {{ return f64::NAN; }}");
            line!("    if x == f64::INFINITY {{ return core::f64::consts::FRAC_PI_2; }}");
            line!("    if x == f64::NEG_INFINITY {{ return -core::f64::consts::FRAC_PI_2; }}");
            line!("    if x == 0.0_f64 {{ return x; }}");
            line!("    let y = atan_positive_f64(x.abs());");
            line!("    if x.is_sign_negative() {{ -y }} else {{ y }}");
        }
        ScalarKind::F128 => {
            line!("    if x.is_nan() {{ return f128::NAN; }}");
            line!("    if x == f128::INFINITY {{ return core::f128::consts::FRAC_PI_2; }}");
            line!("    if x == f128::NEG_INFINITY {{ return -core::f128::consts::FRAC_PI_2; }}");
            line!("    if x == 0.0_f128 {{ return x; }}");
            line!("    let y = atan_positive_f128(x.abs());");
            line!("    if x.is_sign_negative() {{ -y }} else {{ y }}");
        }
        ScalarKind::F32 => {
            line!("    if x.is_nan() {{ return f32::NAN; }}");
            line!("    if x == f32::INFINITY {{ return core::f32::consts::FRAC_PI_2; }}");
            line!("    if x == f32::NEG_INFINITY {{ return -core::f32::consts::FRAC_PI_2; }}");
            line!("    if x == 0.0_f32 {{ return x; }}");
            line!("    let y = atan_positive_f32(x.abs());");
            line!("    if x.is_sign_negative() {{ -y }} else {{ y }}");
        }
        ScalarKind::F16 => {
            line!("    let x32 = x as f32;");
            line!("    if x32.is_nan() {{ return f16::NAN; }}");
            line!("    if x32 == f32::INFINITY {{ return core::f32::consts::FRAC_PI_2 as f16; }}");
            line!("    if x32 == f32::NEG_INFINITY {{ return -(core::f32::consts::FRAC_PI_2 as f16); }}");
            line!("    if x32 == 0.0_f32 {{ return x; }}");
            line!("    atanf(x32) as f16");
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
    let _ = (ctx, fn_name);

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
        "    if x32 == f32::INFINITY {{ return {}::from(core::f32::consts::FRAC_PI_2); }}",
        spec.ty
    );
    line!(
        "    if x32 == f32::NEG_INFINITY {{ return {}::from(-core::f32::consts::FRAC_PI_2); }}",
        spec.ty
    );
    line!("    if x32 == 0.0_f32 {{ return x; }}");
    line!("    {}::from(atanf(x32))", spec.ty);

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
