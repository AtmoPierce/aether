use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(LogGenerator, "log");

pub(super) fn write_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    super::common::write_simple_kernel_helpers(w, plan, upper, true, true)
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
        "#[inline(always)] fn log_reduced_f64(x: f64) -> f64 {{ let (m, e) = reduce_log_f64(x); let base = k_simple_f64(m); let ef = e as f64; (base + ef * LN_2_LO_F64) + ef * LN_2_HI_F64 }}"
    );
    line!(
        "#[inline(always)] fn log_reduced_f32(x: f32) -> f32 {{ let (m, e) = reduce_log_f32(x); let base = k_simple_f32(m); let ef = e as f64; ((base as f64 + ef * LN_2_LO_F64) + ef * LN_2_HI_F64) as f32 }}"
    );
    line!(
        "#[inline(always)] fn log_reduced_f128(x: f128) -> f128 {{ let (m, e) = reduce_log_f128(x); let base = k_simple_f128(m); let ef = e as f128; ((base + ef * LN_2_LO_F128) + ef * LN_2_MID_F128) + ef * LN_2_HI_F128 }}"
    );
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
            line!("    if x == 0.0_f64 {{ return f64::NEG_INFINITY; }}");
            line!("    if x < 0.0_f64 {{ return f64::NAN; }}");
            line!("    if x == f64::INFINITY {{ return f64::INFINITY; }}");
            line!("    if x == 1.0_f64 {{ return 0.0_f64; }}");
            line!("    log_reduced_f64(x)");
        }
        ScalarKind::F128 => {
            line!("    if x.is_nan() {{ return f128::NAN; }}");
            line!("    if x == 0.0_f128 {{ return f128::NEG_INFINITY; }}");
            line!("    if x < 0.0_f128 {{ return f128::NAN; }}");
            line!("    if x == f128::INFINITY {{ return f128::INFINITY; }}");
            line!("    if x == 1.0_f128 {{ return 0.0_f128; }}");
            line!("    log_reduced_f128(x)");
        }
        ScalarKind::F32 => {
            line!("    if x.is_nan() {{ return f32::NAN; }}");
            line!("    if x == 0.0_f32 {{ return f32::NEG_INFINITY; }}");
            line!("    if x < 0.0_f32 {{ return f32::NAN; }}");
            line!("    if x == f32::INFINITY {{ return f32::INFINITY; }}");
            line!("    if x == 1.0_f32 {{ return 0.0_f32; }}");
            line!("    log_reduced_f32(x)");
        }
        ScalarKind::F16 => {
            line!("    let x = x as f32;");
            line!("    if x.is_nan() {{ return f16::NAN; }}");
            line!("    if x == 0.0_f32 {{ return f16::NEG_INFINITY; }}");
            line!("    if x < 0.0_f32 {{ return f16::NAN; }}");
            line!("    if x == f32::INFINITY {{ return f16::INFINITY; }}");
            line!("    let (m, e) = reduce_log_f32(x);");
            line!("    let base = k_simple_f16(m as f16) as f32;");
            line!("    (base + (e as f32) * LN_2_F32) as f16");
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
        "    if x32.is_nan() {{ return {}::from(f32::NAN); }}",
        spec.ty
    );
    line!(
        "    if x32 == 0.0_f32 {{ return {}::from(f32::NEG_INFINITY); }}",
        spec.ty
    );
    line!(
        "    if x32 < 0.0_f32 {{ return {}::from(f32::NAN); }}",
        spec.ty
    );
    line!(
        "    if x32 == f32::INFINITY {{ return {}::from(f32::INFINITY); }}",
        spec.ty
    );
    line!("    {}::from(log_reduced_f32(x32))", spec.ty);

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn public_log_scalars_do_not_bypass_the_generated_kernel() {
        let ctx = PublicFnCtx {
            ident: "log",
            upper: "LOG",
        };

        for kind in [
            ScalarKind::F16,
            ScalarKind::F32,
            ScalarKind::F64,
            ScalarKind::F128,
        ] {
            let mut body = Vec::new();
            write_approx_scalar_body(&mut body, &ctx, kind).unwrap();
            let body = String::from_utf8(body).unwrap();

            assert!(body.contains("k_simple") || body.contains("log_reduced"));
            assert!(!body.contains("log1p"), "generated log bypass:\n{body}");
            assert!(
                !body.contains("EPSILON"),
                "generated log near-one bypass:\n{body}"
            );
        }
    }
}
