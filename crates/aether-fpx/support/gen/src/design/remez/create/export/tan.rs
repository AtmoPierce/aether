use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(TanGenerator, "tan");

pub(super) fn write_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    super::common::write_tan_kernel_helpers(w, plan, upper)
}

pub(super) fn write_scalar_body<W: Write>(
    w: &mut W,
    _ctx: &PublicFnCtx<'_>,
    kind: ScalarKind,
) -> io::Result<()> {
    match kind {
        ScalarKind::F64 => {
            writeln!(w, "    if !x.is_finite() {{ return f64::NAN; }}")?;
            writeln!(w, "    let (r, q) = rem_pio2_tan_f64(x);")?;
            writeln!(
                w,
                "    if q & 1 == 0 && r.abs() <= 0.003_f64 {{ return reconstruct_tan_f64(q, k_tan_f64(apply_policy_tan_f64(r))); }}"
            )?;
            writeln!(w, "    let s = crate::trig::sin::sin(r);")?;
            writeln!(w, "    let c = crate::trig::cos::cos(r);")?;
            writeln!(w, "    if q & 1 == 0 {{ s / c }} else {{ -c / s }}")
        }
        ScalarKind::F128 => {
            writeln!(w, "    if !x.is_finite() {{ return f128::NAN; }}")?;
            writeln!(w, "    let (r, q) = rem_pio2_tan_f128(x);")?;
            writeln!(w, "    let s = crate::trig::sin::sinl(r);")?;
            writeln!(w, "    let c = crate::trig::cos::cosl(r);")?;
            writeln!(
                w,
                "    let quotient = if q & 1 == 0 {{ s / c }} else {{ -c / s }};"
            )?;
            writeln!(w, "    let doubled = 2.0_f128 * r;")?;
            writeln!(w, "    let sin_doubled = crate::trig::sin::sinl(doubled);")?;
            writeln!(
                w,
                "    let one_plus_cos_doubled = 1.0_f128 + crate::trig::cos::cosl(doubled);"
            )?;
            writeln!(
                w,
                "    let identity = if q & 1 == 0 {{ sin_doubled / one_plus_cos_doubled }} else {{ -one_plus_cos_doubled / sin_doubled }};"
            )?;
            writeln!(
                w,
                "    if !quotient.is_finite() || !identity.is_finite() {{ quotient }} else {{ quotient + 0.5_f128 * (identity - quotient) }}"
            )
        }
        ScalarKind::F32 => {
            writeln!(w, "    if !x.is_finite() {{ return f32::NAN; }}")?;
            writeln!(w, "    let (r, q) = rem_pio2_tan_f32(x);")?;
            writeln!(
                w,
                "    reconstruct_tan_f32(q, k_tan_f32(apply_policy_tan_f32(r)))"
            )
        }
        ScalarKind::F16 => {
            writeln!(w, "    if !x.is_finite() {{ return f16::NAN; }}")?;
            writeln!(w, "    let (r, q) = rem_pio2_tan_f16(x);")?;
            writeln!(w, "    reconstruct_tan_f16(q, k_tan_f16(r))")
        }
    }
}

pub(super) fn write_mx_body<W: Write>(
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
    writeln!(w, "    {}::from(tanf(x32))", spec.ty)
}
