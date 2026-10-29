use std::io::{self, Write};

use super::common::{PublicFnCtx, ScalarKind};

super::common::define_function_generator!(CosGenerator, "cos");

pub(super) fn write_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    super::common::write_cos_kernel_helpers(w, plan, upper)
}

pub(super) fn write_scalar_body<W: Write>(
    w: &mut W,
    _ctx: &PublicFnCtx<'_>,
    kind: ScalarKind,
) -> io::Result<()> {
    match kind {
        ScalarKind::F64 => {
            writeln!(w, "    if !x.is_finite() {{ return f64::NAN; }}")?;
            writeln!(w, "    let (r, q) = rem_pio2_f64(x);")?;
            writeln!(w, "    match q {{ 0 => k_cos_f64(r), 1 => -k_sin_f64(r), 2 => -k_cos_f64(r), 3 => k_sin_f64(r), _ => unreachable!() }}")
        }
        ScalarKind::F128 => {
            writeln!(w, "    if !x.is_finite() {{ return f128::NAN; }}")?;
            writeln!(w, "    let (r, q) = rem_pio2_f128(x);")?;
            writeln!(w, "    match q {{ 0 => k_cos_f128(r), 1 => -k_sin_f128(r), 2 => -k_cos_f128(r), 3 => k_sin_f128(r), _ => unreachable!() }}")
        }
        ScalarKind::F32 => {
            writeln!(w, "    if !x.is_finite() {{ return f32::NAN; }}")?;
            writeln!(w, "    let (r, q) = rem_pio2_f32(x);")?;
            writeln!(w, "    match q {{ 0 => k_cos_f32(r), 1 => -k_sin_f32(r), 2 => -k_cos_f32(r), 3 => k_sin_f32(r), _ => unreachable!() }}")
        }
        ScalarKind::F16 => {
            writeln!(w, "    if !x.is_finite() {{ return f16::NAN; }}")?;
            writeln!(w, "    let (r, q) = rem_pio2_f16(x);")?;
            writeln!(w, "    match q {{ 0 => k_cos_f16(r), 1 => -k_sin_f16(r), 2 => -k_cos_f16(r), 3 => k_sin_f16(r), _ => unreachable!() }}")
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
    writeln!(w, "    let (r, q) = rem_pio2_f32(x32);")?;
    writeln!(w, "    let r = {}::from(r);", spec.ty)?;
    writeln!(w, "    match q {{ 0 => k_cos_{}(r), 1 => -k_sin_{}(r), 2 => -k_cos_{}(r), 3 => k_sin_{}(r), _ => unreachable!() }}", spec.suffix, spec.suffix, spec.suffix, spec.suffix)
}
