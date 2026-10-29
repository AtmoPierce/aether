use std::io::{self, Write};

#[derive(Clone, Copy)]
pub(super) enum ScalarKind {
    F64,
    F128,
    F32,
    F16,
}

#[derive(Clone, Copy)]
pub(super) struct PublicFnCtx<'a> {
    pub ident: &'a str,
    pub upper: &'a str,
}

#[derive(Clone, Copy)]
pub(super) enum HelperTarget<'a> {
    Scalar(ScalarKind),
    Mx {
        fn_name: &'a str,
        spec: &'a super::MxSpec,
    },
}

#[allow(dead_code)]
#[derive(Clone, Copy)]
pub(super) struct FunctionHelpers {
    pub func_name: &'static str,
}

impl FunctionHelpers {
    pub const fn new(func_name: &'static str) -> Self {
        Self { func_name }
    }

    pub fn ident(&self) -> String {
        super::sanitize_ident(self.func_name)
    }

    pub fn upper(&self) -> String {
        self.ident().to_uppercase()
    }
}

macro_rules! define_function_generator {
    ($struct_name:ident, $func_name:literal) => {
        #[derive(Clone, Copy, Default)]
        pub(super) struct $struct_name;

        pub(super) const GENERATOR: $struct_name = $struct_name;

        impl $struct_name {
            pub(super) fn helpers(&self) -> super::common::FunctionHelpers {
                super::common::FunctionHelpers::new($func_name)
            }

            pub(super) fn write_rust_body<W: std::io::Write>(
                &self,
                w: &mut W,
                ctx: &super::common::PublicFnCtx<'_>,
                target: super::common::HelperTarget<'_>,
            ) -> std::io::Result<()> {
                match target {
                    super::common::HelperTarget::Scalar(kind) => write_scalar_body(w, ctx, kind),
                    super::common::HelperTarget::Mx { fn_name, spec } => {
                        write_mx_body(w, ctx, fn_name, spec)
                    }
                }
            }

            pub(super) fn rust_generate(
                &self,
                func_name: &str,
                func: &super::FunctionSpec,
                plan: &super::FunctionPlan,
            ) {
                super::write_function_to_rs_shared(func_name, func, plan);
            }

            pub(super) fn c_generate<const N: usize, const K: usize, F>(
                &self,
                rz: &crate::remez::RemezApprox<f128, N, K>,
                func_name: &str,
                f_x: F,
                interval_x: (f128, f128),
                max_error: f128,
                file_prefix: &str,
            ) where
                F: Fn(f128) -> f128,
            {
                super::write_coeffs_to_c_shared(
                    rz,
                    func_name,
                    f_x,
                    interval_x,
                    max_error,
                    file_prefix,
                );
            }

            pub(super) fn cuda_generate(&self, export: &super::FunctionRemezExport) -> String {
                super::render_gpu_backend_file(export, "cuda")
            }

            pub(super) fn rocm_generate(&self, export: &super::FunctionRemezExport) -> String {
                super::render_gpu_backend_file(export, "rocm")
            }
        }
    };
}

pub(super) use define_function_generator;

#[derive(Clone, Copy)]
pub(super) enum FunctionGenerator {
    Sin,
    Cos,
    Tan,
    Exp,
    Expm1,
    Tanh,
    Sinh,
    Cosh,
    Log,
    Log1p,
    Asin,
    Acos,
    Atan,
}

impl FunctionGenerator {
    pub(super) fn helpers(&self) -> FunctionHelpers {
        match self {
            Self::Sin => super::sin::GENERATOR.helpers(),
            Self::Cos => super::cos::GENERATOR.helpers(),
            Self::Tan => super::tan::GENERATOR.helpers(),
            Self::Exp => super::exp::GENERATOR.helpers(),
            Self::Expm1 => super::expm1::GENERATOR.helpers(),
            Self::Tanh => super::tanh::GENERATOR.helpers(),
            Self::Sinh => super::sinh::GENERATOR.helpers(),
            Self::Cosh => super::cosh::GENERATOR.helpers(),
            Self::Log => super::log::GENERATOR.helpers(),
            Self::Log1p => super::log1p::GENERATOR.helpers(),
            Self::Asin => super::asin::GENERATOR.helpers(),
            Self::Acos => super::acos::GENERATOR.helpers(),
            Self::Atan => super::atan::GENERATOR.helpers(),
        }
    }

    pub(super) fn write_rust_body<W: Write>(
        &self,
        w: &mut W,
        ctx: &PublicFnCtx<'_>,
        target: HelperTarget<'_>,
    ) -> io::Result<()> {
        match self {
            Self::Sin => super::sin::GENERATOR.write_rust_body(w, ctx, target),
            Self::Cos => super::cos::GENERATOR.write_rust_body(w, ctx, target),
            Self::Tan => super::tan::GENERATOR.write_rust_body(w, ctx, target),
            Self::Exp => super::exp::GENERATOR.write_rust_body(w, ctx, target),
            Self::Expm1 => super::expm1::GENERATOR.write_rust_body(w, ctx, target),
            Self::Tanh => super::tanh::GENERATOR.write_rust_body(w, ctx, target),
            Self::Sinh => super::sinh::GENERATOR.write_rust_body(w, ctx, target),
            Self::Cosh => super::cosh::GENERATOR.write_rust_body(w, ctx, target),
            Self::Log => super::log::GENERATOR.write_rust_body(w, ctx, target),
            Self::Log1p => super::log1p::GENERATOR.write_rust_body(w, ctx, target),
            Self::Asin => super::asin::GENERATOR.write_rust_body(w, ctx, target),
            Self::Acos => super::acos::GENERATOR.write_rust_body(w, ctx, target),
            Self::Atan => super::atan::GENERATOR.write_rust_body(w, ctx, target),
        }
    }

    pub(super) fn rust_generate(
        &self,
        func_name: &str,
        func: &super::FunctionSpec,
        plan: &super::FunctionPlan,
    ) {
        match self {
            Self::Sin => super::sin::GENERATOR.rust_generate(func_name, func, plan),
            Self::Cos => super::cos::GENERATOR.rust_generate(func_name, func, plan),
            Self::Tan => super::tan::GENERATOR.rust_generate(func_name, func, plan),
            Self::Exp => super::exp::GENERATOR.rust_generate(func_name, func, plan),
            Self::Expm1 => super::expm1::GENERATOR.rust_generate(func_name, func, plan),
            Self::Tanh => super::tanh::GENERATOR.rust_generate(func_name, func, plan),
            Self::Sinh => super::sinh::GENERATOR.rust_generate(func_name, func, plan),
            Self::Cosh => super::cosh::GENERATOR.rust_generate(func_name, func, plan),
            Self::Log => super::log::GENERATOR.rust_generate(func_name, func, plan),
            Self::Log1p => super::log1p::GENERATOR.rust_generate(func_name, func, plan),
            Self::Asin => super::asin::GENERATOR.rust_generate(func_name, func, plan),
            Self::Acos => super::acos::GENERATOR.rust_generate(func_name, func, plan),
            Self::Atan => super::atan::GENERATOR.rust_generate(func_name, func, plan),
        }
    }

    pub(super) fn write_helpers<W: Write>(&self, w: &mut W) -> io::Result<()> {
        match self {
            Self::Exp => super::exp::write_helpers(w),
            Self::Expm1 => super::expm1::write_helpers(w),
            Self::Log => super::log::write_helpers(w),
            Self::Log1p => super::log1p::write_helpers(w),
            _ => Ok(()),
        }
    }

    pub(super) fn write_kernel_helpers<W: Write>(
        &self,
        w: &mut W,
        plan: &super::FunctionPlan,
        upper: &str,
    ) -> io::Result<()> {
        match self {
            Self::Sin => super::sin::write_kernel_helpers(w, plan, upper),
            Self::Cos => super::cos::write_kernel_helpers(w, plan, upper),
            Self::Tan => super::tan::write_kernel_helpers(w, plan, upper),
            Self::Exp => super::exp::write_kernel_helpers(w, plan, upper),
            Self::Expm1 => super::expm1::write_kernel_helpers(w, plan, upper),
            Self::Tanh => super::tanh::write_kernel_helpers(w, plan, upper),
            Self::Sinh => super::sinh::write_kernel_helpers(w, plan, upper),
            Self::Cosh => super::cosh::write_kernel_helpers(w, plan, upper),
            Self::Log => super::log::write_kernel_helpers(w, plan, upper),
            Self::Log1p => super::log1p::write_kernel_helpers(w, plan, upper),
            Self::Asin => super::asin::write_kernel_helpers(w, plan, upper),
            Self::Acos => super::acos::write_kernel_helpers(w, plan, upper),
            Self::Atan => super::atan::write_kernel_helpers(w, plan, upper),
        }
    }

    pub(super) fn c_generate<const N: usize, const K: usize, F>(
        &self,
        rz: &crate::remez::RemezApprox<f128, N, K>,
        func_name: &str,
        f_x: F,
        interval_x: (f128, f128),
        max_error: f128,
        file_prefix: &str,
    ) where
        F: Fn(f128) -> f128,
    {
        match self {
            Self::Sin => super::sin::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Cos => super::cos::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Tan => super::tan::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Exp => super::exp::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Expm1 => super::expm1::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Tanh => super::tanh::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Sinh => super::sinh::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Cosh => super::cosh::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Log => super::log::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Log1p => super::log1p::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Asin => super::asin::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Acos => super::acos::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
            Self::Atan => super::atan::GENERATOR.c_generate(
                rz,
                func_name,
                f_x,
                interval_x,
                max_error,
                file_prefix,
            ),
        }
    }

    pub(super) fn cuda_generate(&self, export: &super::FunctionRemezExport) -> String {
        match self {
            Self::Sin => super::sin::GENERATOR.cuda_generate(export),
            Self::Cos => super::cos::GENERATOR.cuda_generate(export),
            Self::Tan => super::tan::GENERATOR.cuda_generate(export),
            Self::Exp => super::exp::GENERATOR.cuda_generate(export),
            Self::Expm1 => super::expm1::GENERATOR.cuda_generate(export),
            Self::Tanh => super::tanh::GENERATOR.cuda_generate(export),
            Self::Sinh => super::sinh::GENERATOR.cuda_generate(export),
            Self::Cosh => super::cosh::GENERATOR.cuda_generate(export),
            Self::Log => super::log::GENERATOR.cuda_generate(export),
            Self::Log1p => super::log1p::GENERATOR.cuda_generate(export),
            Self::Asin => super::asin::GENERATOR.cuda_generate(export),
            Self::Acos => super::acos::GENERATOR.cuda_generate(export),
            Self::Atan => super::atan::GENERATOR.cuda_generate(export),
        }
    }

    pub(super) fn rocm_generate(&self, export: &super::FunctionRemezExport) -> String {
        match self {
            Self::Sin => super::sin::GENERATOR.rocm_generate(export),
            Self::Cos => super::cos::GENERATOR.rocm_generate(export),
            Self::Tan => super::tan::GENERATOR.rocm_generate(export),
            Self::Exp => super::exp::GENERATOR.rocm_generate(export),
            Self::Expm1 => super::expm1::GENERATOR.rocm_generate(export),
            Self::Tanh => super::tanh::GENERATOR.rocm_generate(export),
            Self::Sinh => super::sinh::GENERATOR.rocm_generate(export),
            Self::Cosh => super::cosh::GENERATOR.rocm_generate(export),
            Self::Log => super::log::GENERATOR.rocm_generate(export),
            Self::Log1p => super::log1p::GENERATOR.rocm_generate(export),
            Self::Asin => super::asin::GENERATOR.rocm_generate(export),
            Self::Acos => super::acos::GENERATOR.rocm_generate(export),
            Self::Atan => super::atan::GENERATOR.rocm_generate(export),
        }
    }
}

pub(super) fn generator_for(func_name: &str) -> io::Result<FunctionGenerator> {
    match func_name {
        "sin" => Ok(FunctionGenerator::Sin),
        "cos" => Ok(FunctionGenerator::Cos),
        "tan" => Ok(FunctionGenerator::Tan),
        "exp" => Ok(FunctionGenerator::Exp),
        "expm1" => Ok(FunctionGenerator::Expm1),
        "tanh" => Ok(FunctionGenerator::Tanh),
        "sinh" => Ok(FunctionGenerator::Sinh),
        "cosh" => Ok(FunctionGenerator::Cosh),
        "log" => Ok(FunctionGenerator::Log),
        "log1p" => Ok(FunctionGenerator::Log1p),
        "asin" => Ok(FunctionGenerator::Asin),
        "acos" => Ok(FunctionGenerator::Acos),
        "atan" => Ok(FunctionGenerator::Atan),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("unsupported function {func_name}"),
        )),
    }
}

pub(super) fn write_scalar_public_body<W: Write>(
    w: &mut W,
    func_name: &str,
    ctx: &PublicFnCtx<'_>,
    kind: ScalarKind,
) -> io::Result<()> {
    generator_for(func_name)?.write_rust_body(w, ctx, HelperTarget::Scalar(kind))
}

pub(super) fn write_mx_public_body<W: Write>(
    w: &mut W,
    func_name: &str,
    ctx: &PublicFnCtx<'_>,
    fn_name: &str,
    spec: &super::MxSpec,
) -> io::Result<()> {
    generator_for(func_name)?.write_rust_body(w, ctx, HelperTarget::Mx { fn_name, spec })
}

fn write_trig_dispatch_aliases<W: Write>(w: &mut W) -> io::Result<()> {
    writeln!(w, "#[inline(always)] fn k_sin_f64(x: f64) -> f64 {{ crate::trig::sin::sin_kernel_estrin_f64(x) }}")?;
    writeln!(w, "#[inline(always)] fn k_cos_f64(x: f64) -> f64 {{ crate::trig::cos::cos_kernel_estrin_f64(x) }}")?;
    writeln!(
        w,
        "#[inline(always)] fn k_sin_f32(x: f32) -> f32 {{ crate::trig::sin::sin_kernel_f32(x) }}"
    )?;
    writeln!(
        w,
        "#[inline(always)] fn k_cos_f32(x: f32) -> f32 {{ crate::trig::cos::cos_kernel_f32(x) }}"
    )?;
    writeln!(
        w,
        "#[inline(always)] fn k_sin_f16(x: f16) -> f16 {{ crate::trig::sin::sin_kernel_f16(x) }}"
    )?;
    writeln!(
        w,
        "#[inline(always)] fn k_cos_f16(x: f16) -> f16 {{ crate::trig::cos::cos_kernel_f16(x) }}"
    )?;
    writeln!(w, "#[inline(always)] fn k_sin_f128(x: f128) -> f128 {{ crate::trig::sin::sin_kernel_f128(x) }}")?;
    writeln!(w, "#[inline(always)] fn k_cos_f128(x: f128) -> f128 {{ crate::trig::cos::cos_kernel_f128(x) }}")?;
    writeln!(w)?;
    Ok(())
}

pub(crate) fn write_sin_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    if plan.segments.len() == 1 {
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn sin_kernel_f64(x: f64) -> f64 {{ let z = x * x; x * eval_horner_f64_slice(z, &COEFFS64) }}"
		)?;
        writeln!(
            w,
            "#[inline(always)] pub(crate) fn sin_kernel_estrin_f64(x: f64) -> f64 {{"
        )?;
        super::write_straight_estrin_z_f64(
            w,
            "    ",
            "x",
            "COEFFS64",
            plan.segments[0].c64.as_slice().len(),
            super::TrigKernel::Sin,
        )?;
        writeln!(w, "}}")?;
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn sin_kernel_f32(x: f32) -> f32 {{ let z = x * x; x * eval_horner_f32_slice(z, &COEFFS32) }}"
		)?;
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn sin_kernel_estrin_f32(x: f32) -> f32 {{ let z = x * x; x * eval_estrin_f32_slice(z, &COEFFS32) }}"
		)?;
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn sin_kernel_f16(x: f16) -> f16 {{ let z = x * x; x * eval_horner_f16_slice(z, &COEFFS16) }}"
		)?;
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn sin_kernel_f128(x: f128) -> f128 {{ let z = x * x; x * eval_horner_f128_slice(z, &COEFFS128) }}"
		)?;
    } else {
        writeln!(
            w,
            "#[inline(always)] pub(crate) fn sin_kernel_f64(x: f64) -> f64 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{ let z = x * x; x * eval_horner_f64_slice(z, &COEFFS64_SEG{idx}) }}")?;
        }
        let last = plan.segments.len() - 1;
        writeln!(
            w,
            "    else {{ let z = x * x; x * eval_horner_f64_slice(z, &COEFFS64_SEG{last}) }}"
        )?;
        writeln!(w, "}}")?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn sin_kernel_estrin_f64(x: f64) -> f64 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{")?;
            super::write_straight_estrin_z_f64(
                w,
                "        ",
                "x",
                &format!("COEFFS64_SEG{idx}"),
                plan.segments[idx].c64.as_slice().len(),
                super::TrigKernel::Sin,
            )?;
            writeln!(w, "    }}")?;
        }
        writeln!(w, "    else {{")?;
        super::write_straight_estrin_z_f64(
            w,
            "        ",
            "x",
            &format!("COEFFS64_SEG{last}"),
            plan.segments[last].c64.as_slice().len(),
            super::TrigKernel::Sin,
        )?;
        writeln!(w, "    }}")?;
        writeln!(w, "}}")?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn sin_kernel_f32(x: f32) -> f32 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = x * x; x * eval_horner_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
        }
        writeln!(
            w,
            "    else {{ let z = x * x; x * eval_horner_f32_slice(z, &COEFFS32_SEG{last}) }}"
        )?;
        writeln!(w, "}}")?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn sin_kernel_estrin_f32(x: f32) -> f32 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = x * x; x * eval_estrin_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
        }
        writeln!(
            w,
            "    else {{ let z = x * x; x * eval_estrin_f32_slice(z, &COEFFS32_SEG{last}) }}"
        )?;
        writeln!(w, "}}")?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn sin_kernel_f16(x: f16) -> f16 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= ({upper}_SEG_BOUNDS_F32[{idx}] as f16) {{ let z = x * x; x * eval_horner_f16_slice(z, &COEFFS16_SEG{idx}) }}")?;
        }
        writeln!(
            w,
            "    else {{ let z = x * x; x * eval_horner_f16_slice(z, &COEFFS16_SEG{last}) }}"
        )?;
        writeln!(w, "}}")?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn sin_kernel_f128(x: f128) -> f128 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F128[{idx}] {{ let z = x * x; x * eval_horner_f128_slice(z, &COEFFS128_SEG{idx}) }}")?;
        }
        writeln!(
            w,
            "    else {{ let z = x * x; x * eval_horner_f128_slice(z, &COEFFS128_SEG{last}) }}"
        )?;
        writeln!(w, "}}")?;
    }

    write_trig_dispatch_aliases(w)
}

pub(crate) fn write_cos_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    if plan.segments.len() == 1 {
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn cos_kernel_f64(x: f64) -> f64 {{ let z = x * x; eval_horner_f64_slice(z, &COEFFS64) }}"
		)?;
        writeln!(
            w,
            "#[inline(always)] pub(crate) fn cos_kernel_estrin_f64(x: f64) -> f64 {{"
        )?;
        super::write_straight_estrin_z_f64(
            w,
            "    ",
            "x",
            "COEFFS64",
            plan.segments[0].c64.as_slice().len(),
            super::TrigKernel::Cos,
        )?;
        writeln!(w, "}}")?;
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn cos_kernel_f32(x: f32) -> f32 {{ let z = x * x; eval_horner_f32_slice(z, &COEFFS32) }}"
		)?;
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn cos_kernel_estrin_f32(x: f32) -> f32 {{ let z = x * x; eval_estrin_f32_slice(z, &COEFFS32) }}"
		)?;
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn cos_kernel_f16(x: f16) -> f16 {{ let z = x * x; eval_horner_f16_slice(z, &COEFFS16) }}"
		)?;
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn cos_kernel_f128(x: f128) -> f128 {{ let z = x * x; eval_horner_f128_slice(z, &COEFFS128) }}"
		)?;
    } else {
        writeln!(
            w,
            "#[inline(always)] pub(crate) fn cos_kernel_f64(x: f64) -> f64 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{ let z = x * x; eval_horner_f64_slice(z, &COEFFS64_SEG{idx}) }}")?;
        }
        let last = plan.segments.len() - 1;
        writeln!(
            w,
            "    else {{ let z = x * x; eval_horner_f64_slice(z, &COEFFS64_SEG{last}) }}"
        )?;
        writeln!(w, "}}")?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn cos_kernel_estrin_f64(x: f64) -> f64 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{")?;
            super::write_straight_estrin_z_f64(
                w,
                "        ",
                "x",
                &format!("COEFFS64_SEG{idx}"),
                plan.segments[idx].c64.as_slice().len(),
                super::TrigKernel::Cos,
            )?;
            writeln!(w, "    }}")?;
        }
        writeln!(w, "    else {{")?;
        super::write_straight_estrin_z_f64(
            w,
            "        ",
            "x",
            &format!("COEFFS64_SEG{last}"),
            plan.segments[last].c64.as_slice().len(),
            super::TrigKernel::Cos,
        )?;
        writeln!(w, "    }}")?;
        writeln!(w, "}}")?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn cos_kernel_f32(x: f32) -> f32 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = x * x; eval_horner_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
        }
        writeln!(
            w,
            "    else {{ let z = x * x; eval_horner_f32_slice(z, &COEFFS32_SEG{last}) }}"
        )?;
        writeln!(w, "}}")?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn cos_kernel_estrin_f32(x: f32) -> f32 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = x * x; eval_estrin_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
        }
        writeln!(
            w,
            "    else {{ let z = x * x; eval_estrin_f32_slice(z, &COEFFS32_SEG{last}) }}"
        )?;
        writeln!(w, "}}")?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn cos_kernel_f16(x: f16) -> f16 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= ({upper}_SEG_BOUNDS_F32[{idx}] as f16) {{ let z = x * x; eval_horner_f16_slice(z, &COEFFS16_SEG{idx}) }}")?;
        }
        writeln!(
            w,
            "    else {{ let z = x * x; eval_horner_f16_slice(z, &COEFFS16_SEG{last}) }}"
        )?;
        writeln!(w, "}}")?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn cos_kernel_f128(x: f128) -> f128 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= {upper}_SEG_BOUNDS_F128[{idx}] {{ let z = x * x; eval_horner_f128_slice(z, &COEFFS128_SEG{idx}) }}")?;
        }
        writeln!(
            w,
            "    else {{ let z = x * x; eval_horner_f128_slice(z, &COEFFS128_SEG{last}) }}"
        )?;
        writeln!(w, "}}")?;
    }

    write_trig_dispatch_aliases(w)
}

pub(crate) fn write_tan_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
) -> io::Result<()> {
    if plan.segments.len() == 1 {
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn tan_kernel_f64(x: f64) -> f64 {{ let ax = x.abs(); let z = ax * ax; let t = ax + ax * z * eval_horner_f64_slice(z, &COEFFS64); if x.is_sign_negative() {{ -t }} else {{ t }} }}"
		)?;
        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_estrin_f64(x: f64) -> f64 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        writeln!(w, "    let t = {{")?;
        super::write_straight_estrin_z_f64(
            w,
            "        ",
            "ax",
            "COEFFS64",
            plan.segments[0].c64.as_slice().len(),
            super::TrigKernel::Tan,
        )?;
        writeln!(w, "    }};")?;
        writeln!(w, "    if x.is_sign_negative() {{ -t }} else {{ t }}")?;
        writeln!(w, "}}")?;
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn tan_kernel_f32(x: f32) -> f32 {{ let ax = x.abs(); let z = ax * ax; let t = ax + ax * z * eval_horner_f32_slice(z, &COEFFS32); if x.is_sign_negative() {{ -t }} else {{ t }} }}"
		)?;
        writeln!(
			w,
			"#[inline(always)] pub(crate) fn tan_kernel_estrin_f32(x: f32) -> f32 {{ let ax = x.abs(); let z = ax * ax; let t = ax + ax * z * eval_estrin_f32_slice(z, &COEFFS32); if x.is_sign_negative() {{ -t }} else {{ t }} }}"
		)?;
        writeln!(
            w,
            "#[inline(always)] fn k_tan_f64(x: f64) -> f64 {{ tan_kernel_estrin_f64(x) }}"
        )?;
        writeln!(
            w,
            "#[inline(always)] fn k_tan_f32(x: f32) -> f32 {{ tan_kernel_f32(x) }}"
        )?;
        writeln!(w, "#[inline(always)] fn k_tan_f16(x: f16) -> f16 {{ let z = x * x; x + x * z * eval_horner_f16_slice(z, &COEFFS16) }}")?;
        writeln!(w, "#[inline(always)] fn k_tan_f128(x: f128) -> f128 {{ let z = x * x; x + x * z * eval_horner_f128_slice(z, &COEFFS128) }}")?;
        writeln!(w)?;
    } else {
        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_f64(x: f64) -> f64 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        writeln!(w, "    let t = if ax <= {upper}_SEG_BOUNDS_F64[0] {{ let z = ax * ax; ax + ax * z * eval_horner_f64_slice(z, &COEFFS64_SEG0) }}")?;
        for idx in 1..plan.segments.len() {
            writeln!(w, "    else if ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{ let z = ax * ax; ax + ax * z * eval_horner_f64_slice(z, &COEFFS64_SEG{idx}) }}")?;
        }
        let last = plan.segments.len() - 1;
        writeln!(w, "    else {{ let z = ax * ax; ax + ax * z * eval_horner_f64_slice(z, &COEFFS64_SEG{last}) }};")?;
        writeln!(w, "    if x.is_sign_negative() {{ -t }} else {{ t }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_estrin_f64(x: f64) -> f64 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        writeln!(w, "    let t = if ax <= {upper}_SEG_BOUNDS_F64[0] {{")?;
        super::write_straight_estrin_z_f64(
            w,
            "        ",
            "ax",
            "COEFFS64_SEG0",
            plan.segments[0].c64.as_slice().len(),
            super::TrigKernel::Tan,
        )?;
        writeln!(w, "    }}")?;
        for idx in 1..plan.segments.len() {
            writeln!(w, "    else if ax <= {upper}_SEG_BOUNDS_F64[{idx}] {{")?;
            super::write_straight_estrin_z_f64(
                w,
                "        ",
                "ax",
                &format!("COEFFS64_SEG{idx}"),
                plan.segments[idx].c64.as_slice().len(),
                super::TrigKernel::Tan,
            )?;
            writeln!(w, "    }}")?;
        }
        writeln!(w, "    else {{")?;
        super::write_straight_estrin_z_f64(
            w,
            "        ",
            "ax",
            &format!("COEFFS64_SEG{last}"),
            plan.segments[last].c64.as_slice().len(),
            super::TrigKernel::Tan,
        )?;
        writeln!(w, "    }};")?;
        writeln!(w, "    if x.is_sign_negative() {{ -t }} else {{ t }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_f32(x: f32) -> f32 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        writeln!(w, "    let t = if ax <= {upper}_SEG_BOUNDS_F32[0] {{ let z = ax * ax; ax + ax * z * eval_horner_f32_slice(z, &COEFFS32_SEG0) }}")?;
        for idx in 1..plan.segments.len() {
            writeln!(w, "    else if ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = ax * ax; ax + ax * z * eval_horner_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
        }
        writeln!(w, "    else {{ let z = ax * ax; ax + ax * z * eval_horner_f32_slice(z, &COEFFS32_SEG{last}) }};")?;
        writeln!(w, "    if x.is_sign_negative() {{ -t }} else {{ t }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;

        writeln!(
            w,
            "#[inline(always)] pub(crate) fn tan_kernel_estrin_f32(x: f32) -> f32 {{"
        )?;
        writeln!(w, "    let ax = x.abs();")?;
        writeln!(w, "    let t = if ax <= {upper}_SEG_BOUNDS_F32[0] {{ let z = ax * ax; ax + ax * z * eval_estrin_f32_slice(z, &COEFFS32_SEG0) }}")?;
        for idx in 1..plan.segments.len() {
            writeln!(w, "    else if ax <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let z = ax * ax; ax + ax * z * eval_estrin_f32_slice(z, &COEFFS32_SEG{idx}) }}")?;
        }
        writeln!(w, "    else {{ let z = ax * ax; ax + ax * z * eval_estrin_f32_slice(z, &COEFFS32_SEG{last}) }};")?;
        writeln!(w, "    if x.is_sign_negative() {{ -t }} else {{ t }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;

        writeln!(
            w,
            "#[inline(always)] fn k_tan_f64(x: f64) -> f64 {{ tan_kernel_estrin_f64(x) }}"
        )?;
        writeln!(
            w,
            "#[inline(always)] fn k_tan_f32(x: f32) -> f32 {{ tan_kernel_f32(x) }}"
        )?;
        writeln!(w, "#[inline(always)] fn k_tan_f16(x: f16) -> f16 {{")?;
        writeln!(w, "    let ax = x.abs();")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} ax <= ({upper}_SEG_BOUNDS_F32[{idx}] as f16) {{ let z = x * x; x + x * z * eval_horner_f16_slice(z, &COEFFS16_SEG{idx}) }}")?;
        }
        writeln!(w, "    else {{ let z = x * x; x + x * z * eval_horner_f16_slice(z, &COEFFS16_SEG{last}) }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;

        writeln!(w, "#[inline(always)] fn k_tan_f128(x: f128) -> f128 {{")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            writeln!(w, "    {kw} x <= {upper}_SEG_BOUNDS_F128[{idx}] {{ let z = x * x; x + x * z * eval_horner_f128_slice(z, &COEFFS128_SEG{idx}) }}")?;
        }
        writeln!(w, "    else {{ let z = x * x; x + x * z * eval_horner_f128_slice(z, &COEFFS128_SEG{last}) }}")?;
        writeln!(w, "}}")?;
        writeln!(w)?;
    }

    writeln!(
        w,
        "#[inline(always)] fn reconstruct_tan_f64(q: u8, t: f64) -> f64 {{"
    )?;
    writeln!(w, "    match q & 1 {{")?;
    writeln!(w, "        0 => t,")?;
    writeln!(w, "        1 => -1.0_f64 / t,")?;
    writeln!(w, "        _ => unreachable!(),")?;
    writeln!(w, "    }}")?;
    writeln!(w, "}}")?;
    writeln!(w)?;
    writeln!(
        w,
        "#[inline(always)] fn reconstruct_tan_f32(q: u8, t: f32) -> f32 {{"
    )?;
    writeln!(w, "    match q & 1 {{")?;
    writeln!(w, "        0 => t,")?;
    writeln!(w, "        1 => -1.0_f32 / t,")?;
    writeln!(w, "        _ => unreachable!(),")?;
    writeln!(w, "    }}")?;
    writeln!(w, "}}")?;
    writeln!(w)?;
    writeln!(
        w,
        "#[inline(always)] fn reconstruct_tan_f16(q: u8, t: f16) -> f16 {{"
    )?;
    writeln!(w, "    let tf = t as f32;")?;
    writeln!(w, "    let y = match q & 1 {{")?;
    writeln!(w, "        0 => tf,")?;
    writeln!(w, "        1 => -1.0_f32 / tf,")?;
    writeln!(w, "        _ => unreachable!(),")?;
    writeln!(w, "    }};")?;
    writeln!(w, "    y as f16")?;
    writeln!(w, "}}")?;
    writeln!(w)?;
    writeln!(
        w,
        "#[inline(always)] fn reconstruct_tan_f128(q: u8, t: f128) -> f128 {{"
    )?;
    writeln!(w, "    match q & 1 {{")?;
    writeln!(w, "        0 => t,")?;
    writeln!(w, "        1 => -1.0_f128 / t,")?;
    writeln!(w, "        _ => unreachable!(),")?;
    writeln!(w, "    }}")?;
    writeln!(w, "}}")?;
    writeln!(w)?;
    Ok(())
}

pub(crate) fn write_simple_kernel_helpers<W: Write>(
    w: &mut W,
    plan: &super::FunctionPlan,
    upper: &str,
    is_log: bool,
    prefer_f32_for_f16: bool,
) -> io::Result<()> {
    if plan.segments.len() == 1 {
        if is_log {
            writeln!(w, "#[inline(always)] fn k_simple_f64(x: f64) -> f64 {{ let u = x - 1.0_f64; let s = u / (x + 1.0_f64); let z = s * s; let q = eval_horner_f64_slice(z, &COEFFS64); 2.0_f64 * s * q.mul_add(z, 1.0_f64) }}")?;
            writeln!(w, "#[inline(always)] fn k_simple_f32(x: f32) -> f32 {{ let u = x - 1.0_f32; let s = u / (x + 1.0_f32); let z = s * s; let q = eval_horner_f32_slice(z, &COEFFS32); 2.0_f32 * s * q.mul_add(z, 1.0_f32) }}")?;
        } else {
            writeln!(w, "#[inline(always)] fn k_simple_f64(x: f64) -> f64 {{ eval_horner_f64_slice(x, &COEFFS64) }}")?;
            writeln!(w, "#[inline(always)] fn k_simple_f32(x: f32) -> f32 {{ eval_horner_f32_slice(x, &COEFFS32) }}")?;
        }
        if prefer_f32_for_f16 {
            if is_log {
                writeln!(w, "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{ let xf = x as f32; let u = xf - 1.0_f32; let s = u / (xf + 1.0_f32); let z = s * s; let q = eval_horner_f32_slice(z, &COEFFS32); (2.0_f32 * s * q.mul_add(z, 1.0_f32)) as f16 }}")?;
            } else {
                writeln!(w, "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{ k_simple_f32(x as f32) as f16 }}")?;
            }
        } else if is_log {
            writeln!(w, "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{ let u = x - 1.0_f16; let s = u / (x + 1.0_f16); let z = s * s; let q = eval_horner_f16_slice(z, &COEFFS16); 2.0_f16 * s * (1.0_f16 + z * q) }}")?;
        } else {
            writeln!(w, "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{ eval_horner_f16_slice(x, &COEFFS16) }}")?;
        }
        if is_log {
            writeln!(w, "#[inline(always)] fn k_simple_f128(x: f128) -> f128 {{ let u = x - 1.0_f128; let s = u / (x + 1.0_f128); let z = s * s; let q = eval_horner_f128_slice(z, &COEFFS128); 2.0_f128 * s * (1.0_f128 + z * q) }}")?;
        } else {
            writeln!(w, "#[inline(always)] fn k_simple_f128(x: f128) -> f128 {{ eval_horner_f128_slice(x, &COEFFS128) }}")?;
        }
        writeln!(w)?;
    } else {
        writeln!(w, "#[inline(always)] fn k_simple_f64(x: f64) -> f64 {{")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            if is_log {
                writeln!(w, "    {kw} x <= {upper}_SEG_BOUNDS_F64[{idx}] {{ let u = x - 1.0_f64; let s = u / (x + 1.0_f64); let z = s * s; let q = eval_horner_f64_slice(z, &COEFFS64_SEG{idx}); 2.0_f64 * s * q.mul_add(z, 1.0_f64) }}")?;
            } else {
                writeln!(w, "    {kw} x <= {upper}_SEG_BOUNDS_F64[{idx}] {{ eval_horner_f64_slice(x, &COEFFS64_SEG{idx}) }}")?;
            }
        }
        let last = plan.segments.len() - 1;
        if is_log {
            writeln!(w, "    else {{ let u = x - 1.0_f64; let s = u / (x + 1.0_f64); let z = s * s; let q = eval_horner_f64_slice(z, &COEFFS64_SEG{last}); 2.0_f64 * s * q.mul_add(z, 1.0_f64) }}")?;
        } else {
            writeln!(
                w,
                "    else {{ eval_horner_f64_slice(x, &COEFFS64_SEG{last}) }}"
            )?;
        }
        writeln!(w, "}}")?;

        writeln!(w, "#[inline(always)] fn k_simple_f32(x: f32) -> f32 {{")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            if is_log {
                writeln!(w, "    {kw} x <= {upper}_SEG_BOUNDS_F32[{idx}] {{ let u = x - 1.0_f32; let s = u / (x + 1.0_f32); let z = s * s; let q = eval_horner_f32_slice(z, &COEFFS32_SEG{idx}); 2.0_f32 * s * q.mul_add(z, 1.0_f32) }}")?;
            } else {
                writeln!(w, "    {kw} x <= {upper}_SEG_BOUNDS_F32[{idx}] {{ eval_horner_f32_slice(x, &COEFFS32_SEG{idx}) }}")?;
            }
        }
        if is_log {
            writeln!(w, "    else {{ let u = x - 1.0_f32; let s = u / (x + 1.0_f32); let z = s * s; let q = eval_horner_f32_slice(z, &COEFFS32_SEG{last}); 2.0_f32 * s * q.mul_add(z, 1.0_f32) }}")?;
        } else {
            writeln!(
                w,
                "    else {{ eval_horner_f32_slice(x, &COEFFS32_SEG{last}) }}"
            )?;
        }
        writeln!(w, "}}")?;

        if prefer_f32_for_f16 {
            writeln!(w, "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{ k_simple_f32(x as f32) as f16 }}")?;
        } else {
            writeln!(w, "#[inline(always)] fn k_simple_f16(x: f16) -> f16 {{")?;
            for idx in 0..plan.segments.len() {
                let kw = if idx == 0 { "if" } else { "else if" };
                if is_log {
                    writeln!(w, "    {kw} x <= ({upper}_SEG_BOUNDS_F32[{idx}] as f16) {{ let u = x - 1.0_f16; let s = u / (x + 1.0_f16); let z = s * s; let q = eval_horner_f16_slice(z, &COEFFS16_SEG{idx}); 2.0_f16 * s * (1.0_f16 + z * q) }}")?;
                } else {
                    writeln!(w, "    {kw} x <= ({upper}_SEG_BOUNDS_F32[{idx}] as f16) {{ eval_horner_f16_slice(x, &COEFFS16_SEG{idx}) }}")?;
                }
            }
            if is_log {
                writeln!(w, "    else {{ let u = x - 1.0_f16; let s = u / (x + 1.0_f16); let z = s * s; let q = eval_horner_f16_slice(z, &COEFFS16_SEG{last}); 2.0_f16 * s * (1.0_f16 + z * q) }}")?;
            } else {
                writeln!(
                    w,
                    "    else {{ eval_horner_f16_slice(x, &COEFFS16_SEG{last}) }}"
                )?;
            }
            writeln!(w, "}}")?;
        }

        writeln!(w, "#[inline(always)] fn k_simple_f128(x: f128) -> f128 {{")?;
        for idx in 0..plan.segments.len() {
            let kw = if idx == 0 { "if" } else { "else if" };
            if is_log {
                writeln!(w, "    {kw} x <= {upper}_SEG_BOUNDS_F128[{idx}] {{ let u = x - 1.0_f128; let s = u / (x + 1.0_f128); let z = s * s; let q = eval_horner_f128_slice(z, &COEFFS128_SEG{idx}); 2.0_f128 * s * (1.0_f128 + z * q) }}")?;
            } else {
                writeln!(w, "    {kw} x <= {upper}_SEG_BOUNDS_F128[{idx}] {{ eval_horner_f128_slice(x, &COEFFS128_SEG{idx}) }}")?;
            }
        }
        if is_log {
            writeln!(w, "    else {{ let u = x - 1.0_f128; let s = u / (x + 1.0_f128); let z = s * s; let q = eval_horner_f128_slice(z, &COEFFS128_SEG{last}); 2.0_f128 * s * (1.0_f128 + z * q) }}")?;
        } else {
            writeln!(
                w,
                "    else {{ eval_horner_f128_slice(x, &COEFFS128_SEG{last}) }}"
            )?;
        }
        writeln!(w, "}}")?;
        writeln!(w)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{write_scalar_public_body, PublicFnCtx, ScalarKind};

    #[test]
    fn generated_public_scalars_do_not_delegate_to_platform_transcendentals() {
        const FUNCTIONS: [&str; 13] = [
            "sin", "cos", "tan", "exp", "expm1", "tanh", "sinh", "cosh", "log", "log1p", "asin",
            "acos", "atan",
        ];
        const FORBIDDEN: [&str; 13] = [
            ".sin()",
            ".cos()",
            ".tan()",
            ".exp()",
            ".exp_m1()",
            ".tanh()",
            ".sinh()",
            ".cosh()",
            ".ln()",
            ".ln_1p()",
            ".asin()",
            ".acos()",
            ".atan()",
        ];

        for function in FUNCTIONS {
            let ident = function.to_owned();
            let upper = function.to_ascii_uppercase();
            let ctx = PublicFnCtx {
                ident: &ident,
                upper: &upper,
            };

            for kind in [
                ScalarKind::F16,
                ScalarKind::F32,
                ScalarKind::F64,
                ScalarKind::F128,
            ] {
                let mut body = Vec::new();
                write_scalar_public_body(&mut body, function, &ctx, kind).unwrap();
                let body = String::from_utf8(body).unwrap();

                for forbidden in FORBIDDEN {
                    assert!(
                        !body.contains(forbidden),
                        "{function} emitted platform delegation `{forbidden}`:\n{body}"
                    );
                }
            }
        }
    }
}
