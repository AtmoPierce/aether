// aether-fpx-gen/src/design/remez/create/coefficients.rs

use core::{f128, f16};

/* ------------------------------------------------------------------------------------------------
 * Polynomial domain & kernel
 * ------------------------------------------------------------------------------------------------ */

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PolyDomain {
    /// Polynomial in x
    X,
    /// Polynomial in z = x²
    Z,
}

impl PolyDomain {
    #[inline(always)]
    pub fn from_function_name(name: &str) -> Self {
        match name {
            // z = x² domain
            "sin" | "cos" | "tan" => PolyDomain::Z,

            // direct x-domain polynomial
            "sinh" | "cosh" | "asin" | "acos" | "atan" | "tanh" | "log" | "log1p" | "exp"
            | "expm1" => PolyDomain::X,

            _ => panic!("Unknown polynomial domain for function: {name}"),
        }
    }

    #[inline(always)]
    pub fn table_domain_str(self) -> &'static str {
        match self {
            PolyDomain::X => "x",
            PolyDomain::Z => "z",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrigKernel {
    Sin,
    Cos,
    Tan,
}

/// Per-segment runtime shape.
/// - X: evaluate P(x)
/// - Z(k): evaluate P(z=x^2) and reconstruct via kernel:
///     Cos: P(z)
///     Sin: x * P(z)
///     Tan: x + x^3 * P(z)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SegmentShape {
    X,
    Z(TrigKernel),
}

/// Canonical domain selection aligned with supported functions.
pub fn poly_domain(name: &str) -> PolyDomain {
    PolyDomain::from_function_name(name)
}

/// Kernel selection for Z-domain functions.
/// Panics if called for a non-(sin|cos) function name.
pub fn trig_kernel(name: &str) -> TrigKernel {
    match name {
        "sin" => TrigKernel::Sin,
        "cos" => TrigKernel::Cos,
        "tan" => TrigKernel::Tan,
        _ => panic!("Unknown trig kernel: {}", name),
    }
}

/* ------------------------------------------------------------------------------------------------
 * Polynomial shape
 * ------------------------------------------------------------------------------------------------ */

#[derive(Copy, Clone, Debug)]
pub struct Nk {
    pub n: usize,
    pub k: usize,
}

/* ------------------------------------------------------------------------------------------------
 * Coefficient sets
 * ------------------------------------------------------------------------------------------------ */

#[derive(Copy, Clone, Debug)]
pub enum CoeffSet64 {
    N2([f64; 2]),
    N3([f64; 3]),
    N4([f64; 4]),
    N5([f64; 5]),
    N6([f64; 6]),
    N7([f64; 7]),
    N8([f64; 8]),
    N9([f64; 9]),
    N10([f64; 10]),
    N11([f64; 11]),
    N12([f64; 12]),
    N13([f64; 13]),
    N14([f64; 14]),
    N15([f64; 15]),
    N16([f64; 16]),
    N17([f64; 17]),
    N18([f64; 18]),
    N19([f64; 19]),
    N20([f64; 20]),
    N21([f64; 21]),
    N22([f64; 22]),
    N23([f64; 23]),
    N24([f64; 24]),
    N25([f64; 25]),
}

#[derive(Copy, Clone, Debug)]
pub enum CoeffSet128 {
    N2([f128; 2]),
    N3([f128; 3]),
    N4([f128; 4]),
    N5([f128; 5]),
    N6([f128; 6]),
    N7([f128; 7]),
    N8([f128; 8]),
    N9([f128; 9]),
    N10([f128; 10]),
    N11([f128; 11]),
    N12([f128; 12]),
    N13([f128; 13]),
    N14([f128; 14]),
    N15([f128; 15]),
    N16([f128; 16]),
    N17([f128; 17]),
    N18([f128; 18]),
    N19([f128; 19]),
    N20([f128; 20]),
    N21([f128; 21]),
    N22([f128; 22]),
    N23([f128; 23]),
    N24([f128; 24]),
    N25([f128; 25]),
}

#[derive(Copy, Clone, Debug)]
pub enum CoeffSet32 {
    N2([f32; 2]),
    N3([f32; 3]),
    N4([f32; 4]),
    N5([f32; 5]),
    N6([f32; 6]),
    N7([f32; 7]),
    N8([f32; 8]),
    N9([f32; 9]),
    N10([f32; 10]),
    N11([f32; 11]),
    N12([f32; 12]),
    N13([f32; 13]),
    N14([f32; 14]),
    N15([f32; 15]),
    N16([f32; 16]),
    N17([f32; 17]),
    N18([f32; 18]),
    N19([f32; 19]),
    N20([f32; 20]),
    N21([f32; 21]),
    N22([f32; 22]),
    N23([f32; 23]),
    N24([f32; 24]),
    N25([f32; 25]),
}

#[derive(Copy, Clone, Debug)]
pub enum CoeffSet16 {
    N2([f16; 2]),
    N3([f16; 3]),
    N4([f16; 4]),
    N5([f16; 5]),
    N6([f16; 6]),
    N7([f16; 7]),
    N8([f16; 8]),
    N9([f16; 9]),
    N10([f16; 10]),
    N11([f16; 11]),
    N12([f16; 12]),
    N13([f16; 13]),
    N14([f16; 14]),
    N15([f16; 15]),
    N16([f16; 16]),
    N17([f16; 17]),
    N18([f16; 18]),
    N19([f16; 19]),
    N20([f16; 20]),
    N21([f16; 21]),
    N22([f16; 22]),
    N23([f16; 23]),
    N24([f16; 24]),
    N25([f16; 25]),
}

macro_rules! impl_coeff_set {
    ($name:ident, $ty:ty) => {
        impl $name {
            #[inline(always)]
            pub fn as_slice(&self) -> &[$ty] {
                match self {
                    Self::N2(cs) => cs,
                    Self::N3(cs) => cs,
                    Self::N4(cs) => cs,
                    Self::N5(cs) => cs,
                    Self::N6(cs) => cs,
                    Self::N7(cs) => cs,
                    Self::N8(cs) => cs,
                    Self::N9(cs) => cs,
                    Self::N10(cs) => cs,
                    Self::N11(cs) => cs,
                    Self::N12(cs) => cs,
                    Self::N13(cs) => cs,
                    Self::N14(cs) => cs,
                    Self::N15(cs) => cs,
                    Self::N16(cs) => cs,
                    Self::N17(cs) => cs,
                    Self::N18(cs) => cs,
                    Self::N19(cs) => cs,
                    Self::N20(cs) => cs,
                    Self::N21(cs) => cs,
                    Self::N22(cs) => cs,
                    Self::N23(cs) => cs,
                    Self::N24(cs) => cs,
                    Self::N25(cs) => cs,
                }
            }

            pub fn from_slice(cs: &[$ty]) -> Self {
                match cs.len() {
                    2 => Self::N2(cs.try_into().unwrap()),
                    3 => Self::N3(cs.try_into().unwrap()),
                    4 => Self::N4(cs.try_into().unwrap()),
                    5 => Self::N5(cs.try_into().unwrap()),
                    6 => Self::N6(cs.try_into().unwrap()),
                    7 => Self::N7(cs.try_into().unwrap()),
                    8 => Self::N8(cs.try_into().unwrap()),
                    9 => Self::N9(cs.try_into().unwrap()),
                    10 => Self::N10(cs.try_into().unwrap()),
                    11 => Self::N11(cs.try_into().unwrap()),
                    12 => Self::N12(cs.try_into().unwrap()),
                    13 => Self::N13(cs.try_into().unwrap()),
                    14 => Self::N14(cs.try_into().unwrap()),
                    15 => Self::N15(cs.try_into().unwrap()),
                    16 => Self::N16(cs.try_into().unwrap()),
                    17 => Self::N17(cs.try_into().unwrap()),
                    18 => Self::N18(cs.try_into().unwrap()),
                    19 => Self::N19(cs.try_into().unwrap()),
                    20 => Self::N20(cs.try_into().unwrap()),
                    21 => Self::N21(cs.try_into().unwrap()),
                    22 => Self::N22(cs.try_into().unwrap()),
                    23 => Self::N23(cs.try_into().unwrap()),
                    24 => Self::N24(cs.try_into().unwrap()),
                    25 => Self::N25(cs.try_into().unwrap()),
                    _ => panic!(
                        "unsupported coefficient length {}; supported range is N2..N25",
                        cs.len()
                    ),
                }
            }
        }
    };
}

impl_coeff_set!(CoeffSet64, f64);
impl_coeff_set!(CoeffSet128, f128);
impl_coeff_set!(CoeffSet32, f32);
impl_coeff_set!(CoeffSet16, f16);

/* ------------------------------------------------------------------------------------------------
 * Plans
 * ------------------------------------------------------------------------------------------------ */

#[derive(Clone, Debug)]
pub struct SegmentPlan {
    /// Lower bound of the segment
    pub a128: f128,
    pub a64: f64,
    pub a32: f32,

    /// Upper bound of the segment
    pub b128: f128,
    pub b64: f64,
    pub b32: f32,

    /// Polynomial coefficients in solve-domain
    pub c128: CoeffSet128,
    pub c64: CoeffSet64,
    pub c32: CoeffSet32,
    pub c16: CoeffSet16,

    pub nk: Nk,

    /// Runtime evaluation shape (X vs Z + kernel)
    pub shape: SegmentShape,
}

#[derive(Clone, Debug)]
pub struct FunctionPlan {
    pub name: String,

    /// Primary domain intent (informational; runtime uses per-segment `shape`)
    pub dom: PolyDomain,

    /// Function-wide kernel (informational; runtime uses per-segment `shape`)
    pub kernel: Option<TrigKernel>,

    pub segments: Vec<SegmentPlan>, // sorted by b64/b32

    pub interval_min_f128: f128,
    pub interval_max_f128: f128,
    pub interval_min_f64: f64,
    pub interval_max_f64: f64,
    pub interval_min_f32: f32,
    pub interval_max_f32: f32,
}

/* ------------------------------------------------------------------------------------------------
 * Horner evaluators
 * ------------------------------------------------------------------------------------------------ */

#[inline(always)]
fn eval_horner_f64_slice(x: f64, cs: &[f64]) -> f64 {
    let mut p: f64 = 0.0;
    for &c in cs.iter().rev() {
        // Keep planning and acceptance bit-for-bit aligned with
        // aether-fpx::trig::common::eval_horner_f64_slice.
        p = p.mul_add(x, c);
    }
    p
}

#[inline(always)]
fn eval_horner_f32_slice(x: f32, cs: &[f32]) -> f32 {
    let mut p: f32 = 0.0;
    for &c in cs.iter().rev() {
        p = p.mul_add(x, c);
    }
    p
}

#[inline(always)]
fn eval_horner_f128_slice(x: f128, cs: &[f128]) -> f128 {
    let mut p: f128 = 0.0_f128;
    for &c in cs.iter().rev() {
        p = p * x + c;
    }
    p
}

/* ------------------------------------------------------------------------------------------------
 * Domain-specialized evaluation
 * ------------------------------------------------------------------------------------------------ */

#[inline(always)]
pub fn eval_coeff64_x(x: f64, c: &CoeffSet64) -> f64 {
    eval_horner_f64_slice(x, c.as_slice())
}

#[inline(always)]
pub fn eval_coeff32_x(x: f32, c: &CoeffSet32) -> f32 {
    eval_horner_f32_slice(x, c.as_slice())
}

#[inline(always)]
pub fn eval_coeff128_x(x: f128, c: &CoeffSet128) -> f128 {
    eval_horner_f128_slice(x, c.as_slice())
}

#[inline(always)]
pub fn eval_coeff64_log_atanh(x: f64, c: &CoeffSet64) -> f64 {
    let u = x - 1.0_f64;
    let s = u / (x + 1.0_f64);
    let z = s * s;
    let q = eval_coeff64_x(z, c);
    2.0_f64 * s * q.mul_add(z, 1.0_f64)
}

#[inline(always)]
pub fn eval_coeff32_log_atanh(x: f32, c: &CoeffSet32) -> f32 {
    let u = x - 1.0_f32;
    let s = u / (x + 1.0_f32);
    let z = s * s;
    let q = eval_coeff32_x(z, c);
    2.0_f32 * s * q.mul_add(z, 1.0_f32)
}

#[inline(always)]
pub fn eval_coeff128_log_atanh(x: f128, c: &CoeffSet128) -> f128 {
    let u = x - 1.0_f128;
    let s = u / (x + 1.0_f128);
    let z = s * s;
    let q = eval_coeff128_x(z, c);
    2.0_f128 * s * (1.0_f128 + z * q)
}

#[inline(always)]
pub fn eval_coeff64_z(x: f64, c: &CoeffSet64, k: TrigKernel) -> f64 {
    let z = x * x;

    match k {
        TrigKernel::Cos => eval_coeff64_x(z, c),
        TrigKernel::Sin => x * eval_coeff64_x(z, c),
        TrigKernel::Tan => x + x * z * eval_coeff64_x(z, c),
    }
}

#[inline(always)]
pub fn eval_coeff32_z(x: f32, c: &CoeffSet32, k: TrigKernel) -> f32 {
    let z = x * x;

    match k {
        TrigKernel::Cos => eval_coeff32_x(z, c),
        TrigKernel::Sin => x * eval_coeff32_x(z, c),
        TrigKernel::Tan => x + x * z * eval_coeff32_x(z, c),
    }
}

#[inline(always)]
pub fn eval_coeff128_z(x: f128, c: &CoeffSet128, k: TrigKernel) -> f128 {
    let z = x * x;

    match k {
        TrigKernel::Cos => eval_coeff128_x(z, c),
        TrigKernel::Sin => x * eval_coeff128_x(z, c),
        TrigKernel::Tan => x + x * z * eval_coeff128_x(z, c),
    }
}

#[inline(always)]
pub fn eval_coeff64_shape(x: f64, c: &CoeffSet64, shape: SegmentShape) -> f64 {
    match shape {
        SegmentShape::X => eval_coeff64_x(x, c),
        SegmentShape::Z(k) => eval_coeff64_z(x, c, k),
    }
}

#[inline(always)]
pub fn eval_coeff32_shape(x: f32, c: &CoeffSet32, shape: SegmentShape) -> f32 {
    match shape {
        SegmentShape::X => eval_coeff32_x(x, c),
        SegmentShape::Z(k) => eval_coeff32_z(x, c, k),
    }
}

#[inline(always)]
pub fn eval_coeff128_shape(x: f128, c: &CoeffSet128, shape: SegmentShape) -> f128 {
    match shape {
        SegmentShape::X => eval_coeff128_x(x, c),
        SegmentShape::Z(k) => eval_coeff128_z(x, c, k),
    }
}

/* ------------------------------------------------------------------------------------------------
 * Segment lookup (single-bound binary search)
 * ------------------------------------------------------------------------------------------------ */

#[inline(always)]
fn find_segment_f64(segs: &[SegmentPlan], x: f64) -> usize {
    debug_assert!(!segs.is_empty());

    let mut lo = 0usize;
    let mut hi = segs.len();

    while lo < hi {
        let mid = (lo + hi) >> 1;

        if x <= segs[mid].b64 {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }

    lo.min(segs.len() - 1)
}

#[inline(always)]
fn find_segment_f32(segs: &[SegmentPlan], x: f32) -> usize {
    debug_assert!(!segs.is_empty());

    let mut lo = 0usize;
    let mut hi = segs.len();

    while lo < hi {
        let mid = (lo + hi) >> 1;

        if x <= segs[mid].b32 {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }

    lo.min(segs.len() - 1)
}

/* ------------------------------------------------------------------------------------------------
 * Public evaluation entry points
 * ------------------------------------------------------------------------------------------------ */

#[inline(always)]
pub fn eval_plan_f64(plan: &FunctionPlan, x: f64) -> f64 {
    let i = find_segment_f64(&plan.segments, x);
    let seg = unsafe { plan.segments.get_unchecked(i) };

    if plan.name == "log" {
        eval_coeff64_log_atanh(x, &seg.c64)
    } else {
        eval_coeff64_shape(x, &seg.c64, seg.shape)
    }
}

#[inline(always)]
pub fn eval_plan_f32(plan: &FunctionPlan, x: f32) -> f32 {
    let i = find_segment_f32(&plan.segments, x);
    let seg = unsafe { plan.segments.get_unchecked(i) };

    if plan.name == "log" {
        eval_coeff32_log_atanh(x, &seg.c32)
    } else {
        eval_coeff32_shape(x, &seg.c32, seg.shape)
    }
}

#[cfg(test)]
mod tests {
    use super::{eval_coeff32_x, eval_coeff64_x, CoeffSet32, CoeffSet64};

    #[test]
    fn planner_f64_horner_matches_runtime_fma_semantics() {
        // (1 + d) * (1 - d) - 1 is -d^2 exactly.  The product rounds to 1
        // before the subtraction when multiply and add are evaluated
        // separately, making this a regression case for accidental non-FMA
        // planner evaluation.
        let d = 2.0_f64.powi(-28);
        let leading = 1.0 + d;
        let x = 1.0 - d;
        let coefficients = CoeffSet64::N2([-1.0, leading]);

        let planned = eval_coeff64_x(x, &coefficients);
        let runtime_fma = leading.mul_add(x, -1.0);
        let unfused = leading * x - 1.0;

        assert_eq!(planned.to_bits(), runtime_fma.to_bits());
        assert_eq!(runtime_fma.to_bits(), (-(d * d)).to_bits());
        assert_ne!(planned.to_bits(), unfused.to_bits());
    }

    #[test]
    fn planner_f32_horner_matches_runtime_fma_semantics() {
        let d = 2.0_f32.powi(-13);
        let leading = 1.0 + d;
        let x = 1.0 - d;
        let coefficients = CoeffSet32::N2([-1.0, leading]);

        let planned = eval_coeff32_x(x, &coefficients);
        let runtime_fma = leading.mul_add(x, -1.0);
        let unfused = leading * x - 1.0;

        assert_eq!(planned.to_bits(), runtime_fma.to_bits());
        assert_eq!(runtime_fma.to_bits(), (-(d * d)).to_bits());
        assert_ne!(planned.to_bits(), unfused.to_bits());
    }
}
