use core::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Copy)]
struct FormatSpec {
    name: &'static str,
    exp_bits: u8,
    man_bits: u8,
    bias: i32,
    has_inf: bool,
    nan_exp: u8,
    nan_man: u8,
    ei_low_bits: u8,
    mi_high_bits: u8,
}

const E5M2: FormatSpec = FormatSpec {
    name: "E5M2",
    exp_bits: 5,
    man_bits: 2,
    bias: 15,
    has_inf: true,
    nan_exp: 0x1F,
    nan_man: 0,
    ei_low_bits: 2,
    mi_high_bits: 1,
};

const E4M3: FormatSpec = FormatSpec {
    name: "E4M3",
    exp_bits: 4,
    man_bits: 3,
    bias: 7,
    has_inf: false,
    nan_exp: 0x0F,
    nan_man: 0x07,
    ei_low_bits: 2,
    mi_high_bits: 2,
};

fn decode_fp8(bits: u8, spec: FormatSpec) -> f32 {
    let sign = if ((bits >> 7) & 1) == 1 { -1.0 } else { 1.0 };
    let exp_mask = (1u8 << spec.exp_bits) - 1;
    let man_mask = (1u8 << spec.man_bits) - 1;

    let exp = (bits >> spec.man_bits) & exp_mask;
    let man = bits & man_mask;

    if spec.name == "E5M2" && exp == 0x1F {
        if man == 0 {
            return sign * f32::INFINITY;
        }
        return f32::NAN;
    }

    if spec.name == "E4M3" && exp == spec.nan_exp && man == spec.nan_man {
        return f32::NAN;
    }

    if exp == 0 {
        if man == 0 {
            return sign * 0.0;
        }
        let frac = (man as f32) / ((1u32 << spec.man_bits) as f32);
        return sign * frac * 2.0f32.powi(1 - spec.bias);
    }

    let frac = 1.0 + ((man as f32) / ((1u32 << spec.man_bits) as f32));
    sign * frac * 2.0f32.powi(exp as i32 - spec.bias)
}

fn encode_fp8(v: f32, spec: FormatSpec) -> u8 {
    if v.is_nan() {
        return 0x7F;
    }

    let mut x = v;
    if x.is_infinite() {
        if spec.has_inf {
            return if x.is_sign_positive() { 0x7C } else { 0xFC };
        }
        x = if x.is_sign_positive() {
            f32::MAX
        } else {
            -f32::MAX
        };
    }

    let mut best_bits = 0u8;
    let mut best_err = f32::INFINITY;

    for bits in 0u8..=255 {
        let decoded = decode_fp8(bits, spec);
        if decoded.is_nan() || decoded.is_infinite() {
            continue;
        }
        let err = (decoded - x).abs();
        if err < best_err {
            best_err = err;
            best_bits = bits;
        } else if err == best_err {
            let cur_lsb = best_bits & 1;
            let new_lsb = bits & 1;
            if new_lsb == 0 && cur_lsb == 1 {
                best_bits = bits;
            }
        }
    }

    best_bits
}

fn is_finite_encoding(bits: u8, spec: FormatSpec) -> bool {
    let v = decode_fp8(bits, spec);
    !v.is_nan() && !v.is_infinite()
}

fn k_index(exp: u8, man: u8, spec: FormatSpec) -> usize {
    let e_mask = (1u8 << spec.ei_low_bits) - 1;
    let e_low = exp & e_mask;
    let m_hi = man >> (spec.man_bits - spec.mi_high_bits);
    ((e_low as usize) << spec.mi_high_bits) | (m_hi as usize)
}

fn fit_tables(spec: FormatSpec, t1: f32, t2: f32) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let table_len = 1usize << (spec.ei_low_bits + spec.mi_high_bits);
    let mut buckets: Vec<Vec<(u8, u8)>> = vec![Vec::new(); table_len];

    let exp_mask = (1u8 << spec.exp_bits) - 1;
    let man_mask = (1u8 << spec.man_bits) - 1;

    for bits in 0u8..=255 {
        if !is_finite_encoding(bits, spec) {
            continue;
        }

        let x = decode_fp8(bits, spec);
        let ax = x.abs();
        if !(ax >= t1 && ax <= t2) {
            continue;
        }

        let exp = (bits >> spec.man_bits) & exp_mask;
        let man = bits & man_mask;
        let idx = k_index(exp, man, spec);

        let target_bits = encode_fp8(ax.tanh(), spec);
        buckets[idx].push((man, target_bits));
    }

    let max_man = (1u8 << spec.man_bits) - 1;
    let max_exp = (1u8 << spec.exp_bits) - 1;

    let mut table_e = vec![0u8; table_len];
    let mut table_r = vec![0u8; table_len];
    let mut table_b = vec![0u8; table_len];

    for idx in 0..table_len {
        if buckets[idx].is_empty() {
            table_e[idx] = spec.bias as u8;
            table_r[idx] = 0;
            table_b[idx] = 0;
            continue;
        }

        let mut best_score = f32::INFINITY;
        let mut best_et = spec.bias as u8;
        let mut best_rt = 0u8;
        let mut best_bt = 0u8;

        for et in 0u8..=max_exp {
            if spec.has_inf && et == max_exp {
                continue;
            }
            for rt in 0u8..=spec.man_bits {
                for bt in 0u8..=max_man {
                    let mut score = 0.0f32;
                    let mut valid = true;

                    for (man, target_bits) in &buckets[idx] {
                        let mut mo = (*man >> rt).saturating_add(bt);
                        if mo > max_man {
                            mo = max_man;
                        }

                        if spec.name == "E4M3" && et == spec.nan_exp && mo == spec.nan_man {
                            valid = false;
                            break;
                        }

                        let out_bits = (et << spec.man_bits) | mo;
                        let out = decode_fp8(out_bits, spec);
                        let tgt = decode_fp8(*target_bits & 0x7F, spec);
                        score += (out - tgt).abs();
                    }

                    if !valid {
                        continue;
                    }

                    if score < best_score
                        || (score == best_score && (et, rt, bt) < (best_et, best_rt, best_bt))
                    {
                        best_score = score;
                        best_et = et;
                        best_rt = rt;
                        best_bt = bt;
                    }
                }
            }
        }

        table_e[idx] = best_et;
        table_r[idx] = best_rt;
        table_b[idx] = best_bt;
    }

    (table_e, table_r, table_b)
}

fn vec_literal(values: &[u8]) -> String {
    let mut s = String::new();
    for (i, value) in values.iter().enumerate() {
        if i > 0 {
            s.push_str(", ");
        }
        let _ = write!(s, "{}", value);
    }
    s
}

fn render_rust(
    t1: f32,
    t2: f32,
    e5: (Vec<u8>, Vec<u8>, Vec<u8>),
    e4: (Vec<u8>, Vec<u8>, Vec<u8>),
) -> String {
    let (e5_e, e5_r, e5_b) = e5;
    let (e4_e, e4_r, e4_b) = e4;

    format!(
        "// @generated by fpx tanh table generator
pub const E5M2_T1: f32 = {t1};
pub const E5M2_T2: f32 = {t2};
pub const E5M2_MI_HIGH_BITS: usize = 1;
pub const E5M2_TABLE_E: [u8; 8] = [{e5_e_lit}];
pub const E5M2_TABLE_R: [u8; 8] = [{e5_r_lit}];
pub const E5M2_TABLE_B: [u8; 8] = [{e5_b_lit}];
pub const E4M3_T1: f32 = {t1};
pub const E4M3_T2: f32 = {t2};
pub const E4M3_MI_HIGH_BITS: usize = 2;
pub const E4M3_TABLE_E: [u8; 16] = [{e4_e_lit}];
pub const E4M3_TABLE_R: [u8; 16] = [{e4_r_lit}];
pub const E4M3_TABLE_B: [u8; 16] = [{e4_b_lit}];
// @end generatred by fpx tanh table generator",
        t1 = t1,
        t2 = t2,
        e5_e_lit = vec_literal(&e5_e),
        e5_r_lit = vec_literal(&e5_r),
        e5_b_lit = vec_literal(&e5_b),
        e4_e_lit = vec_literal(&e4_e),
        e4_r_lit = vec_literal(&e4_r),
        e4_b_lit = vec_literal(&e4_b),
    )
}

fn parse_arg_f32(args: &[String], name: &str, default: f32) -> f32 {
    let key = format!("--{}", name);
    for i in 0..args.len() {
        if args[i] == key {
            if i + 1 < args.len() {
                if let Ok(v) = args[i + 1].parse::<f32>() {
                    return v;
                }
            }
        }
    }
    default
}

fn parse_arg_path(args: &[String], name: &str, default: &str) -> PathBuf {
    let key = format!("--{}", name);
    for i in 0..args.len() {
        if args[i] == key {
            if i + 1 < args.len() {
                return PathBuf::from(&args[i + 1]);
            }
        }
    }
    PathBuf::from(default)
}

fn resolve_output_path(path: &Path) -> PathBuf {
    if path.is_absolute() {
        return path.to_path_buf();
    }

    let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(4)
        .unwrap()
        .to_path_buf();
    repo_root.join(path)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let t1 = parse_arg_f32(&args, "t1", 0.5);
    let t2 = parse_arg_f32(&args, "t2", 2.0);
    let out_arg = parse_arg_path(&args, "output", "crates/aether-fpx/src/ocp/trig/tanh.rs");

    let e5 = fit_tables(E5M2, t1, t2);
    let e4 = fit_tables(E4M3, t1, t2);

    let output = resolve_output_path(&out_arg);
    if let Some(parent) = output.parent() {
        let _ = fs::create_dir_all(parent);
    }
    let generated = render_rust(t1, t2, e5, e4);
    let existing = fs::read_to_string(&output).expect("failed to read target file");

    let start = "// @generated by fpx tanh table generator";
    let end = "// @end generated by fpx tanh table generator";
    let start_idx = existing
        .find(start)
        .expect("missing BEGIN GENERATED K-TANH TABLES marker");
    let end_idx = existing
        .find(end)
        .expect("missing END GENERATED K-TANH TABLES marker");
    let end_tail = end_idx + end.len();

    let mut merged = String::new();
    merged.push_str(&existing[..start_idx]);
    merged.push_str(&generated);
    merged.push_str(&existing[end_tail..]);

    fs::write(&output, merged).expect("failed to write updated tanh file");

    println!("wrote {}", output.display());
}
