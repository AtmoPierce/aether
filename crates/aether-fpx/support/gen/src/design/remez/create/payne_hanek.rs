use num_bigint::{BigInt, BigUint};
use std::fmt::Write;

const TWO_OVER_PI_BITS: usize = 16_640;
const TWO_OVER_PI_LIMBS: usize = TWO_OVER_PI_BITS / 64;
const CHUDNOVSKY_C3_OVER_24: u64 = 10_939_058_860_032_000;
const CHUDNOVSKY_BITS_PER_TERM: usize = 46;
const LOW_GUARD_BITS: usize = 256;
const HIGH_GUARD_BITS: usize = 512;

const CONSTANTS_START: &str = "const TWO_OVER_PI_BITS:";
const CONSTANTS_END: &str = "\n\nconst PRODUCT_LIMBS:";

struct ChudnovskySplit {
    p: BigInt,
    q: BigInt,
    t: BigInt,
}

fn chudnovsky_split(a: usize, b: usize) -> ChudnovskySplit {
    if b - a == 1 {
        let (p, q) = if a == 0 {
            (BigInt::from(1_u8), BigInt::from(1_u8))
        } else {
            let a = a as u64;
            let p = BigInt::from(6 * a - 5) * BigInt::from(2 * a - 1) * BigInt::from(6 * a - 1);
            let q = BigInt::from(a).pow(3) * BigInt::from(CHUDNOVSKY_C3_OVER_24);
            (p, q)
        };

        let mut t = &p * BigInt::from(13_591_409_u64 + 545_140_134_u64 * a as u64);
        if a & 1 != 0 {
            t = -t;
        }
        return ChudnovskySplit { p, q, t };
    }

    let mid = (a + b) / 2;
    let left = chudnovsky_split(a, mid);
    let right = chudnovsky_split(mid, b);
    let p = &left.p * &right.p;
    let q = &left.q * &right.q;
    let t = left.t * right.q + left.p * right.t;
    ChudnovskySplit { p, q, t }
}

fn scaled_two_over_pi(bits: usize, guard_bits: usize) -> BigUint {
    let work_bits = bits + guard_bits;
    let terms = work_bits.div_ceil(CHUDNOVSKY_BITS_PER_TERM) + 2;
    let split = chudnovsky_split(0, terms);
    let q = split.q.to_biguint().expect("Chudnovsky Q must be positive");
    let t = split.t.to_biguint().expect("Chudnovsky T must be positive");

    // pi = 426880 * sqrt(10005) * Q / T. Compute sqrt(10005)
    // as a fixed-point integer, then form floor((2/pi) * 2^bits).
    let sqrt_radicand = BigUint::from(10_005_u32) << (2 * work_bits);
    let sqrt_scaled = sqrt_radicand.sqrt();
    let numerator = t << (bits + 1 + work_bits);
    let denominator = q * BigUint::from(426_880_u32) * sqrt_scaled;
    numerator / denominator
}

fn two_over_pi_limbs() -> Vec<u64> {
    let low = scaled_two_over_pi(TWO_OVER_PI_BITS, LOW_GUARD_BITS);
    let high = scaled_two_over_pi(TWO_OVER_PI_BITS, HIGH_GUARD_BITS);
    assert_eq!(
        low, high,
        "2/pi table did not stabilize across independent guard precisions"
    );

    let mut limbs = high.to_u64_digits();
    assert!(
        limbs.len() <= TWO_OVER_PI_LIMBS,
        "2/pi table exceeded its requested bit width"
    );
    limbs.resize(TWO_OVER_PI_LIMBS, 0);
    limbs
}

fn render_constants() -> String {
    let limbs = two_over_pi_limbs();
    let mut output = String::new();
    writeln!(output, "const TWO_OVER_PI_BITS: usize = 16_640;").unwrap();
    writeln!(output, "const TWO_OVER_PI: [u64; {}] = [", limbs.len()).unwrap();
    for limb in limbs {
        writeln!(output, "    0x{limb:016x},").unwrap();
    }
    output.push_str("];");
    output
}

pub(super) fn render_source(template: &str) -> String {
    let start = template
        .find(CONSTANTS_START)
        .expect("Payne-Hanek template is missing the 2/pi constants");
    let end = template[start..]
        .find(CONSTANTS_END)
        .map(|offset| start + offset)
        .expect("Payne-Hanek template is missing PRODUCT_LIMBS");

    let constants = render_constants();
    let mut output = String::with_capacity(template.len());
    output.push_str(&template[..start]);
    output.push_str(&constants);
    output.push_str(&template[end..]);
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generated_table_has_the_known_leading_hex_digits() {
        let limbs = two_over_pi_limbs();
        assert_eq!(
            &limbs[TWO_OVER_PI_LIMBS - 4..],
            &[
                0xfe51_63ab_debb_c561,
                0xdb62_9599_3c43_9041,
                0xfc27_57d1_f534_ddc0,
                0xa2f9_836e_4e44_1529,
            ]
        );
    }

    #[test]
    fn generated_source_matches_aether_fpx() {
        let manifest = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let template = std::fs::read_to_string(manifest.join("templates/payne_hanek.rs"))
            .expect("failed to read Payne-Hanek template");
        let checked_in =
            std::fs::read_to_string(manifest.join("../../src/trig/payne_hanek.rs"))
                .expect("failed to read checked-in Payne-Hanek source");
        assert_eq!(render_source(&template), checked_in);
    }
}
