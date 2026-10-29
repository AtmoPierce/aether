pub(crate) fn ordered_u32(x: f32) -> u32 {
    let bits = x.to_bits();
    if bits & 0x8000_0000 != 0 {
        !bits
    } else {
        bits | 0x8000_0000
    }
}

pub(crate) fn ordered_u64(x: f64) -> u64 {
    let bits = x.to_bits();
    if bits & 0x8000_0000_0000_0000 != 0 {
        !bits
    } else {
        bits | 0x8000_0000_0000_0000
    }
}

pub(crate) fn ordered_u128(x: f128) -> u128 {
    let bits = x.to_bits();
    if bits & (1_u128 << 127) != 0 {
        !bits
    } else {
        bits | (1_u128 << 127)
    }
}

pub(crate) fn ulp_diff_f32(a: f32, b: f32) -> u32 {
    ordered_u32(a).abs_diff(ordered_u32(b))
}

pub(crate) fn ulp_diff_f64(a: f64, b: f64) -> u64 {
    ordered_u64(a).abs_diff(ordered_u64(b))
}

pub(crate) fn ulp_diff_f128(a: f128, b: f128) -> u128 {
    ordered_u128(a).abs_diff(ordered_u128(b))
}
