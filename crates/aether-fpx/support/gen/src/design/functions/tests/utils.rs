#[inline]
pub fn ulp_diff_f64(a: f64, b: f64) -> u64 {
    let ia = a.to_bits();
    let ib = b.to_bits();

    fn ordered(u: u64) -> i64 {
        if u & 0x8000_0000_0000_0000 != 0 {
            (!u) as i64
        } else {
            u as i64
        }
    }

    let da = ordered(ia);
    let db = ordered(ib);
    da.wrapping_sub(db).abs() as u64
}

#[cfg(feature = "f128")]
#[inline]
pub fn ulp_diff_f128(a: f128, b: f128) -> u128 {
    let ia = a.to_bits();
    let ib = b.to_bits();

    fn ordered(u: u128) -> i128 {
        if u & 0x8000_0000_0000_0000_0000_0000_0000_0000 != 0 {
            (!u) as i128
        } else {
            u as i128
        }
    }

    let da = ordered(ia);
    let db = ordered(ib);
    let d = da.wrapping_sub(db);
    if d < 0 {
        (-d) as u128
    } else {
        d as u128
    }
}
