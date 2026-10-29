#[inline]
pub fn exp2i(k: i32) -> f32 {
    let e = (k + 127).clamp(1, 254);
    f32::from_bits((e as u32) << 23)
}
