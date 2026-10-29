#[inline]
fn mirror_f64(x: f64) -> f64 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}
#[inline]
fn mirror_f32(x: f32) -> f32 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}
