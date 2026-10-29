/// Trait to convert an angle-like value to radians.
pub trait ToRadians {
    type Output;
    fn to_radians(self) -> Self::Output;
}

/// Trait to convert an angle-like value to degrees.
pub trait ToDegrees {
    type Output;
    fn to_degrees(self) -> Self::Output;
}

// Implementations for primitive Reals

impl ToRadians for f32 {
    type Output = f32;
    fn to_radians(self) -> Self::Output {
        self * core::f32::consts::PI / 180.0
    }
}

impl ToDegrees for f32 {
    type Output = f32;
    fn to_degrees(self) -> Self::Output {
        self * 180.0 / core::f32::consts::PI
    }
}

impl ToRadians for f64 {
    type Output = f64;
    fn to_radians(self) -> Self::Output {
        self * core::f64::consts::PI / 180.0
    }
}

impl ToDegrees for f64 {
    type Output = f64;
    fn to_degrees(self) -> Self::Output {
        self * 180.0 / core::f64::consts::PI
    }
}

impl<T, const N: usize> ToRadians for [T; N]
where
    T: ToRadians + Copy,
{
    type Output = [T::Output; N];
    fn to_radians(self) -> Self::Output {
        self.map(|x| x.to_radians())
    }
}

impl<T, const N: usize> ToDegrees for [T; N]
where
    T: ToDegrees + Copy,
{
    type Output = [T::Output; N];
    fn to_degrees(self) -> Self::Output {
        self.map(|x| x.to_degrees())
    }
}
