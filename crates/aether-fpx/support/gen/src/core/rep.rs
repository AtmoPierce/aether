// I wanted to make these common but makes a little more sense this way.
#[derive(Copy, Clone, Debug)]
pub struct UnpackedE2 {
    pub sign: i8, // +1 or -1
    pub exp: i8,  // unbiased exponent
    pub mant: u8, // Representation is dependent on fp type
    pub is_zero: bool,
}

#[derive(Copy, Clone, Debug)]
pub struct UnpackedE3 {
    pub sign: i8,  // +1 or -1
    pub exp: i16,  // unbiased exponent
    pub mant: u16, // Representation is dependent on fp type
    pub is_zero: bool,
    pub is_subnormal: bool,
}
#[derive(Copy, Clone, Debug)]
pub struct UnpackedE4 {
    pub sign: i8,           // +1 or -1
    pub exp: i16,           // unbiased exponent
    pub mant: u16,          // Representation is dependent on fp type
    pub is_zero: bool,      // Convenience
    pub is_subnormal: bool, // Convenience
    pub is_nan: bool,       // Convenience
}
#[derive(Copy, Clone, Debug)]
pub struct UnpackedE5 {
    pub sign: i8,           // +1 or -1
    pub exp: i16,           // unbiased exponent
    pub mant: u16,          // Representation is dependent on fp type
    pub is_zero: bool,      // Convenience
    pub is_subnormal: bool, // Convenience
    pub is_inf: bool,       // Convenience
    pub is_nan: bool,       // Convenience
}
#[derive(Copy, Clone, Debug)]
pub struct Unpacked6 {
    pub sign: i8,           // +1 or -1
    pub exp: i8,            // unbiased exponent
    pub mant: u16,          // Representation is dependent on fp type
    pub is_zero: bool,      // Convenience
    pub is_subnormal: bool, // Convenience
}
