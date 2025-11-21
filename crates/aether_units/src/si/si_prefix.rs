// ----------------- SI prefixes & prefixed units -----------------
use aether_core::real::Real;
use crate::si::si_base::*;

// SI prefix factors as f32
pub const YOCTO: f32 = 1.0e-24_f32;
pub const ZEPTO: f32 = 1.0e-21_f32;
pub const ATTO:  f32 = 1.0e-18_f32;
pub const FEMTO: f32 = 1.0e-15_f32;
pub const PICO:  f32 = 1.0e-12_f32;
pub const NANO:  f32 = 1.0e-9_f32;
pub const MICRO: f32 = 1.0e-6_f32;
pub const MILLI: f32 = 1.0e-3_f32;
pub const CENTI: f32 = 1.0e-2_f32;
pub const DECI:  f32 = 1.0e-1_f32;
pub const DECA:  f32 = 1.0e1_f32;
pub const HECTO: f32 = 1.0e2_f32;
pub const KILO:  f32 = 1.0e3_f32;
pub const MEGA:  f32 = 1.0e6_f32;
pub const GIGA:  f32 = 1.0e9_f32;
pub const TERA:  f32 = 1.0e12_f32;
pub const PETA:  f32 = 1.0e15_f32;
pub const EXA:   f32 = 1.0e18_f32;
pub const ZETTA: f32 = 1.0e21_f32;
pub const YOTTA: f32 = 1.0e24_f32;

macro_rules! si_prefixed_unit {
    ($PrefixUnit:ident, $BaseUnit:ident, $factor:expr,
     $to_prefixed:ident, $to_base:ident) =>
    {
        #[derive(Copy, Clone, Debug, PartialEq, Eq)]
        pub enum $PrefixUnit {}

        // Base -> Prefixed
        impl<T: Real> $crate::Quantity<$BaseUnit, T> {
            #[inline]
            pub fn $to_prefixed(self) -> $crate::Quantity<$PrefixUnit, T> {
                let k = T::from_f32($factor);
                $crate::Quantity::new(self.raw() / k)
            }
        }

        // Prefixed -> Base
        impl<T: Real> $crate::Quantity<$PrefixUnit, T> {
            #[inline]
            pub fn $to_base(self) -> $crate::Quantity<$BaseUnit, T> {
                let k = T::from_f32($factor);
                $crate::Quantity::new(self.raw() * k)
            }
        }
    };
}

/// Define a family of prefixed units for a given base.
macro_rules! si_prefixed_family {
    (
        base = $BaseUnit:ident;
        to_base = $to_base:ident;
        $(
            ($PrefixUnit:ident, $factor:expr, $to_prefixed:ident)
        ),* $(,)?
    ) => {
        $(
            si_prefixed_unit!($PrefixUnit, $BaseUnit, $factor, $to_prefixed, $to_base);
        )*
    };
}

// Length (meter) – yocto to yotta
si_prefixed_family! {
    base = Meter;
    to_base = to_meter;

    // tiny
    (Yoctometer,  YOCTO, to_yoctometer),
    (Zeptometer,  ZEPTO, to_zeptometer),
    (Attometer,   ATTO,  to_attometer),
    (Femtometer,  FEMTO, to_femtometer),
    (Picometer,   PICO,  to_picometer),
    (Nanometer,   NANO,  to_nanometer),
    (Micrometer,  MICRO, to_micrometer),
    (Millimeter,  MILLI, to_millimeter),
    (Centimeter,  CENTI, to_centimeter),
    (Decimeter,   DECI,  to_decimeter),

    // > 1 m
    (Decameter,   DECA,  to_decameter),
    (Hectometer,  HECTO, to_hectometer),
    (Kilometer,   KILO,  to_kilometer),
    (Megameter,   MEGA,  to_megameter),
    (Gigameter,   GIGA,  to_gigameter),
    (Terameter,   TERA,  to_terameter),
    (Petameter,   PETA,  to_petameter),
    (Exameter,    EXA,   to_exameter),
    (Zettameter,  ZETTA, to_zettameter),
    (Yottameter,  YOTTA, to_yottameter),
}

// Time (second) – yocto to yotta
si_prefixed_family! {
    base = Second;
    to_base = to_second;

    (Yoctosecond,  YOCTO, to_yoctosecond),
    (Zeptosecond,  ZEPTO, to_zeptosecond),
    (Attosecond,   ATTO,  to_attosecond),
    (Femtosecond,  FEMTO, to_femtosecond),
    (Picosecond,   PICO,  to_picosecond),
    (Nanosecond,   NANO,  to_nanosecond),
    (Microsecond,  MICRO, to_microsecond),
    (Millisecond,  MILLI, to_millisecond),
    (Centisecond,  CENTI, to_centisecond),
    (Decisecond,   DECI,  to_decisecond),

    (Decasecond,   DECA,  to_decasecond),
    (Hectosecond,  HECTO, to_hectosecond),
    (Kilosecond,   KILO,  to_kilosecond),
    (Megasecond,   MEGA,  to_megasecond),
    (Gigasecond,   GIGA,  to_gigasecond),
    (Terasecond,   TERA,  to_terasecond),
    (Petasecond,   PETA,  to_petasecond),
    (Exasecond,    EXA,   to_exasecond),
    (Zettasecond,  ZETTA, to_zettasecond),
    (Yottasecond,  YOTTA, to_yottasecond),
}

// Mass (kilogram) – prefixes conceptually apply to gram, so factor * 1e-3
si_prefixed_family! {
    base = Kilogram;
    to_base = to_kilogram;

    (Yoctogram,  YOCTO * 1.0e-3_f32, to_yoctogram),
    (Zeptogram,  ZEPTO * 1.0e-3_f32, to_zeptogram),
    (Attogram,   ATTO  * 1.0e-3_f32, to_attogram),
    (Femtogram,  FEMTO * 1.0e-3_f32, to_femtogram),
    (Picogram,   PICO  * 1.0e-3_f32, to_picogram),
    (Nanogram,   NANO  * 1.0e-3_f32, to_nanogram),
    (Microgram,  MICRO * 1.0e-3_f32, to_microgram),
    (Milligram,  MILLI * 1.0e-3_f32, to_milligram),
    (Centigram,  CENTI * 1.0e-3_f32, to_centigram),
    (Decigram,   DECI  * 1.0e-3_f32, to_decigram),

    (Decagram,   DECA  * 1.0e-3_f32, to_decagram),
    (Hectogram,  HECTO * 1.0e-3_f32, to_hectogram),
    (Megagram,   MEGA  * 1.0e-3_f32, to_megagram),
    (Gigagram,   GIGA  * 1.0e-3_f32, to_gigagram),
    (Teragram,   TERA  * 1.0e-3_f32, to_teragram),
    (Petagram,   PETA  * 1.0e-3_f32, to_petagram),
    (Exagram,    EXA   * 1.0e-3_f32, to_exagram),
    (Zettagram,  ZETTA * 1.0e-3_f32, to_zettagram),
    (Yottagram,  YOTTA * 1.0e-3_f32, to_yottagram),
}
