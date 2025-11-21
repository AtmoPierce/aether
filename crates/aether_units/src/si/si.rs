use aether_core::real::Real;
use crate::Quantity;
use crate::si::si_base::*;
use crate::si::si_composites::*;
use crate::si::si_derived::*;
use crate::si::si_prefix::*;

// Single trait for all SI extensions
pub trait SiExt: Real {
    // base
    fn m(self)   -> Quantity<Meter, Self>;
    fn kg(self)  -> Quantity<Kilogram, Self>;
    fn s(self)   -> Quantity<Second, Self>;
    fn A(self)   -> Quantity<Ampere, Self>;
    fn K(self)   -> Quantity<Kelvin, Self>;
    fn mol(self) -> Quantity<Mole, Self>;
    fn cd(self)  -> Quantity<Candela, Self>;

    // derived
    fn Hz(self)  -> Quantity<Hertz, Self>;
    fn N(self)   -> Quantity<Newton, Self>;
    fn Pa(self)  -> Quantity<Pascal, Self>;
    fn J(self)   -> Quantity<Joule, Self>;
    fn W(self)   -> Quantity<Watt, Self>;

    fn rad(self) -> Quantity<Radian, Self>;

    // composites
    fn m_per_s(self)    -> Quantity<MeterPerSecond, Self>;
    fn m_per_s2(self)   -> Quantity<MeterPerSecondSquared, Self>;
    fn rad_per_s(self)  -> Quantity<RadianPerSecond, Self>;
    fn rad_per_s2(self) -> Quantity<RadianPerSecondSquared, Self>;
}

// Below here builds all of the derived units (mill->mega)
macro_rules! si_ext_ctor {
    ($fn:ident, $Unit:ident) => {
        #[inline]
        fn $fn(self) -> Quantity<$Unit, Self> {
            Quantity::new(self)
        }
    };
}

macro_rules! si_ext_ctors {
    ( $( ($fn:ident, $Unit:ident) ),* $(,)? ) => {
        $(
            si_ext_ctor!($fn, $Unit);
        )*
    };
}

impl<T: Real> SiExt for T {
    // base units
    si_ext_ctors! {
        (m,   Meter),
        (kg,  Kilogram),
        (s,   Second),
        (A,   Ampere),
        (K,   Kelvin),
        (mol, Mole),
        (cd,  Candela),
    }

    // derived
    si_ext_ctors! {
        (Hz,  Hertz),
        (N,   Newton),
        (Pa,  Pascal),
        (J,   Joule),
        (W,   Watt),
        (rad, Radian),
    }

    // composites
    si_ext_ctors! {
        (m_per_s,    MeterPerSecond),
        (m_per_s2,   MeterPerSecondSquared),
        (rad_per_s,  RadianPerSecond),
        (rad_per_s2, RadianPerSecondSquared),
    }
}


pub trait SiPrefixExt: Real {
    // Length prefixes
    fn yoctometer(self) -> Quantity<Yoctometer, Self>;
    fn zeptometer(self) -> Quantity<Zeptometer, Self>;
    fn attometer(self)  -> Quantity<Attometer,  Self>;
    fn femtometer(self) -> Quantity<Femtometer, Self>;
    fn picometer(self)  -> Quantity<Picometer,  Self>;
    fn nanometer(self)  -> Quantity<Nanometer,  Self>;
    fn micrometer(self) -> Quantity<Micrometer, Self>;
    fn millimeter(self) -> Quantity<Millimeter, Self>;
    fn centimeter(self) -> Quantity<Centimeter, Self>;
    fn decimeter(self)  -> Quantity<Decimeter,  Self>;

    fn decameter(self)  -> Quantity<Decameter,  Self>;
    fn hectometer(self) -> Quantity<Hectometer, Self>;
    fn kilometer(self)  -> Quantity<Kilometer,  Self>;
    fn megameter(self)  -> Quantity<Megameter,  Self>;
    fn gigameter(self)  -> Quantity<Gigameter,  Self>;
    fn terameter(self)  -> Quantity<Terameter,  Self>;
    fn petameter(self)  -> Quantity<Petameter,  Self>;
    fn exameter(self)   -> Quantity<Exameter,   Self>;
    fn zettameter(self) -> Quantity<Zettameter, Self>;
    fn yottameter(self) -> Quantity<Yottameter, Self>;

    // Time prefixes
    fn yoctosecond(self) -> Quantity<Yoctosecond, Self>;
    fn zeptosecond(self) -> Quantity<Zeptosecond, Self>;
    fn attosecond(self)  -> Quantity<Attosecond,  Self>;
    fn femtosecond(self) -> Quantity<Femtosecond, Self>;
    fn picosecond(self)  -> Quantity<Picosecond,  Self>;
    fn nanosecond(self)  -> Quantity<Nanosecond,  Self>;
    fn microsecond(self) -> Quantity<Microsecond, Self>;
    fn millisecond(self) -> Quantity<Millisecond, Self>;
    fn centisecond(self) -> Quantity<Centisecond, Self>;
    fn decisecond(self)  -> Quantity<Decisecond,  Self>;

    fn decasecond(self)  -> Quantity<Decasecond,  Self>;
    fn hectosecond(self) -> Quantity<Hectosecond, Self>;
    fn kilosecond(self)  -> Quantity<Kilosecond,  Self>;
    fn megasecond(self)  -> Quantity<Megasecond,  Self>;
    fn gigasecond(self)  -> Quantity<Gigasecond,  Self>;
    fn terasecond(self)  -> Quantity<Terasecond,  Self>;
    fn petasecond(self)  -> Quantity<Petasecond,  Self>;
    fn exasecond(self)   -> Quantity<Exasecond,   Self>;
    fn zettasecond(self) -> Quantity<Zettasecond, Self>;
    fn yottasecond(self) -> Quantity<Yottasecond, Self>;

    // Mass prefixes (gram-based, mapped onto Kilogram base)
    fn yoctogram(self)  -> Quantity<Yoctogram,  Self>;
    fn zeptogram(self)  -> Quantity<Zeptogram,  Self>;
    fn attogram(self)   -> Quantity<Attogram,   Self>;
    fn femtogram(self)  -> Quantity<Femtogram,  Self>;
    fn picogram(self)   -> Quantity<Picogram,   Self>;
    fn nanogram(self)   -> Quantity<Nanogram,   Self>;
    fn microgram(self)  -> Quantity<Microgram,  Self>;
    fn milligram(self)  -> Quantity<Milligram,  Self>;
    fn centigram(self)  -> Quantity<Centigram,  Self>;
    fn decigram(self)   -> Quantity<Decigram,   Self>;

    fn decagram(self)   -> Quantity<Decagram,   Self>;
    fn hectogram(self)  -> Quantity<Hectogram,  Self>;
    fn megagram(self)   -> Quantity<Megagram,   Self>;
    fn gigagram(self)   -> Quantity<Gigagram,   Self>;
    fn teragram(self)   -> Quantity<Teragram,   Self>;
    fn petagram(self)   -> Quantity<Petagram,   Self>;
    fn exagram(self)    -> Quantity<Exagram,    Self>;
    fn zettagram(self)  -> Quantity<Zettagram,  Self>;
    fn yottagram(self)  -> Quantity<Yottagram,  Self>;
}

macro_rules! si_prefix_ctor {
    ($fn:ident, $Unit:ident) => {
        #[inline]
        fn $fn(self) -> Quantity<$Unit, Self> {
            Quantity::new(self)
        }
    };
}

macro_rules! si_prefix_ctors {
    ( $( ($fn:ident, $Unit:ident) ),* $(,)? ) => {
        $(
            si_prefix_ctor!($fn, $Unit);
        )*
    };
}

impl<T: Real> SiPrefixExt for T {
    si_prefix_ctors! {
        // Small length prefixes
        (yoctometer, Yoctometer),
        (zeptometer, Zeptometer),
        (attometer,  Attometer),
        (femtometer, Femtometer),
        (picometer,  Picometer),
        (nanometer,  Nanometer),
        (micrometer, Micrometer),
        (millimeter, Millimeter),
        (centimeter, Centimeter),
        (decimeter,  Decimeter),

        // Large length prefixes
        (decameter,   Decameter),
        (hectometer,  Hectometer),
        (kilometer,   Kilometer),
        (megameter,   Megameter),
        (gigameter,   Gigameter),
        (terameter,   Terameter),
        (petameter,   Petameter),
        (exameter,    Exameter),
        (zettameter,  Zettameter),
        (yottameter,  Yottameter),
    }

    si_prefix_ctors! {
        // Small time prefixes
        (yoctosecond, Yoctosecond),
        (zeptosecond, Zeptosecond),
        (attosecond,  Attosecond),
        (femtosecond, Femtosecond),
        (picosecond,  Picosecond),
        (nanosecond,  Nanosecond),
        (microsecond, Microsecond),
        (millisecond, Millisecond),
        (centisecond, Centisecond),
        (decisecond,  Decisecond),

        // Large time prefixes
        (decasecond,   Decasecond),
        (hectosecond,  Hectosecond),
        (kilosecond,   Kilosecond),
        (megasecond,   Megasecond),
        (gigasecond,   Gigasecond),
        (terasecond,   Terasecond),
        (petasecond,   Petasecond),
        (exasecond,    Exasecond),
        (zettasecond,  Zettasecond),
        (yottasecond,  Yottasecond),
    }

    si_prefix_ctors! {
        // Mass prefixes (gram-based prefixes on Kilogram base)
        (yoctogram,  Yoctogram),
        (zeptogram,  Zeptogram),
        (attogram,   Attogram),
        (femtogram,  Femtogram),
        (picogram,   Picogram),
        (nanogram,   Nanogram),
        (microgram,  Microgram),
        (milligram,  Milligram),
        (centigram,  Centigram),
        (decigram,   Decigram),

        (decagram,   Decagram),
        (hectogram,  Hectogram),
        (megagram,   Megagram),
        (gigagram,   Gigagram),
        (teragram,   Teragram),
        (petagram,   Petagram),
        (exagram,    Exagram),
        (zettagram,  Zettagram),
        (yottagram,  Yottagram),
    }
}
