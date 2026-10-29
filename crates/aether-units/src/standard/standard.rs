// src/units/standard.rs
#![cfg(feature = "units")]

use super::{Unit};
use super::metric::{Meter, KgPerCubicMeter, Kelvin, Celsius};

/* ---------- Standard / US customary unit markers ---------- */

// Length
pub enum Inch {}
pub enum Foot {}
pub enum Yard {}
pub enum Mile {}

// Pressure
pub enum Psi {}        // pounds per square inch

// Temperature
pub enum Fahrenheit {}

// Density
pub enum SlugPerCubicFoot {}

/* ---------- Construction on f64 ---------- */

pub trait StandardExt {
    fn inch(self) -> Unit<Inch>;
    fn ft(self)   -> Unit<Foot>;
    fn yd(self)   -> Unit<Yard>;
    fn mile(self) -> Unit<Mile>;
}

impl StandardExt for f64 {
    #[inline] fn inch(self) -> Unit<Inch> { Unit::new(self) }
    #[inline] fn ft(self)   -> Unit<Foot> { Unit::new(self) }
    #[inline] fn yd(self)   -> Unit<Yard> { Unit::new(self) }
    #[inline] fn mile(self) -> Unit<Mile> { Unit::new(self) }
}

/* ---------- Standard <-> metric conversions ---------- */

impl Unit<Inch> {
    #[inline]
    pub fn to_feet(self) -> Unit<Foot> {
        Unit::new(self.raw() / 12.0)
    }

    #[inline]
    pub fn to_meters(self) -> Unit<Meter> {
        // 1 inch = 0.0254 m
        Unit::new(self.raw() * 0.0254)
    }
    #[inline] fn psi(self)  -> Unit<Psi>{ 
        Unit::new(self) 
    }
    #[inline] fn F(self)    -> Unit<Fahrenheit> { 
        Unit::new(self) 
    }
    #[inline] fn slug_per_ft3(self) -> Unit<SlugPerCubicFoot> {
        Unit::new(self)
    }
}

impl Unit<Foot> {
    #[inline]
    pub fn to_inches(self) -> Unit<Inch> {
        Unit::new(self.raw() * 12.0)
    }

    #[inline]
    pub fn to_yards(self) -> Unit<Yard> {
        Unit::new(self.raw() / 3.0)
    }

    #[inline]
    pub fn to_meters(self) -> Unit<Meter> {
        // 1 ft = 0.3048 m
        Unit::new(self.raw() * 0.3048)
    }
}

impl Unit<Yard> {
    #[inline]
    pub fn to_feet(self) -> Unit<Foot> {
        Unit::new(self.raw() * 3.0)
    }

    #[inline]
    pub fn to_meters(self) -> Unit<Meter> {
        Unit::new(self.raw() * 0.3048*3)
    }
}

impl Unit<Mile> {
    #[inline]
    pub fn to_feet(self) -> Unit<Foot> {
        // 1 mile = 5280 ft
        Unit::new(self.raw() * 5280.0)
    }

    #[inline]
    pub fn to_meters(self) -> Unit<Meter> {
        // 1 mile = 1609.344 m
        Unit::new(self.raw() * 1609.344)
    }
}

impl Unit<Psi> {
    #[inline]
    pub fn to_pascal(self) -> Unit<super::metric::Pascal> {
        // 1 psi = 6894.757293168 Pa (approx)
        Unit::new(self.raw() * 6894.757_293_168)
    }

    #[inline]
    pub fn to_kilopascal(self) -> Unit<super::metric::KiloPascal> {
        self.to_pascal().to_kilopascal()
    }

    #[inline]
    pub fn to_bar(self) -> Unit<super::metric::Bar> {
        self.to_pascal().to_bar()
    }

    #[inline]
    pub fn to_atm(self) -> Unit<super::metric::Atmosphere> {
        self.to_pascal().to_atm()
    }
}

impl Unit<Fahrenheit> {
    #[inline]
    pub fn to_celsius(self) -> Unit<Celsius> {
        // C = (F - 32) * 5/9
        Unit::new((self.raw() - 32.0) * (5.0 / 9.0))
    }

    #[inline]
    pub fn to_kelvin(self) -> Unit<Kelvin> {
        self.to_celsius().to_kelvin()
    }
}

impl Unit<Celsius> {
    #[inline]
    pub fn to_fahrenheit(self) -> Unit<Fahrenheit> {
        // F = C * 9/5 + 32
        Unit::new(self.raw() * (9.0 / 5.0) + 32.0)
    }
}

impl Unit<Kelvin> {
    #[inline]
    pub fn to_fahrenheit(self) -> Unit<Fahrenheit> {
        self.to_celsius().to_fahrenheit()
    }
}
