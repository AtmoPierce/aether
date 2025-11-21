#[cfg(test)]
mod tests {
    use aether_core::real::Real;

    use crate::Quantity;
    use crate::si::si_base::*;
    use crate::si::si_composites::*;
    use crate::si::si_derived::*;
    use crate::si::si_prefix::*;
    use crate::si::{SiExt, SiPrefixExt}; // assuming you export both traits

    // Tiny helper to compare f32 with a loose epsilon
    fn feq(a: f32, b: f32, eps: f32) -> bool {
        (a - b).abs() <= eps * b.abs().max(1.0)
    }

    macro_rules! check_ctor {
        ($val:expr, $method:ident, $Unit:ty) => {{
            let x: f32 = $val;
            let q: Quantity<$Unit, f32> = x.$method();
            assert!(
                feq(q.raw(), x, 1e-6),
                "ctor {}(): raw={} expected={}",
                stringify!($method),
                q.raw(),
                x
            );
        }};
    }

    macro_rules! check_prefix_ctor {
        ($val:expr,
         $method:ident,
         $Unit:ty,
         $to_base:ident,
         $BaseUnit:ty,
         $scale:expr) => {{
            let x: f32 = $val;
            let q: Quantity<$Unit, f32> = x.$method();
            let base: Quantity<$BaseUnit, f32> = q.$to_base();
            let expected = x * $scale;
            assert!(
                feq(base.raw(), expected, 1e-5),
                "prefix {}(): {} -> base raw={} expected={}",
                stringify!($method),
                stringify!($Unit),
                base.raw(),
                expected
            );
        }};
    }

    #[test]
    fn si_ext_and_prefix_ctors_work() {
        let x = 1.0_f32;

        // -------- base units --------
        check_ctor!(x, m,   Meter);
        check_ctor!(x, kg,  Kilogram);
        check_ctor!(x, s,   Second);
        check_ctor!(x, A,   Ampere);
        check_ctor!(x, K,   Kelvin);
        check_ctor!(x, mol, Mole);
        check_ctor!(x, cd,  Candela);

        // -------- derived units --------
        check_ctor!(x, Hz,  Hertz);
        check_ctor!(x, N,   Newton);
        check_ctor!(x, Pa,  Pascal);
        check_ctor!(x, J,   Joule);
        check_ctor!(x, W,   Watt);
        check_ctor!(x, rad, Radian);

        // -------- composite units --------
        check_ctor!(x, m_per_s,    MeterPerSecond);
        check_ctor!(x, m_per_s2,   MeterPerSecondSquared);
        check_ctor!(x, rad_per_s,  RadianPerSecond);
        check_ctor!(x, rad_per_s2, RadianPerSecondSquared);

        // =====================================================================
        // Prefix units: meter family (-> Meter, factor = prefix)
        // =====================================================================

        // tiny
        check_prefix_ctor!(x, yoctometer, Yoctometer, to_meter, Meter, YOCTO);
        check_prefix_ctor!(x, zeptometer, Zeptometer, to_meter, Meter, ZEPTO);
        check_prefix_ctor!(x, attometer,  Attometer,  to_meter, Meter, ATTO);
        check_prefix_ctor!(x, femtometer, Femtometer, to_meter, Meter, FEMTO);
        check_prefix_ctor!(x, picometer,  Picometer,  to_meter, Meter, PICO);
        check_prefix_ctor!(x, nanometer,  Nanometer,  to_meter, Meter, NANO);
        check_prefix_ctor!(x, micrometer, Micrometer, to_meter, Meter, MICRO);
        check_prefix_ctor!(x, millimeter, Millimeter, to_meter, Meter, MILLI);
        check_prefix_ctor!(x, centimeter, Centimeter, to_meter, Meter, CENTI);
        check_prefix_ctor!(x, decimeter,  Decimeter,  to_meter, Meter, DECI);

        // > 1 m
        check_prefix_ctor!(x, decameter,   Decameter,   to_meter, Meter, DECA);
        check_prefix_ctor!(x, hectometer,  Hectometer,  to_meter, Meter, HECTO);
        check_prefix_ctor!(x, kilometer,   Kilometer,   to_meter, Meter, KILO);
        check_prefix_ctor!(x, megameter,   Megameter,   to_meter, Meter, MEGA);
        check_prefix_ctor!(x, gigameter,   Gigameter,   to_meter, Meter, GIGA);
        check_prefix_ctor!(x, terameter,   Terameter,   to_meter, Meter, TERA);
        check_prefix_ctor!(x, petameter,   Petameter,   to_meter, Meter, PETA);
        check_prefix_ctor!(x, exameter,    Exameter,    to_meter, Meter, EXA);
        check_prefix_ctor!(x, zettameter,  Zettameter,  to_meter, Meter, ZETTA);
        check_prefix_ctor!(x, yottameter,  Yottameter,  to_meter, Meter, YOTTA);

        // =====================================================================
        // Prefix units: second family (-> Second, factor = prefix)
        // =====================================================================

        // tiny
        check_prefix_ctor!(x, yoctosecond, Yoctosecond, to_second, Second, YOCTO);
        check_prefix_ctor!(x, zeptosecond, Zeptosecond, to_second, Second, ZEPTO);
        check_prefix_ctor!(x, attosecond,  Attosecond,  to_second, Second, ATTO);
        check_prefix_ctor!(x, femtosecond, Femtosecond, to_second, Second, FEMTO);
        check_prefix_ctor!(x, picosecond,  Picosecond,  to_second, Second, PICO);
        check_prefix_ctor!(x, nanosecond,  Nanosecond,  to_second, Second, NANO);
        check_prefix_ctor!(x, microsecond, Microsecond, to_second, Second, MICRO);
        check_prefix_ctor!(x, millisecond, Millisecond, to_second, Second, MILLI);
        check_prefix_ctor!(x, centisecond, Centisecond, to_second, Second, CENTI);
        check_prefix_ctor!(x, decisecond,  Decisecond,  to_second, Second, DECI);

        // > 1 s
        check_prefix_ctor!(x, decasecond,   Decasecond,   to_second, Second, DECA);
        check_prefix_ctor!(x, hectosecond,  Hectosecond,  to_second, Second, HECTO);
        check_prefix_ctor!(x, kilosecond,   Kilosecond,   to_second, Second, KILO);
        check_prefix_ctor!(x, megasecond,   Megasecond,   to_second, Second, MEGA);
        check_prefix_ctor!(x, gigasecond,   Gigasecond,   to_second, Second, GIGA);
        check_prefix_ctor!(x, terasecond,   Terasecond,   to_second, Second, TERA);
        check_prefix_ctor!(x, petasecond,   Petasecond,   to_second, Second, PETA);
        check_prefix_ctor!(x, exasecond,    Exasecond,    to_second, Second, EXA);
        check_prefix_ctor!(x, zettasecond,  Zettasecond,  to_second, Second, ZETTA);
        check_prefix_ctor!(x, yottasecond,  Yottasecond,  to_second, Second, YOTTA);

        // =====================================================================
        // Prefix units: mass family (Kilogram base, gram prefixes => factor * 1e-3)
        // =====================================================================

        check_prefix_ctor!(x, yoctogram,  Yoctogram,  to_kilogram, Kilogram, YOCTO * 1.0e-3_f32);
        check_prefix_ctor!(x, zeptogram,  Zeptogram,  to_kilogram, Kilogram, ZEPTO * 1.0e-3_f32);
        check_prefix_ctor!(x, attogram,   Attogram,   to_kilogram, Kilogram, ATTO  * 1.0e-3_f32);
        check_prefix_ctor!(x, femtogram,  Femtogram,  to_kilogram, Kilogram, FEMTO * 1.0e-3_f32);
        check_prefix_ctor!(x, picogram,   Picogram,   to_kilogram, Kilogram, PICO  * 1.0e-3_f32);
        check_prefix_ctor!(x, nanogram,   Nanogram,   to_kilogram, Kilogram, NANO  * 1.0e-3_f32);
        check_prefix_ctor!(x, microgram,  Microgram,  to_kilogram, Kilogram, MICRO * 1.0e-3_f32);
        check_prefix_ctor!(x, milligram,  Milligram,  to_kilogram, Kilogram, MILLI * 1.0e-3_f32);
        check_prefix_ctor!(x, centigram,  Centigram,  to_kilogram, Kilogram, CENTI * 1.0e-3_f32);
        check_prefix_ctor!(x, decigram,   Decigram,   to_kilogram, Kilogram, DECI  * 1.0e-3_f32);

        check_prefix_ctor!(x, decagram,   Decagram,   to_kilogram, Kilogram, DECA  * 1.0e-3_f32);
        check_prefix_ctor!(x, hectogram,  Hectogram,  to_kilogram, Kilogram, HECTO * 1.0e-3_f32);
        check_prefix_ctor!(x, megagram,   Megagram,   to_kilogram, Kilogram, MEGA  * 1.0e-3_f32);
        check_prefix_ctor!(x, gigagram,   Gigagram,   to_kilogram, Kilogram, GIGA  * 1.0e-3_f32);
        check_prefix_ctor!(x, teragram,   Teragram,   to_kilogram, Kilogram, TERA  * 1.0e-3_f32);
        check_prefix_ctor!(x, petagram,   Petagram,   to_kilogram, Kilogram, PETA  * 1.0e-3_f32);
        check_prefix_ctor!(x, exagram,    Exagram,    to_kilogram, Kilogram, EXA   * 1.0e-3_f32);
        check_prefix_ctor!(x, zettagram,  Zettagram,  to_kilogram, Kilogram, ZETTA * 1.0e-3_f32);
        check_prefix_ctor!(x, yottagram,  Yottagram,  to_kilogram, Kilogram, YOTTA * 1.0e-3_f32);
    }

    #[test]
    fn addition(){
        let meter = 1.0.m();
        let km = 1.0.kilometer();
        let one_one = km + meter.to_kilometer();
        assert_eq!(one_one.value, 1.001);
        assert_eq!(1.001.kilometer(), one_one);
    }
    #[test]
    fn subtraction(){
        let meter = 1.0.m();
        let km = 1.0.kilometer();
        let one_one = km - meter.to_kilometer();
        assert_eq!(one_one.value, 0.999);
        assert_eq!(0.999.kilometer(), one_one);
    }
}
