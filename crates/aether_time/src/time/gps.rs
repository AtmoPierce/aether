use crate::time::Time;
use crate::time::GPS;

// GPS
pub const GPS_WEEK_SECONDS: i64 = 604_800;
pub const GPS_WEEK_MODULUS: u32 = 1024;

#[derive(Copy, Clone, Debug)]
pub struct GpsWeekTime {
    pub raw_week: u16,   // 0..1023
    pub sow: i64,        // seconds of week
    pub full_week: u32,  // rollover-corrected
}

impl Time<GPS> {
    pub fn to_gps_week(self, current_full_week_hint: u32) -> GpsWeekTime {
        let total = self.sec;
        let raw_week = ((total / GPS_WEEK_SECONDS) % GPS_WEEK_MODULUS as i64) as u16;
        let sow = total % GPS_WEEK_SECONDS;

        // Compute the nearest rollover-corrected week
        let raw = raw_week as i64;
        let hint = current_full_week_hint as i64;
        let base = (hint / GPS_WEEK_MODULUS as i64) * GPS_WEEK_MODULUS as i64;

        // choose nearest rollover window
        let mut full = base + raw;
        if (full - hint).abs() > (full + GPS_WEEK_MODULUS as i64 - hint).abs() {
            full += GPS_WEEK_MODULUS as i64;
        } else if (full - GPS_WEEK_MODULUS as i64 - hint).abs() < (full - hint).abs() {
            full -= GPS_WEEK_MODULUS as i64;
        }

        GpsWeekTime {
            raw_week,
            sow,
            full_week: full as u32,
        }
    }
}

#[cfg(feature = "defmt")]
impl defmt::Format for GpsWeekTime {
    fn format(&self, f: defmt::Formatter) {
        defmt::write!(
            f,
            "W{} (raw {}), +{}s",
            self.full_week, self.raw_week, self.sow
        );
    }
}