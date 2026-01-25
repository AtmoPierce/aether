use core::cmp::Ordering;

use super::time::Time;
use super::duration::Duration;
use super::UTC;


#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct UtcDateTime {
    pub year: i32,
    pub month: u8,   // 1-12
    pub day: u8,     // 1-31
    pub hour: u8,
    pub minute: u8,
    pub second: u8,
    pub nanos: i32,  // 0..1_000_000_000
}

impl UtcDateTime {
    #[inline]
    pub fn new(
        year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: u8,
        nanos: i32,
    ) -> Self {
        UtcDateTime { year, month, day, hour, minute, second, nanos }
    }
}

// Generic timing...
fn days_to_ymd(days: i64, hour: u8, minute: u8, second: u8) -> UtcDateTime {
    // Days since Unix epoch
    let mut z = days + 719_468; // shift to civil-based epoch
    let era = if z >= 0 {
        z / 146_097
    } else {
        (z - 146_096) / 146_097
    };
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let mut y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = mp + if mp < 10 { 3 } else { -9 };
    y += (m <= 2) as i64;

    UtcDateTime {
        year:  y as i32,
        month: m as u8,
        day:   d as u8,
        hour,
        minute,
        second,
        nanos: 0, // days_to_ymd works at 1-second granularity
    }
}

#[derive(Copy, Clone, Debug)]
pub struct UtcAlarm {
    target: Time<UTC>,
    fired: bool,
}

impl UtcAlarm {
    pub fn new(target: Time<UTC>) -> Self {
        Self { target, fired: false }
    }

    /// Returns `true` exactly once when `now >= target`.
    pub fn check(&mut self, now: Time<UTC>) -> bool {
        if !self.fired && now >= self.target {
            self.fired = true;
            true
        } else {
            false
        }
    }

    pub fn is_done(&self) -> bool {
        self.fired
    }
}


#[derive(Copy, Clone, Debug)]
pub enum Recurrence {
    Once,
    Interval(Duration), // e.g. every 5s, 100ms, etc.
    Hourly { minute: u8, second: u8 },
    Daily  { hour: u8, minute: u8, second: u8 },
    Yearly { month: u8, day: u8, hour: u8, minute: u8, second: u8 },
}

#[derive(Copy, Clone, Debug)]
pub struct UtcSchedule {
    next_fire: Time<UTC>,
    pattern: Recurrence,
}


const SECS_PER_MIN: i64 = 60;
const SECS_PER_HOUR: i64 = 60 * SECS_PER_MIN;
const SECS_PER_DAY:  i64 = 24 * SECS_PER_HOUR;

fn next_interval(after: Time<UTC>, period: Duration) -> Time<UTC> {
    // basic: step forward at least one period, catching up if we’re far behind
    let mut t = after + period;
    // optional: if you want to catch up multiple periods, loop:
    // while t <= after { t = t + period; }
    t
}

fn next_hourly(after: Time<UTC>, minute: u8, second: u8) -> Time<UTC> {
    let now_sec = after.sec;
    let hour_start = (now_sec / SECS_PER_HOUR) * SECS_PER_HOUR;
    let target_in_hour =
        hour_start + (minute as i64) * SECS_PER_MIN + (second as i64);

    let mut next_sec = target_in_hour;
    if next_sec <= now_sec {
        next_sec += SECS_PER_HOUR;
    }
    Time::<UTC>::new(next_sec, 0)
}

fn next_daily(after: Time<UTC>, hour: u8, minute: u8, second: u8) -> Time<UTC> {
    let now_sec = after.sec;
    let day_start = (now_sec / SECS_PER_DAY) * SECS_PER_DAY;
    let target_today =
        day_start + (hour as i64) * SECS_PER_HOUR
                  + (minute as i64) * SECS_PER_MIN
                  + (second as i64);

    let mut next_sec = target_today;
    if next_sec <= now_sec {
        next_sec += SECS_PER_DAY;
    }
    Time::<UTC>::new(next_sec, 0)
}

fn next_yearly(after: Time<UTC>, month: u8, day: u8,
               hour: u8, minute: u8, second: u8) -> Time<UTC> {
    // Use your existing YMD conversion:
    let dt = after.to_utc_datetime();

    // first try this year
    let mut year = dt.year;
    let mut candidate_dt = UtcDateTime::new(
        year, month, day, hour, minute, second, 0,
    );
    let mut candidate = Time::<UTC>::from_utc_datetime(candidate_dt);

    if candidate <= after {
        // move to next year
        year += 1;
        candidate_dt = UtcDateTime::new(
            year, month, day, hour, minute, second, 0,
        );
        candidate = Time::<UTC>::from_utc_datetime(candidate_dt);
    }

    candidate
}


#[cfg(feature = "defmt")]
impl defmt::Format for UtcDateTime {
    fn format(&self, f: defmt::Formatter) {
        defmt::write!(
            f,
            "{:04}-{:02}-{:02} {:02}:{:02}:{:02}.{:09} UTC",
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
            self.nanos
        );
    }
}

impl UtcSchedule {
    /// Initialize from a reference `start_from` (usually "now").
    pub fn new(start_from: Time<UTC>, pattern: Recurrence) -> Self {
        let next_fire = match pattern {
            Recurrence::Once => start_from,
            Recurrence::Interval(period) => next_interval(start_from, period),
            Recurrence::Hourly { minute, second } =>
                next_hourly(start_from, minute, second),
            Recurrence::Daily { hour, minute, second } =>
                next_daily(start_from, hour, minute, second),
            Recurrence::Yearly { month, day, hour, minute, second } =>
                next_yearly(start_from, month, day, hour, minute, second),
        };

        Self { next_fire, pattern }
    }

    /// Check if this schedule should fire at `now`.
    /// Returns `true` if it fired (and advances the next_fire accordingly).
    pub fn check(&mut self, now: Time<UTC>) -> bool {
        if now < self.next_fire {
            return false;
        }

        // It fires now
        match self.pattern {
            Recurrence::Once => {
                // For Once, you might choose to never fire again;
                // leave next_fire as-is or set to max.
                self.pattern = Recurrence::Once;
            }
            Recurrence::Interval(period) => {
                // Catch up if we missed multiple ticks.
                while self.next_fire <= now {
                    self.next_fire = self.next_fire + period;
                }
            }
            Recurrence::Hourly { minute, second } => {
                self.next_fire = next_hourly(now, minute, second);
            }
            Recurrence::Daily { hour, minute, second } => {
                self.next_fire = next_daily(now, hour, minute, second);
            }
            Recurrence::Yearly { month, day, hour, minute, second } => {
                self.next_fire = next_yearly(now, month, day, hour, minute, second);
            }
        }

        true
    }

    pub fn next_fire(&self) -> Time<UTC> {
        self.next_fire
    }
}

impl Time<UTC> {
    pub fn to_utc_datetime(self) -> UtcDateTime {
        // Integer division / modulo: sec + nanos -> Y-M-D h:m:s.nnnnnnnnn
        let mut secs = self.sec;
        let second = (secs % 60) as u8;
        secs /= 60;
        let minute = (secs % 60) as u8;
        secs /= 60;
        let hour = (secs % 24) as u8;
        let days = secs / 24;

        let mut dt = days_to_ymd(days, hour, minute, second);
        dt.nanos = self.nano_sec; // preserve nanos
        dt
    }

    pub fn from_utc_datetime(dt: UtcDateTime) -> Time<UTC> {
        // Convert Y/M/D to days since Unix epoch
        let year = dt.year as i64;
        let month = dt.month as i64;
        let day = dt.day as i64;

        // Shift months so March = 1, Feb = 12 of previous year
        let (y, m) = if month <= 2 {
            (year - 1, month + 12)
        } else {
            (year, month)
        };

        // Days since civil epoch 0000-03-01
        let era = y / 400;
        let yoe = y - era * 400;                         // [0, 399]
        let doy = (153 * (m - 3) + 2) / 5 + day - 1;     // [0, 365]
        let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy; // [0, 146096]

        // Shift back to Unix epoch
        let days = era * 146_097 + doe - 719_468;

        // Convert to seconds since Unix epoch
        let secs =
            days * 24 * 3600 +
            (dt.hour as i64) * 3600 +
            (dt.minute as i64) * 60 +
            (dt.second as i64);

        // Carry nanos straight into Time<UTC>
        Time::<UTC>::new(secs, dt.nanos)
    }
}
