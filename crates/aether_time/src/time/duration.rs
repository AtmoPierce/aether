use core::cmp::Ordering;
use core::ops::{Add, Sub};

use super::time::Time;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Duration {
    pub sec: i64,
    pub nano_sec: i32,
}

// -------------------------
// Duration utilities
// -------------------------

impl Duration {
    pub const fn from_nanos(total: i64) -> Self {
        let sec = total / 1_000_000_000;
        let nano_sec = (total % 1_000_000_000) as i32;
        Self { sec, nano_sec }
    }

    pub fn normalize(mut self) -> Self {
        if self.nano_sec >= 1_000_000_000 {
            self.sec += (self.nano_sec / 1_000_000_000) as i64;
            self.nano_sec %= 1_000_000_000;
        } else if self.nano_sec < 0 {
            self.sec -= 1;
            self.nano_sec += 1_000_000_000;
        }
        self
    }

    pub const fn zero() -> Self {
        Duration { sec: 0, nano_sec: 0 }
    }

    pub const fn from_secs(sec: i64) -> Self {
        Duration { sec, nano_sec: 0 }
    }

    pub fn from_millis(ms: i64) -> Self {
        Duration {
            sec: ms / 1000,
            nano_sec: ((ms % 1000) * 1_000_000) as i32,
        }.normalize()
    }

    pub fn from_micros(us: i64) -> Self {
        Duration {
            sec: us / 1_000_000,
            nano_sec: ((us % 1_000_000) * 1_000) as i32,
        }.normalize()
    }

    pub const fn from_nanos_i64(ns: i64) -> Self {
        Duration::from_nanos(ns)
    }
}

// -------------------------
// Arithmetic
// -------------------------

impl<TS> Add<Duration> for Time<TS> {
    type Output = Time<TS>;

    fn add(self, rhs: Duration) -> Self::Output {
        let mut sec = self.sec + rhs.sec;
        let mut nano_sec = self.nano_sec + rhs.nano_sec;
        if nano_sec >= 1_000_000_000 {
            sec += 1;
            nano_sec -= 1_000_000_000;
        } else if nano_sec < 0 {
            sec -= 1;
            nano_sec += 1_000_000_000;
        }
        Time::new(sec, nano_sec)
    }
}

impl<TS> Sub<Time<TS>> for Time<TS> {
    type Output = Duration;

    fn sub(self, rhs: Time<TS>) -> Self::Output {
        let mut sec = self.sec - rhs.sec;
        let mut nano_sec = self.nano_sec - rhs.nano_sec;
        if nano_sec < 0 {
            sec -= 1;
            nano_sec += 1_000_000_000;
        }
        Duration { sec, nano_sec }
    }
}

// References
impl<TS> core::ops::Sub<&Time<TS>> for &Time<TS> {
    type Output = Duration;

    fn sub(self, rhs: &Time<TS>) -> Duration {
        let mut sec = self.sec - rhs.sec;
        let mut nano_sec = self.nano_sec - rhs.nano_sec;
        if nano_sec < 0 {
            sec -= 1;
            nano_sec += 1_000_000_000;
        }
        Duration { sec, nano_sec }
    }
}

impl<TS> core::ops::Add<Duration> for &Time<TS> {
    type Output = Time<TS>;

    fn add(self, rhs: Duration) -> Time<TS> {
        let mut sec = self.sec + rhs.sec;
        let mut nano_sec = self.nano_sec + rhs.nano_sec;
        if nano_sec >= 1_000_000_000 {
            sec += 1;
            nano_sec -= 1_000_000_000;
        } else if nano_sec < 0 {
            sec -= 1;
            nano_sec += 1_000_000_000;
        }
        Time::new(sec, nano_sec)
    }
}

impl core::ops::Add for Duration {
    type Output = Duration;

    fn add(self, rhs: Duration) -> Duration {
        Duration {
            sec: self.sec + rhs.sec,
            nano_sec: self.nano_sec + rhs.nano_sec,
        }.normalize()
    }
}

impl core::ops::Sub for Duration {
    type Output = Duration;

    fn sub(self, rhs: Duration) -> Duration {
        Duration {
            sec: self.sec - rhs.sec,
            nano_sec: self.nano_sec - rhs.nano_sec,
        }.normalize()
    }
}

impl PartialOrd for Duration {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Duration {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.sec.cmp(&other.sec) {
            Ordering::Equal => self.nano_sec.cmp(&other.nano_sec),
            ord => ord,
        }
    }
}