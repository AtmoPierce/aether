#![no_std]

use core::marker::PhantomData;
use core::ops::{Add, Sub};
use core::sync::atomic::{AtomicI32, Ordering as AtomicOrdering};
use core::cmp::Ordering;


#[cfg(feature = "telemetry-defmt")]
use defmt::{Format};

// -------------------------
// Time scale markers
// -------------------------

#[derive(Copy, Clone, Debug)]
pub enum Unix {} // Seconds since 1970-01-01 UTC
#[derive(Copy, Clone, Debug)]
pub enum UTC  {} // Same epoch as Unix
#[derive(Copy, Clone, Debug)]
pub enum TAI  {} // TAI since Unix epoch
#[derive(Copy, Clone, Debug)]
pub enum GPS  {} // GPS seconds since GPS epoch

// -------------------------
// Leap second model
// -------------------------

/// The historical fixed offset at the Unix epoch:
/// At 1970-01-01: TAI - UTC = 10s
pub const TAI_UTC_BASE_OFFSET: i32 = 10;

/// Runtime leap second count (UTC lags TAI by BASE + leap_seconds)
static LEAP_SECONDS: AtomicI32 = AtomicI32::new(27); // 27 after 1972 -> 10+27 = 37s today

/// Current (TAI - UTC)
#[inline]
pub fn tai_minus_utc() -> i32 {
    TAI_UTC_BASE_OFFSET + LEAP_SECONDS.load(AtomicOrdering::Relaxed)
}

/// Set leap seconds at runtime (e.g., via GPS message)
#[inline]
pub fn set_leap_seconds(leaps: i32) {
    LEAP_SECONDS.store(leaps, AtomicOrdering::Relaxed);
}

// -------------------------
// GPS epoch
// -------------------------

/// Unix seconds at GPS epoch (1980-01-06 00:00:00 UTC)
pub const UNIX_GPS_EPOCH: i64 = 315_964_800;

// -------------------------
// Core Types
// -------------------------

#[derive(Clone, Copy, Debug)]
pub struct Time<TS> {
    pub sec: i64,
    pub nano_sec: i32,
    _ts: PhantomData<TS>,
}

impl<TS> Time<TS> {
    pub const fn new(sec: i64, nano_sec: i32) -> Self {
        Self { sec, nano_sec, _ts: PhantomData }
    }
}

// Equality: compare only the numeric fields, ignore the PhantomData.
impl<TS> PartialEq for Time<TS> {
    fn eq(&self, other: &Self) -> bool {
        self.sec == other.sec && self.nano_sec == other.nano_sec
    }
}
impl<TS> Eq for Time<TS> {}

// Ordering: same idea, only sec / nano_sec.
impl<TS> PartialOrd for Time<TS> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<TS> Ord for Time<TS> {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.sec.cmp(&other.sec) {
            Ordering::Equal => self.nano_sec.cmp(&other.nano_sec),
            ord => ord,
        }
    }
}

// -------------------------
// Unix <-> UTC
// -------------------------

impl Time<Unix> {
    pub fn to_utc(self) -> Time<UTC> {
        Time::<UTC>::new(self.sec, self.nano_sec)
    }
}

impl Time<UTC> {
    pub fn to_unix(self) -> Time<Unix> {
        Time::<Unix>::new(self.sec, self.nano_sec)
    }
}

// -------------------------
// UTC <-> TAI
// -------------------------

impl Time<UTC> {
    pub fn to_tai(self) -> Time<TAI> {
        Time::<TAI>::new(self.sec + tai_minus_utc() as i64, self.nano_sec)
    }
}

impl Time<TAI> {
    pub fn to_utc(self) -> Time<UTC> {
        Time::<UTC>::new(self.sec - tai_minus_utc() as i64, self.nano_sec)
    }
}

// -------------------------
// UTC <-> GPS
// -------------------------

impl Time<UTC> {
    pub fn to_gps(self) -> Time<GPS> {
        let tai_minus_utc = tai_minus_utc() as i64;
        let gps_minus_utc = tai_minus_utc - 19; // always 19 at epoch
        Time::<GPS>::new(self.sec - UNIX_GPS_EPOCH + gps_minus_utc, self.nano_sec)
    }
}

impl Time<GPS> {
    pub fn to_utc(self) -> Time<UTC> {
        let tai_minus_utc = tai_minus_utc() as i64;
        let gps_minus_utc = tai_minus_utc - 19;
        Time::<UTC>::new(self.sec + UNIX_GPS_EPOCH - gps_minus_utc, self.nano_sec)
    }
}


// -------------------------
// GPS <-> TAI
// -------------------------

impl Time<GPS> {
    pub fn to_tai(self) -> Time<TAI> {
        self.to_utc().to_tai()
    }
}

impl Time<TAI> {
    pub fn to_gps(self) -> Time<GPS> {
        self.to_utc().to_gps()
    }
}

