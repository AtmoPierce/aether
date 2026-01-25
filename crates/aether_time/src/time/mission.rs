use core::marker::PhantomData;
use super::time::{Time, TAI, UTC, GPS};
use super::duration::Duration;

// -------------------------
// Mission Epoch
// -------------------------

#[derive(Copy, Clone, Debug)]
pub struct MissionEpoch {
    start_tai: Time<TAI>,
}

impl MissionEpoch {
    pub const fn new_tai(start_tai: Time<TAI>) -> Self {
        Self { start_tai }
    }

    pub fn new_utc(start_utc: Time<UTC>) -> Self {
        Self { start_tai: start_utc.to_tai() }
    }

    pub fn start(&self) -> Time<TAI> {
        self.start_tai
    }
}

// -------------------------
// TimeKeeper
// -------------------------

#[derive(Copy, Clone, Debug)]
pub struct Epoch<TS> {
    mission: MissionEpoch,
    ticks0: u64,
    ticks_per_sec: u64,
    _ts: PhantomData<TS>,
}

impl<TS> Epoch<TS> {
    pub const fn new(mission: MissionEpoch, ticks0: u64, ticks_per_sec: u64) -> Self {
        Self { mission, ticks0, ticks_per_sec, _ts: PhantomData }
    }

    fn duration_from_ticks(&self, now: u64) -> Duration {
        let dt = now.wrapping_sub(self.ticks0);
        let sec = (dt / self.ticks_per_sec) as i64;
        let frac = (dt % self.ticks_per_sec) as u128;
        let nano = ((frac * 1_000_000_000u128) / (self.ticks_per_sec as u128)) as i32;
        Duration { sec, nano_sec: nano }
    }

    pub fn t_plus(&self, now: u64) -> Time<TS> {
        let dt = self.duration_from_ticks(now);
        Time::<TS>::new(dt.sec, dt.nano_sec)
    }

    pub fn mission(&self) -> MissionEpoch {
        self.mission
    }
}

#[derive(Copy, Clone, Debug)]
pub struct TimeKeeper<TS> {
    epoch: Epoch<TS>,
}

impl<TS> TimeKeeper<TS> {
    pub const fn new(mission: MissionEpoch, ticks0: u64, ticks_per_sec: u64) -> Self {
        Self {
            epoch: Epoch::new(mission, ticks0, ticks_per_sec),
        }
    }

    /// T+ time in this scale: seconds since mission epoch.
    pub fn t_plus(&self, ticks: u64) -> Time<TS> {
        self.epoch.t_plus(ticks)
    }

    /// Absolute TAI time "now".
    pub fn now_tai(&self, ticks: u64) -> Time<TAI> {
        let dt = self.epoch.duration_from_ticks(ticks);
        self.epoch.mission().start() + dt
    }

    /// Compute dt based on ticks, with the *caller* holding last_ticks.
    ///
    /// - `last_ticks`: your stored tick count from the previous loop.
    /// - `now_ticks`: current hardware tick count.
    /// Returns duration corresponding to (now_ticks - last_ticks),
    /// and updates `last_ticks` in-place.
    pub fn dt_from_ticks(&self, last_ticks: &mut u64, now_ticks: u64) -> Duration {
        let dt_ticks = now_ticks.wrapping_sub(*last_ticks);
        *last_ticks = now_ticks;

        let sec  = (dt_ticks / self.epoch.ticks_per_sec) as i64;
        let frac = (dt_ticks % self.epoch.ticks_per_sec) as u128;
        let nano = ((frac * 1_000_000_000u128) / (self.epoch.ticks_per_sec as u128)) as i32;

        Duration { sec, nano_sec: nano }
    }
}


impl TimeKeeper<UTC> {
    pub fn now_utc(&self, ticks: u64) -> Time<UTC> {
        self.now_tai(ticks).to_utc()
    }
}

impl TimeKeeper<GPS> {
    pub fn now_gps(&self, ticks: u64) -> Time<GPS> {
        self.now_tai(ticks).to_gps()
    }
}


#[derive(Copy, Clone, Debug)]
pub struct TimerEvent {
    acc: Duration,
    period: Duration,
}

impl TimerEvent {
    /// Create a periodic event with the given period.
    pub const fn new(period: Duration) -> Self {
        Self {
            acc: Duration::zero(),
            period,
        }
    }

    /// Step the timer with a dt.
    ///
    /// Returns `true` if the event fired *at least once* during this step.
    /// If you want to count how many periods elapsed, we can extend this.
    pub fn step(&mut self, dt: Duration) -> bool {
        
        self.acc = self.acc + dt;
        if self.acc >= self.period {
            // preserve overshoot – wrap but keep remainder
            while self.acc >= self.period {
                self.acc = self.acc - self.period;
            }
            true
        } else {
            false
        }
    }
}


#[derive(Debug)]
pub struct TimerSet<const N: usize> {
    events: [TimerEvent; N],
}

impl<const N: usize> TimerSet<N> {
    /// Construct from an array of periods.
    pub const fn new(periods: [Duration; N]) -> Self {
        // const-friendly construction
        let mut events: [TimerEvent; N] = [TimerEvent::new(Duration::zero()); N];
        let mut i = 0;
        while i < N {
            events[i] = TimerEvent::new(periods[i]);
            i += 1;
        }
        Self { events }
    }

    /// Advance all timers by `dt` and return which ones fired.
    pub fn step(&mut self, dt: Duration) -> [bool; N] {
        let mut fired = [false; N];
        let mut i = 0;
        while i < N {
            if self.events[i].step(dt) {
                fired[i] = true;
            }
            i += 1;
        }
        fired
    }

    /// Optional: mutable access to a specific event (for reset/retuning).
    pub fn event_mut(&mut self, idx: usize) -> Option<&mut TimerEvent> {
        self.events.get_mut(idx)
    }
}
