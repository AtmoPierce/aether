// src/performance/timer.rs
#![cfg(feature = "std")]
use crate::performance::clock::now;
use crate::performance::performance::TimingCsv;

#[must_use]
pub struct TimeSpan<'a> {
    name: &'static str,
    t0: minstant::Instant,
    sink: &'a TimingCsv,
}
impl<'a> TimeSpan<'a> {
    #[inline]
    pub fn new(name: &'static str, sink: &'a TimingCsv) -> Self {
        Self {
            name,
            t0: now(),
            sink,
        }
    }
}
impl<'a> Drop for TimeSpan<'a> {
    fn drop(&mut self) {
        let ms = (now() - self.t0).as_secs_f64() * 1e3;
        self.sink.timing_ms(self.name, ms);
    }
}
