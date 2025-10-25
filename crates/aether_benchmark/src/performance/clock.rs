// src/performance/clock.rs
use minstant::{Anchor, Instant};
use once_cell::sync::Lazy;

pub static ANCHOR: Lazy<Anchor> = Lazy::new(Anchor::new);

#[inline]
pub fn now() -> Instant {
    Instant::now()
}
#[inline]
pub fn unix_nanos_now() -> u64 {
    Instant::now().as_unix_nanos(&*ANCHOR)
}
