use super::duration::Duration;
use super::time::Time;
use aether_core::utils::lerp_i128;

#[derive(Clone, Debug)]
pub struct SyncSample<TS> {
    pub local_mono_ns: u64,
    pub remote_time: Time<TS>,
}

#[derive(Clone, Debug)]
pub struct RealtimeCoordinator<TS> {
    prev: Option<SyncSample<TS>>,
    last: Option<SyncSample<TS>>,
    max_extrapolation_ns: u64,
}

impl<TS> RealtimeCoordinator<TS> {
    pub const fn new(max_extrapolation_ns: u64) -> Self {
        Self {
            prev: None,
            last: None,
            max_extrapolation_ns,
        }
    }

    pub fn observe(&mut self, local_mono_ns: u64, remote_time: Time<TS>) {
        let sample = SyncSample {
            local_mono_ns,
            remote_time,
        };
        self.prev = self.last.as_ref().map(|s| SyncSample {
            local_mono_ns: s.local_mono_ns,
            remote_time: Time::new(s.remote_time.sec, s.remote_time.nano_sec),
        });
        self.last = Some(sample);
    }

    pub fn latest(&self) -> Option<SyncSample<TS>> {
        self.last.as_ref().map(|s| SyncSample {
            local_mono_ns: s.local_mono_ns,
            remote_time: Time::new(s.remote_time.sec, s.remote_time.nano_sec),
        })
    }

    pub fn estimate(&self, local_mono_ns: u64) -> Option<Time<TS>> {
        let last = self.last.as_ref()?;

        let dt_since_last = local_mono_ns.abs_diff(last.local_mono_ns);
        if dt_since_last > self.max_extrapolation_ns {
            return None;
        }

        let last_remote_ns = time_to_total_ns(&last.remote_time);

        match self.prev.as_ref() {
            None => {
                let delta_local = local_mono_ns as i128 - last.local_mono_ns as i128;
                let est_remote_ns = last_remote_ns + delta_local;
                Some(total_ns_to_time(est_remote_ns))
            }
            Some(prev) => {
                let prev_local = prev.local_mono_ns as i128;
                let last_local = last.local_mono_ns as i128;
                let prev_remote = time_to_total_ns(&prev.remote_time);
                let est_remote_ns = lerp_i128(
                    prev_local,
                    prev_remote,
                    last_local,
                    last_remote_ns,
                    local_mono_ns as i128,
                );
                Some(total_ns_to_time(est_remote_ns))
            }
        }
    }

    pub fn estimate_with_offset(
        &self,
        local_mono_ns: u64,
        fixed_offset: Duration,
    ) -> Option<Time<TS>> {
        let t = self.estimate(local_mono_ns)?;
        Some(t + fixed_offset)
    }
}

fn time_to_total_ns<TS>(t: &Time<TS>) -> i128 {
    (t.sec as i128) * 1_000_000_000i128 + (t.nano_sec as i128)
}

fn total_ns_to_time<TS>(total_ns: i128) -> Time<TS> {
    let sec = total_ns.div_euclid(1_000_000_000i128) as i64;
    let nano = total_ns.rem_euclid(1_000_000_000i128) as i32;
    Time::new(sec, nano)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time::TAI;

    #[test]
    fn interpolates_midpoint() {
        let mut c = RealtimeCoordinator::<TAI>::new(10_000_000_000);
        c.observe(1_000, Time::new(100, 0));
        c.observe(2_000, Time::new(100, 1_000_000));

        let t = c.estimate(1_500).unwrap();
        assert_eq!(t.sec, 100);
        assert_eq!(t.nano_sec, 500_000);
    }

    #[test]
    fn rejects_far_extrapolation() {
        let mut c = RealtimeCoordinator::<TAI>::new(100);
        c.observe(1_000, Time::new(100, 0));

        assert!(c.estimate(2_000).is_none());
    }
}
