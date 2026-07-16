//! Rubin watch — Layer 1 (MPC/SBDB distant-object poller) per
//! `rubin_watch/design/`. The binary (`main.rs`) wires network sources into
//! the pure sync/battery/report machinery here; everything in this library
//! is offline and fixture-testable.

pub mod battery;
pub mod classify;
pub mod ledger;
pub mod report;
pub mod sources;
pub mod sync;

/// UTC timestamps without a chrono dependency: civil-from-days (Howard
/// Hinnant's algorithm) over `SystemTime`.
pub mod clock {
    /// (ISO-8601 instant, YYYY-MM-DD date) for now.
    pub fn now_utc() -> (String, String) {
        let secs = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock after 1970")
            .as_secs();
        iso_from_unix(secs)
    }

    pub fn iso_from_unix(secs: u64) -> (String, String) {
        let days = (secs / 86_400) as i64;
        let (y, m, d) = civil_from_days(days);
        let rem = secs % 86_400;
        let (hh, mm, ss) = (rem / 3600, (rem % 3600) / 60, rem % 60);
        (
            format!("{y:04}-{m:02}-{d:02}T{hh:02}:{mm:02}:{ss:02}Z"),
            format!("{y:04}-{m:02}-{d:02}"),
        )
    }

    fn civil_from_days(z: i64) -> (i64, u32, u32) {
        let z = z + 719_468;
        let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
        let doe = (z - era * 146_097) as u64;
        let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
        let y = yoe as i64 + era * 400;
        let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
        let mp = (5 * doy + 2) / 153;
        let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
        let m = (if mp < 10 { mp + 3 } else { mp - 9 }) as u32;
        (if m <= 2 { y + 1 } else { y }, m, d)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn epoch_and_known_dates() {
            assert_eq!(iso_from_unix(0).0, "1970-01-01T00:00:00Z");
            // 2026-07-16 00:00:00 UTC = 1784160000
            assert_eq!(iso_from_unix(1_784_160_000).1, "2026-07-16");
            // Leap-year day: 2024-02-29 12:00:00 UTC = 1709208000
            assert_eq!(iso_from_unix(1_709_208_000).0, "2024-02-29T12:00:00Z");
        }
    }
}
