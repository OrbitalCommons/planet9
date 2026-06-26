//! The Q1–Q6 eligibility funnel and its published bin counts.
//!
//! Eldadi & Loeb (2026) build a candidate "bin" for every (KBO × observatory ×
//! photometric band) combination in the MPC archive and pass each through a
//! six-criterion eligibility pipeline. Q1–Q3 cull bins that cannot support a
//! slope fit at all (too few epochs, too small a heliocentric-distance lever
//! arm, unusable photometry); Q4–Q6 are the stricter quality cuts that leave a
//! clean fit. The survivors are then classified by their fitted α.
//!
//! The published funnel (reference constants — the counts are the paper's, not
//! recomputed here):
//!
//! ```text
//!   8557  candidate bins (KBO × observatory × band)
//!     │  Q1–Q3
//!   1089  pass
//!     │  Q4–Q6
//!    186  additionally pass
//!     ├── 53   reflected sunlight   (α ≈ −4, natural)
//!     ├── 24   self-luminous        (α ≈ −2, anomalous)
//!     └── 109  anomalous / other
//! ```
//!
//! Two diagnostic footnotes are also encoded:
//!
//! - Pluto, the brightest TNO, has 22 (observatory × band) bins, and *none* of
//!   them recovers the reflected α ≈ −4 when restricted to a single
//!   instrument+band — the archive can't cleanly run the test even on the
//!   easiest target.
//! - All 24 "self-luminous" detections come from Pan-STARRS PS1/PS2 alone,
//!   pointing at calibration systematics rather than real self-luminosity.

/// Total candidate bins (KBO × observatory × band) in the MPC archive.
pub const CANDIDATE_BINS: u32 = 8557;

/// Bins passing the Q1–Q3 cuts.
pub const PASS_Q1_Q3: u32 = 1089;

/// Bins additionally passing the Q4–Q6 cuts (the clean-fit sample).
pub const PASS_Q4_Q6: u32 = 186;

/// Of [`PASS_Q4_Q6`], bins consistent with reflected sunlight (α ≈ −4).
pub const REFLECTED: u32 = 53;

/// Of [`PASS_Q4_Q6`], bins consistent with self-luminosity (α ≈ −2).
pub const SELF_LUMINOUS: u32 = 24;

/// Of [`PASS_Q4_Q6`], anomalous / other bins (α near neither −4 nor −2).
pub const ANOMALOUS: u32 = 109;

/// Pluto's (observatory × band) bins.
pub const PLUTO_BINS: u32 = 22;

/// Pluto bins that recover the reflected α ≈ −4 under a single instrument+band.
pub const PLUTO_REFLECTED_RECOVERIES: u32 = 0;

/// Whether all 24 self-luminous detections originate solely from Pan-STARRS
/// PS1/PS2 (interpreted as a calibration systematic, not real self-luminosity).
pub const SELF_LUMINOUS_ALL_PANSTARRS: bool = true;

/// The published Q1–Q6 funnel as a single value, with the consistency check in
/// [`Funnel::is_consistent`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Funnel {
    /// Candidate bins (KBO × observatory × band).
    pub candidate_bins: u32,
    /// Bins passing Q1–Q3.
    pub pass_q1_q3: u32,
    /// Bins additionally passing Q4–Q6.
    pub pass_q4_q6: u32,
    /// Reflected-sunlight detections (α ≈ −4).
    pub reflected: u32,
    /// Self-luminous detections (α ≈ −2).
    pub self_luminous: u32,
    /// Anomalous / other detections.
    pub anomalous: u32,
}

/// The published funnel of Eldadi & Loeb (2026).
pub const PUBLISHED_FUNNEL: Funnel = Funnel {
    candidate_bins: CANDIDATE_BINS,
    pass_q1_q3: PASS_Q1_Q3,
    pass_q4_q6: PASS_Q4_Q6,
    reflected: REFLECTED,
    self_luminous: SELF_LUMINOUS,
    anomalous: ANOMALOUS,
};

impl Funnel {
    /// Internal-consistency check on the funnel:
    ///
    /// - the three α-classes partition the Q4–Q6 survivors
    ///   (reflected + self-luminous + anomalous == pass_q4_q6), and
    /// - each cut is a strict subset of the previous stage
    ///   (pass_q4_q6 ≤ pass_q1_q3 ≤ candidate_bins).
    pub fn is_consistent(&self) -> bool {
        self.reflected + self.self_luminous + self.anomalous == self.pass_q4_q6
            && self.pass_q4_q6 <= self.pass_q1_q3
            && self.pass_q1_q3 <= self.candidate_bins
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn published_funnel_is_internally_consistent() {
        assert!(PUBLISHED_FUNNEL.is_consistent());
    }

    #[test]
    fn classes_partition_the_survivors() {
        assert_eq!(REFLECTED + SELF_LUMINOUS + ANOMALOUS, PASS_Q4_Q6);
        assert_eq!(53 + 24 + 109, 186);
    }

    #[test]
    fn cuts_are_nested() {
        assert!(PASS_Q4_Q6 <= PASS_Q1_Q3);
        assert!(PASS_Q1_Q3 <= CANDIDATE_BINS);
    }

    #[test]
    fn pluto_recovers_no_reflected_slope() {
        assert_eq!(PLUTO_BINS, 22);
        assert_eq!(PLUTO_REFLECTED_RECOVERIES, 0);
    }

    #[test]
    fn self_luminous_are_all_panstarrs() {
        assert!(SELF_LUMINOUS_ALL_PANSTARRS);
    }

    #[test]
    fn inconsistent_funnel_is_rejected() {
        let bad = Funnel {
            reflected: 999,
            ..PUBLISHED_FUNNEL
        };
        assert!(!bad.is_consistent());
    }
}
