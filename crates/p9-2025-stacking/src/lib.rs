//! Reproduction of Geringer-Sameth, Golovich & Iwabuchi (2025),
//! "Multi-year stacking searches for solar system bodies" (arXiv:2509.25428).
//!
//! The paper extends digital tracking ("shift-and-stack") to full-sky,
//! multi-year surveys. Three pieces of its method are reproduced here from
//! first principles:
//!
//! - [`p9_core::analysis::stacking::matched_filter`] — the matched-filter stack: co-adding N background-
//!   limited frames along the correct track grows the signal-to-noise as √N
//!   (limiting magnitude as 1.25·log₁₀N), turning a ~20.5-mag single ZTF
//!   exposure into a mid-20s / ~27-mag effective depth.
//! - [`p9_core::analysis::stacking::orbit_metric`] — the heart of the paper: a trial orbit that misses the
//!   true track loses SNR *quadratically* in the parameter error, defining a
//!   Fisher **metric** on orbit space. The metric sets the grid spacing (and
//!   hence the **number of trial orbits**) needed to keep the SNR loss below a
//!   chosen tolerance. A longer baseline sharpens the metric (finer rate
//!   resolution) — which is exactly why a multi-year search needs *billions* of
//!   trial orbits.
//! - [`significance`] — the look-elsewhere penalty: searching N_trials orbits
//!   raises the per-trial detection threshold to ≈ √(2 ln N_trials) σ. The
//!   payoff/cost asymmetry is the headline: the stacking gain is several
//!   magnitudes while the billions-of-trials penalty costs only a few tenths of
//!   a magnitude.
//!
//! All computed quantities are pinned in seeded tests; the paper's headline
//! claims (ZTF: thousands of images over six years, ~10⁹ trial orbits, ~27th
//! magnitude reach) are kept as labelled reference constants.

pub mod significance;

/// Labelled reference numbers from the paper (for regression comparison only;
/// every quantity used in tests is computed, not taken from here).
pub mod published {
    /// ZTF single-exposure r-band limiting magnitude (Bellm et al. 2019).
    pub const ZTF_SINGLE_DEPTH: f64 = 20.5;
    /// Headline stacked depth reached for TNOs (~27th magnitude).
    pub const STACKED_DEPTH_REACH: f64 = 27.0;
    /// Order of magnitude of trial orbits in the multi-year ZTF search.
    pub const TRIAL_ORBITS_ORDER: f64 = 1.0e9;
    /// Survey baseline of the ZTF validation (six years).
    pub const ZTF_BASELINE_YEARS: f64 = 6.0;
}
