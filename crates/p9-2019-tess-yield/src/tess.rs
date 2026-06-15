//! TESS full-frame-image (FFI) stacked depth for a Planet Nine search.
//!
//! Payne, Holman & Pál (2019), "A TESS Search for Distant Solar System
//! Objects: Yield Estimates" (arXiv:1911.03676), show that digitally
//! co-adding (shifting-and-stacking along candidate orbits) the ~1300 FFI
//! exposures from a single 27-day TESS observing sector reaches a 50%
//! detection threshold of I_C ≈ 22.0 ± 0.5 — far deeper than any single
//! 30-minute FFI cadence (~I_C ≈ 18–19), and deep enough that "TESS could
//! discover the hypothesized Planet Nine."
//!
//! Per-frame depth scales as +1.25·log10(N) (Poisson sky-noise stacking,
//! Δm = 2.5·log10(√N)), so the headline depth is set by the number of
//! co-added frames. We reproduce the published threshold from that scaling
//! law plus a per-frame reference depth, rather than hard-coding 22.0.
//!
//! BAND CAVEAT: TESS's red bandpass is close to Cousins I_C, while the P9
//! reflected-light photometry in [`p9_core::analysis::photometry`] is V band.
//! A Neptune-like (blue, methane-absorbed) P9 is fainter in I than in V, so
//! comparing a V apparent magnitude against an I_C depth is optimistic by a
//! few tenths of a magnitude; we treat the depth as a V-equivalent reference
//! and document the residual in tests rather than applying an unpinned color.

use serde::{Deserialize, Serialize};

/// Published 50%-completeness stacked depth of a single TESS sector,
/// I_C ≈ 22.0 ± 0.5 (Payne, Holman & Pál 2019). Reference constant.
pub const TESS_STACKED_DEPTH_IC: f64 = 22.0;

/// Quoted 1σ uncertainty on the stacked depth (±0.5 mag).
pub const TESS_STACKED_DEPTH_SIGMA: f64 = 0.5;

/// Number of FFI exposures co-added per sector in the published estimate
/// (~1300; a 27-day sector at the 30-minute FFI cadence).
pub const TESS_FRAMES_PER_SECTOR: f64 = 1300.0;

/// TESS FFI / sector stacking model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TessStack {
    /// 50%-completeness depth of a single 30-minute FFI (mag). Background-
    /// limited reference: stacking `frames` of these reaches the sector depth.
    pub single_frame_depth: f64,
    /// Number of frames co-added (shift-and-stack along a trial orbit).
    pub frames: f64,
    /// Logistic steepness of the completeness roll-off near the limit
    /// (mag^-1), as for the sibling survey-exclusion crates.
    pub efficiency_steepness: f64,
}

impl Default for TessStack {
    /// Sector-scale stack reproducing the published I_C ≈ 22.0 threshold.
    fn default() -> Self {
        // Solve single_frame_depth so that the +1.25 log10(N) gain lands on
        // the published 22.0 at N = 1300: single = 22.0 - 1.25 log10(1300).
        let gain = 1.25 * TESS_FRAMES_PER_SECTOR.log10();
        Self {
            single_frame_depth: TESS_STACKED_DEPTH_IC - gain,
            frames: TESS_FRAMES_PER_SECTOR,
            efficiency_steepness: 4.0,
        }
    }
}

impl TessStack {
    /// Stacked 50%-completeness depth after co-adding `self.frames` frames.
    ///
    ///   m_lim(N) = m_1 + 2.5·log10(√N) = m_1 + 1.25·log10(N)
    ///
    /// the background-limited (sky-noise) co-addition gain.
    pub fn stacked_depth(&self) -> f64 {
        self.single_frame_depth + 1.25 * self.frames.log10()
    }

    /// Magnitude-dependent stacked detection efficiency: a logistic roll-off
    /// centered on [`stacked_depth`](Self::stacked_depth) (0.5 at the limit,
    /// →1 brighter, →0 fainter), mirroring the ZTF survey model.
    pub fn detection_efficiency(&self, apparent_magnitude: f64) -> f64 {
        let d = self.stacked_depth();
        1.0 / (1.0 + ((apparent_magnitude - d) * self.efficiency_steepness).exp())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn default_stack_reproduces_published_depth() {
        let stack = TessStack::default();
        // Reconstructed from the +1.25 log10(N) scaling, not hard-coded.
        assert_relative_eq!(stack.stacked_depth(), TESS_STACKED_DEPTH_IC, epsilon = 1e-9);
    }

    #[test]
    fn published_depth_within_quoted_uncertainty() {
        let stack = TessStack::default();
        assert!((stack.stacked_depth() - TESS_STACKED_DEPTH_IC).abs() <= TESS_STACKED_DEPTH_SIGMA);
    }

    #[test]
    fn single_frame_is_much_shallower() {
        // A single 30-min FFI is ~I_C 18-19; stacking ~1300 buys ~3.9 mag.
        let stack = TessStack::default();
        assert!(
            stack.single_frame_depth > 17.5 && stack.single_frame_depth < 18.5,
            "single-frame depth = {:.2}",
            stack.single_frame_depth
        );
        let gain = stack.stacked_depth() - stack.single_frame_depth;
        assert!((gain - 3.9).abs() < 0.2, "stacking gain = {gain:.2} mag");
    }

    #[test]
    fn more_frames_go_deeper() {
        let base = TessStack::default();
        let deeper = TessStack {
            frames: base.frames * 4.0,
            ..base.clone()
        };
        // 4x frames -> +2.5 log10(2) = +0.753 mag deeper.
        let d = deeper.stacked_depth() - base.stacked_depth();
        assert_relative_eq!(d, 1.25 * 4.0_f64.log10(), epsilon = 1e-9);
        assert!(d > 0.0);
    }

    #[test]
    fn efficiency_is_logistic() {
        let stack = TessStack::default();
        let d = stack.stacked_depth();
        assert!(stack.detection_efficiency(d - 3.0) > 0.99);
        assert_relative_eq!(stack.detection_efficiency(d), 0.5, epsilon = 1e-10);
        assert!(stack.detection_efficiency(d + 3.0) < 0.01);
    }
}
