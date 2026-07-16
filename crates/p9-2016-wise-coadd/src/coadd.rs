//! Coadd depth gain from √N frame stacking.
//!
//! Meisner, Bromley, Kenyon & Anderson (2016) build *inertial coadds* of W1
//! exposures binned into ∼1-day intervals, then shift-and-stack those daily
//! coadds along trial Planet Nine motions. Stacking N independent,
//! background-limited frames raises the signal-to-noise by √N (signal adds as
//! N, noise as √N), so the limiting flux drops by √N and the limiting
//! *magnitude* deepens by
//!
//!   Δm = 2.5 · log10(√N) = 1.25 · log10(N).
//!
//! This is the same matched-filter / coadd √N law used by the
//! `p9-2025-stacking` crate; here it converts the shared single-frame WISE
//! depth (`p9_core::analysis::surveys`, "WISE W1" = 16.5) into the deeper
//! coadd depth that lets the search reach ~800 AU instead of the
//! single-exposure ~430 AU quoted in the abstract.
//!
//! Per-band depth bookkeeping (Vega), disentangled from a previous version
//! that fed the shared table's 16.5 — itself a single-COVERAGE STACK depth
//! (~8 L1b exposures per sky pass) — into the √N law as if it were one
//! exposure, making `coadd_depth(W1, 100) = 19.0` a depth no WISE coadd
//! achieves and "explaining" the published 16.66 with a nonsensical
//! N ≈ 1.3:
//!
//! - single L1b exposure, 5σ:  W1 ≈ 15.3, W2 ≈ 14.5 (Wright et al. 2010)
//! - one coverage (~8 exp):    W1 ≈ 15.3 + 1.25·log₁₀8 ≈ 16.4 — the shared
//!   table's 16.5 "single-frame-stack" reference
//! - Meisner et al. (2016) coadd 90% completeness: W1 ≈ 16.66, i.e.
//!   N ≈ 12 exposures per ~1-day bin from the single-exposure base — the
//!   physically sensible epoch depth.

use crate::band::{WiseBand, W2};

/// Single L1b EXPOSURE W1 depth (Vega, 5σ; Wright et al. 2010). The base of
/// the √N coadd law. (The shared survey table's 16.5 is the ~8-exposure
/// single-coverage stack, NOT one frame — see the module docs.)
pub const SINGLE_FRAME_W1_DEPTH: f64 = 15.3;

/// The shared survey table's W1 = 16.5 single-coverage stack reference,
/// consistent with ≈8 exposures over the single-exposure base.
pub const STACK_REFERENCE_W1_DEPTH: f64 = 16.5;

/// Single L1b exposure W2 depth (Vega, 5σ). W2 is intrinsically noisier; the
/// familiar 15.6 figure is its ~8-exposure coverage depth, putting one
/// exposure near 14.5.
pub const SINGLE_FRAME_W2_DEPTH: f64 = 14.5;

/// Published 2016 coadd W1 completeness depth (Vega): Meisner et al. (2016)
/// quote 90% completeness at W1 < 16.66 over ∼2000 deg². Labelled reference.
pub const PUBLISHED_COADD_W1_DEPTH: f64 = 16.66;

/// Single-frame depth (Vega) for a WISE band.
pub fn single_frame_depth(band: WiseBand) -> f64 {
    if band == W2 {
        SINGLE_FRAME_W2_DEPTH
    } else {
        SINGLE_FRAME_W1_DEPTH
    }
}

/// Magnitude deepening from coadding `n_frames` background-limited frames:
/// Δm = 2.5·log10(√N) = 1.25·log10(N). Zero for N ≤ 1.
pub fn coadd_depth_gain(n_frames: u32) -> f64 {
    if n_frames <= 1 {
        return 0.0;
    }
    1.25 * (n_frames as f64).log10()
}

/// Coadd limiting magnitude in `band` after stacking `n_frames`: the
/// single-frame depth plus the √N gain (deeper = larger magnitude).
pub fn coadd_depth(band: WiseBand, n_frames: u32) -> f64 {
    single_frame_depth(band) + coadd_depth_gain(n_frames)
}

/// Number of frames implied by a desired magnitude `gain` under the √N law:
/// N = 10^(gain / 1.25). Useful to invert the published coadd depth.
pub fn frames_for_gain(gain: f64) -> f64 {
    10.0_f64.powf(gain / 1.25)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::band::W1;
    use approx::assert_relative_eq;

    #[test]
    fn coadd_is_deeper_than_single_frame() {
        // Any real stack (N>1) is strictly deeper than the single-frame depth.
        let n = 100;
        let deep = coadd_depth(W1, n);
        assert!(
            deep > SINGLE_FRAME_W1_DEPTH,
            "coadd {deep:.3} must exceed single-frame {SINGLE_FRAME_W1_DEPTH}"
        );
    }

    #[test]
    fn gain_follows_root_n_law() {
        // Δm(N) = 1.25 log10 N; quadrupling N adds 1.25·log10(4) ≈ 0.753 mag.
        assert_relative_eq!(
            coadd_depth_gain(4),
            1.25 * 4.0_f64.log10(),
            max_relative = 1e-12
        );
        // 100 frames -> 1.25·2 = 2.5 mag deeper.
        assert_relative_eq!(coadd_depth_gain(100), 2.5, max_relative = 1e-12);
        // Stacking is a √N flux improvement: 4× frames -> 2× SNR -> 0.753 mag.
        assert!((coadd_depth_gain(4) - 0.7526).abs() < 1e-3);
    }

    #[test]
    fn gain_is_monotone_in_n() {
        let mut last = -1.0;
        for n in [1u32, 2, 10, 100, 1000] {
            let g = coadd_depth_gain(n);
            assert!(g >= last, "gain must be non-decreasing in N");
            last = g;
        }
        assert_eq!(coadd_depth_gain(1), 0.0);
        assert_eq!(coadd_depth_gain(0), 0.0);
    }

    #[test]
    fn frames_reproduce_published_coadd_depth() {
        // The published W1 coadd 90% depth (16.66) sits 1.36 mag over the
        // single-EXPOSURE base: N ≈ 12 exposures per ~1-day bin — the
        // physically sensible epoch depth (the old 16.5-as-single-frame
        // bookkeeping "explained" it with N ≈ 1.3).
        let gain = PUBLISHED_COADD_W1_DEPTH - SINGLE_FRAME_W1_DEPTH;
        let n = frames_for_gain(gain);
        assert!(
            (8.0..20.0).contains(&n),
            "gain {gain:.3} -> N {n:.1} exposures/epoch"
        );
        // And the shared table's 16.5 single-coverage reference is the
        // ~8-exposure stack over the same base.
        let n_cov = frames_for_gain(STACK_REFERENCE_W1_DEPTH - SINGLE_FRAME_W1_DEPTH);
        assert!((5.0..13.0).contains(&n_cov), "coverage N = {n_cov:.1}");
    }

    #[test]
    fn w2_single_frame_is_shallower_than_w1() {
        assert!(single_frame_depth(W2) < single_frame_depth(W1));
    }
}
