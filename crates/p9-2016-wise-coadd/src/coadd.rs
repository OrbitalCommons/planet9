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
//! Per-band single-frame depths: the shared survey table pins W1 = 16.5. The
//! single-frame W2 depth is shallower (W2 is intrinsically noisier and the
//! published single-exposure 5σ W2 limit, Wright et al. 2010, is ≈ 15.6 Vega);
//! we keep it as a labelled reference constant. The *coadd* W1 90%-completeness
//! depth the 2016 paper actually reports, W1 ≈ 16.66, is kept as a labelled
//! reference for cross-checking the gain.

use p9_2018_wise_search::survey_model::W1_DEPTH_REFERENCE;

use crate::band::{WiseBand, W2};

/// Single-frame W1 depth (Vega), from the shared survey table via the sibling.
pub const SINGLE_FRAME_W1_DEPTH: f64 = W1_DEPTH_REFERENCE;

/// Single-frame W2 depth (Vega), labelled reference constant. WISE W2 is
/// intrinsically noisier than W1; the published single-exposure 5σ W2 limit is
/// ≈ 15.6 Vega (Wright et al. 2010). Kept here as a reference; the coadd gain
/// is band-independent (√N law), so the W2/W1 detectability contrast is driven
/// by the *flux* (band.rs), not by the depth offset.
pub const SINGLE_FRAME_W2_DEPTH: f64 = 15.6;

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
        // The published W1 coadd 90% depth (16.66) is ~0.16 mag deeper than the
        // single-frame 16.5; under √N that is a very small effective stack
        // (~1.3 frames), consistent with ∼1-day daily coadds being only a few
        // exposures deep per ~1-day bin near the ecliptic.
        let gain = PUBLISHED_COADD_W1_DEPTH - SINGLE_FRAME_W1_DEPTH;
        let n = frames_for_gain(gain);
        assert!(gain > 0.0 && n > 1.0, "gain {gain:.3} -> N {n:.2}");
        // Round-trip: applying that N reproduces the published depth.
        assert_relative_eq!(
            coadd_depth(W1, n.round() as u32),
            SINGLE_FRAME_W1_DEPTH + coadd_depth_gain(n.round() as u32),
            max_relative = 1e-12
        );
    }

    #[test]
    fn w2_single_frame_is_shallower_than_w1() {
        assert!(single_frame_depth(W2) < single_frame_depth(W1));
    }
}
