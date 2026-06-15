//! Matched-filter stacking: SNR and limiting-magnitude gain from co-adding
//! frames along a track.
//!
//! For N background-limited exposures of equal depth, the optimal (matched-
//! filter) combination adds the signal coherently and the noise in quadrature,
//! so the stacked SNR is √N times a single frame. In flux terms the limiting
//! flux falls as 1/√N, i.e. the limiting magnitude deepens by
//! Δm = 2.5·log₁₀(√N) = 1.25·log₁₀(N).

/// Matched-filter stacked SNR for `n_frames` equal-depth exposures, each with
/// single-frame SNR `per_frame_snr` (background/read-noise limited).
pub fn stacked_snr(per_frame_snr: f64, n_frames: usize) -> f64 {
    per_frame_snr * (n_frames as f64).sqrt()
}

/// General matched-filter combination of per-frame SNRs (unequal frames):
/// SNR_stack = √(Σ SNRᵢ²).
pub fn combine_snr(per_frame_snrs: &[f64]) -> f64 {
    per_frame_snrs.iter().map(|s| s * s).sum::<f64>().sqrt()
}

/// Limiting-magnitude gain from stacking `n_frames` background-limited frames:
/// Δm = 1.25·log₁₀(N).
pub fn stack_depth_gain_mag(n_frames: usize) -> f64 {
    if n_frames == 0 {
        return f64::NEG_INFINITY;
    }
    1.25 * (n_frames as f64).log10()
}

/// Stacked limiting magnitude from a single-frame depth and a frame count.
pub fn stacked_limiting_mag(single_frame_depth: f64, n_frames: usize) -> f64 {
    single_frame_depth + stack_depth_gain_mag(n_frames)
}

/// Number of equal frames needed to reach a target stacked depth from a single-
/// frame depth: N = 10^((m_target − m_single)/1.25).
pub fn frames_to_reach_depth(single_frame_depth: f64, target_depth: f64) -> f64 {
    10f64.powf((target_depth - single_frame_depth) / 1.25)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::published;
    use approx::assert_relative_eq;

    #[test]
    fn stacking_grows_snr_as_sqrt_n() {
        assert_relative_eq!(stacked_snr(5.0, 100), 50.0, epsilon = 1e-9);
        // doubling frames raises SNR by √2.
        assert_relative_eq!(
            stacked_snr(5.0, 200) / stacked_snr(5.0, 100),
            std::f64::consts::SQRT_2,
            epsilon = 1e-9
        );
    }

    #[test]
    fn matched_filter_equals_quadrature_sum() {
        let snrs = [5.0; 100];
        assert_relative_eq!(combine_snr(&snrs), stacked_snr(5.0, 100), epsilon = 1e-9);
        // unequal frames still add in quadrature
        assert_relative_eq!(combine_snr(&[3.0, 4.0]), 5.0, epsilon = 1e-9);
    }

    #[test]
    fn depth_gain_is_1_25_log10_n() {
        assert_relative_eq!(stack_depth_gain_mag(10_000), 5.0, epsilon = 1e-9);
        assert_relative_eq!(stack_depth_gain_mag(100), 2.5, epsilon = 1e-9);
    }

    #[test]
    fn ztf_stack_reaches_27th_magnitude() {
        // The paper reaches ~27th mag for TNOs; from a 20.5-mag single ZTF
        // exposure that needs N ≈ 1.6e5 frames (thousands per field accumulated
        // over the six-year baseline).
        let n = frames_to_reach_depth(published::ZTF_SINGLE_DEPTH, published::STACKED_DEPTH_REACH);
        assert!(
            (1.0e5..3.0e5).contains(&n),
            "frames to reach 27th mag = {n:.3e}"
        );
        let depth = stacked_limiting_mag(published::ZTF_SINGLE_DEPTH, n.round() as usize);
        assert_relative_eq!(depth, published::STACKED_DEPTH_REACH, epsilon = 0.05);
        // a few thousand frames already buys ~4 magnitudes.
        assert!(stack_depth_gain_mag(3000) > 4.0);
    }
}
