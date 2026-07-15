//! Matched-filter stacking and the orbit-parameter search metric.
//!
//! Shared primitives lifted from Geringer-Sameth, Golovich & Iwabuchi (2025),
//! "Multi-year stacking searches for solar system bodies" (arXiv:2509.25428):
//!
//! - [`matched_filter`] — the matched-filter stack: co-adding N background-
//!   limited frames along the correct track grows the signal-to-noise as √N
//!   (limiting magnitude as 1.25·log₁₀N).
//! - [`orbit_metric`] — the Fisher metric on orbit space: a trial orbit that
//!   misses the true track loses SNR quadratically in the parameter error,
//!   setting the grid spacing and hence the number of trial orbits.

/// Matched-filter stacking: SNR and limiting-magnitude gain from co-adding
/// frames along a track.
///
/// For N background-limited exposures of equal depth, the optimal (matched-
/// filter) combination adds the signal coherently and the noise in quadrature,
/// so the stacked SNR is √N times a single frame. In flux terms the limiting
/// flux falls as 1/√N, i.e. the limiting magnitude deepens by
/// Δm = 2.5·log₁₀(√N) = 1.25·log₁₀(N).
pub mod matched_filter {
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
        use approx::assert_relative_eq;

        // Labelled reference numbers from the paper (regression comparison only).
        const ZTF_SINGLE_DEPTH: f64 = 20.5;
        const STACKED_DEPTH_REACH: f64 = 27.0;

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
            let n = frames_to_reach_depth(ZTF_SINGLE_DEPTH, STACKED_DEPTH_REACH);
            assert!(
                (1.0e5..3.0e5).contains(&n),
                "frames to reach 27th mag = {n:.3e}"
            );
            let depth = stacked_limiting_mag(ZTF_SINGLE_DEPTH, n.round() as usize);
            assert_relative_eq!(depth, STACKED_DEPTH_REACH, epsilon = 0.05);
            // a few thousand frames already buys ~4 magnitudes.
            assert!(stack_depth_gain_mag(3000) > 4.0);
        }
    }
}

/// The orbit-parameter metric: how fast SNR is lost when a trial orbit misses
/// the true track, and the resulting number of trial orbits to grid the search.
///
/// Model a linear on-sky track over a baseline `T` (days) imaged with a
/// Gaussian PSF of width `θ` (1σ, arcsec). A trial with a rate error `δv`
/// (arcsec/day) leaves the source progressively offset from the template; at
/// time `t` (relative to mid-baseline) the offset is `δv·t`. The matched-filter
/// stacked signal is the mean per-frame Gaussian overlap
/// `⟨exp(−Δ(t)²/4θ²)⟩`, so to leading order the fractional SNR loss is
///
///   L(δv) = ⟨Δ²⟩ / (4θ²) = (δv² T²/12) / (4θ²) = ½ g_vv δv²,
///   g_vv = T² / (24 θ²).
///
/// `g_vv` is one element of the Fisher metric on orbit space (Geringer-Sameth
/// et al. 2025). It sharpens as `T²`, so a multi-year baseline needs an
/// extremely fine rate grid — hence *billions* of trial orbits — while also
/// making each individual track far more constraining.
pub mod orbit_metric {
    use crate::constants::{GM_SUN, PC_AU};
    use crate::types::P9Params;

    /// Rate-error metric coefficient `g_vv = T²/(24 θ²)` (units: (arcsec/day)⁻²).
    pub fn metric_g_vv(psf_arcsec: f64, baseline_days: f64) -> f64 {
        baseline_days * baseline_days / (24.0 * psf_arcsec * psf_arcsec)
    }

    /// Quadratic (leading-order) retained SNR fraction for a rate error.
    pub fn snr_retained_quadratic(psf_arcsec: f64, baseline_days: f64, rate_error: f64) -> f64 {
        let l = 0.5 * metric_g_vv(psf_arcsec, baseline_days) * rate_error * rate_error;
        (1.0 - l).max(0.0)
    }

    /// Exact retained SNR fraction: the mean per-frame Gaussian overlap over a
    /// uniformly-sampled baseline.
    pub fn snr_retained_exact(
        psf_arcsec: f64,
        baseline_days: f64,
        rate_error: f64,
        n_samples: usize,
    ) -> f64 {
        let n = n_samples.max(2);
        let theta2 = psf_arcsec * psf_arcsec;
        let mut sum = 0.0;
        for k in 0..n {
            // t uniform in [-T/2, T/2]
            let t = baseline_days * (k as f64 / (n as f64 - 1.0) - 0.5);
            let delta = rate_error * t;
            sum += (-delta * delta / (4.0 * theta2)).exp();
        }
        sum / n as f64
    }

    /// Rate resolution (cell half-width, arcsec/day) at which the SNR loss reaches
    /// the tolerance `1 − snr_tolerance`: δv = (θ/T)·√(48·(1−η)).
    pub fn rate_resolution(psf_arcsec: f64, baseline_days: f64, snr_tolerance: f64) -> f64 {
        (psf_arcsec / baseline_days) * (48.0 * (1.0 - snr_tolerance)).sqrt()
    }

    /// Number of rate cells along one sky axis spanning `rate_range` (arcsec/day),
    /// at the given SNR tolerance and baseline.
    pub fn n_rate_cells_1d(
        rate_range: f64,
        psf_arcsec: f64,
        baseline_days: f64,
        snr_tolerance: f64,
    ) -> f64 {
        rate_range / (2.0 * rate_resolution(psf_arcsec, baseline_days, snr_tolerance))
    }

    /// Number of trial orbits for a linear two-dimensional on-sky velocity search
    /// (the rate hypotheses; start position is read off the stacked image, not a
    /// separate trial). Equals the square of the per-axis rate-cell count.
    pub fn n_trial_orbits(
        rate_range: f64,
        psf_arcsec: f64,
        baseline_days: f64,
        snr_tolerance: f64,
    ) -> f64 {
        let n1d = n_rate_cells_1d(rate_range, psf_arcsec, baseline_days, snr_tolerance);
        n1d * n1d
    }

    /// Planet Nine's **heliocentric** orbital angular rate projected on the sky
    /// (arcsec/day), from its mean motion n = √(GM☉/a³).
    ///
    /// WARNING: this is NOT the on-sky search rate. The apparent motion of a
    /// distant body seen from Earth is dominated by Earth's parallactic reflex
    /// (~9″/day at 380 AU near opposition, ~20× this heliocentric term).
    /// Sizing a shift-stack rate grid from this value will miss the target;
    /// use [`apparent_sky_rate_at_opposition`] instead.
    pub fn p9_orbital_sky_rate(p9: &P9Params) -> f64 {
        let n_rad_per_day = (GM_SUN / p9.a.powi(3)).sqrt();
        n_rad_per_day * PC_AU // rad/day × (arcsec/rad)
    }

    /// Peak apparent on-sky rate (arcsec/day) of a distant solar-system body
    /// at heliocentric distance `distance_au`, observed near opposition.
    ///
    /// Near opposition the dominant term is Earth's parallactic reflex v_E/Δ,
    /// reduced by the body's own prograde circular-orbit motion
    /// v(d) = v_E/√d (coplanar approximation):
    ///
    ///   μ ≈ 206265 · v_E (1 − 1/√d) / (d − 1)   [arcsec/day]
    ///
    /// with v_E = √(GM☉/1 AU) ≈ 0.0172 AU/day. Examples: ≈45″/day at 70 AU,
    /// ≈8.9″/day at 380 AU, ≈7.5″/day at 450 AU — the scale a shift-stack
    /// trial-rate box must cover.
    pub fn apparent_sky_rate_at_opposition(distance_au: f64) -> f64 {
        let v_earth = GM_SUN.sqrt(); // AU/day on a 1 AU circular orbit
        let v_body = v_earth / distance_au.sqrt();
        let delta = (distance_au - 1.0).max(1.0);
        PC_AU * (v_earth - v_body) / delta
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use approx::assert_relative_eq;

        // ZTF-like single field, six-year baseline.
        const PSF: f64 = 1.0; // arcsec, 1σ
        const T6YR: f64 = 6.0 * 365.25; // days
        const RATE_RANGE: f64 = 40.0; // arcsec/day, covering distant→fast movers

        #[test]
        fn apparent_rate_dominated_by_parallactic_reflex() {
            // At 380 AU (BB21 median) the apparent opposition rate is ~8.9"/day,
            // ~20x the heliocentric mean-motion rate the old API returned.
            let apparent = apparent_sky_rate_at_opposition(380.0);
            assert!((apparent - 8.9).abs() < 0.3, "apparent = {apparent:.2}");
            let helio = p9_orbital_sky_rate(&P9Params::mcmc_2021());
            assert!(
                apparent > 15.0 * helio,
                "apparent {apparent:.2} vs heliocentric {helio:.2} arcsec/day"
            );
            // Nearby TNO distances need a ~±50"/day rate box (TESS regime).
            let near = apparent_sky_rate_at_opposition(70.0);
            assert!((near - 45.3).abs() < 1.0, "70 AU rate = {near:.1}");
        }

        #[test]
        fn metric_quadratic_matches_exact_for_small_error() {
            // pick the rate error that the quadratic metric says costs 5% SNR …
            let dv = rate_resolution(PSF, T6YR, 0.95);
            assert_relative_eq!(snr_retained_quadratic(PSF, T6YR, dv), 0.95, epsilon = 1e-9);
            // … and the exact Gaussian-overlap stack agrees to ~1%.
            let exact = snr_retained_exact(PSF, T6YR, dv, 4001);
            assert!((exact - 0.95).abs() < 0.02, "exact retained = {exact}");
        }

        #[test]
        fn rate_resolution_sharpens_inversely_with_baseline() {
            let r_full = rate_resolution(PSF, T6YR, 0.9);
            let r_half = rate_resolution(PSF, T6YR / 2.0, 0.9);
            assert_relative_eq!(r_half / r_full, 2.0, epsilon = 1e-9);
        }

        #[test]
        fn ztf_six_year_search_needs_order_billion_trial_orbits() {
            let n = n_trial_orbits(RATE_RANGE, PSF, T6YR, 0.9);
            assert!((1.0e7..1.0e11).contains(&n), "trial orbits = {n:.3e}");
        }

        #[test]
        fn trial_orbits_scale_as_baseline_squared() {
            let n_full = n_trial_orbits(RATE_RANGE, PSF, T6YR, 0.9);
            let n_half = n_trial_orbits(RATE_RANGE, PSF, T6YR / 2.0, 0.9);
            assert_relative_eq!(n_full / n_half, 4.0, epsilon = 1e-6);
        }

        #[test]
        fn multi_year_dwarfs_single_night_in_trial_count() {
            let n_6yr = n_trial_orbits(RATE_RANGE, PSF, T6YR, 0.9);
            let n_night = n_trial_orbits(RATE_RANGE, PSF, 0.3, 0.9);
            assert!(n_6yr / n_night > 1.0e6, "ratio = {:.2e}", n_6yr / n_night);
        }

        #[test]
        fn planet_nine_sits_inside_and_is_resolved_by_the_grid() {
            let rate = p9_orbital_sky_rate(&P9Params::mcmc_2021());
            // P9 is a slow mover (~tenths of an arcsec/day) …
            assert!(
                (0.1..1.0).contains(&rate),
                "P9 sky rate = {rate} arcsec/day"
            );
            // … comfortably inside the search's rate range …
            assert!(rate < RATE_RANGE);
            // … and far coarser than the six-year rate resolution, so it is well
            // localised by the metric (many cells across its rate).
            assert!(rate > 10.0 * rate_resolution(PSF, T6YR, 0.9));
        }
    }
}
