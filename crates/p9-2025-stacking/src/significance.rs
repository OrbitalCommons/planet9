//! The look-elsewhere (trials) penalty for a blind stacking search.
//!
//! Searching `n_trials` independent orbit hypotheses inflates the chance of a
//! noise fluctuation passing threshold somewhere: the global false-alarm
//! probability is ≈ n_trials × (per-trial tail). To hold a fixed *global*
//! false-alarm probability the per-trial detection threshold must rise to
//! z ≈ Φ⁻¹(1 − α/n_trials) ≈ √(2 ln n_trials). The point of the paper — and the
//! reason searching billions of orbits is worthwhile — is that this penalty is
//! only *logarithmic* in the trial count: it costs a few tenths of a magnitude,
//! against the several magnitudes the stack itself buys.

use crate::matched_filter::stack_depth_gain_mag;

/// Standard-normal quantile Φ⁻¹(p) (Acklam's rational approximation,
/// |error| ≲ 1e-9 across the tails).
pub fn normal_quantile(p: f64) -> f64 {
    assert!(p > 0.0 && p < 1.0, "p must be in (0,1), got {p}");
    const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    let plow = 0.02425;
    let phigh = 1.0 - plow;
    if p < plow {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= phigh {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

/// Per-trial detection threshold (in σ) needed to keep the *global* false-alarm
/// probability at `global_alpha` while searching `n_trials` hypotheses:
/// z = Φ⁻¹(1 − α/n_trials), computed via the lower tail for accuracy.
pub fn look_elsewhere_threshold_sigma(n_trials: f64, global_alpha: f64) -> f64 {
    let per_trial_tail = global_alpha / n_trials;
    -normal_quantile(per_trial_tail)
}

/// The √(2 ln N) asymptotic form of the trials threshold (no α dependence).
pub fn asymptotic_threshold_sigma(n_trials: f64) -> f64 {
    (2.0 * n_trials.ln()).sqrt()
}

/// Extra limiting-magnitude cost of raising the detection threshold from
/// `baseline_sigma` (e.g. a naive 5σ) to `threshold_sigma`: Δm = 2.5·log₁₀(z/z₀)
/// (limiting flux ∝ threshold SNR).
pub fn trials_penalty_mag(threshold_sigma: f64, baseline_sigma: f64) -> f64 {
    2.5 * (threshold_sigma / baseline_sigma).log10()
}

/// Net limiting-magnitude gain of the search: the √N stacking depth gain minus
/// the look-elsewhere penalty of searching `n_trials` orbits at global FAP
/// `global_alpha` (relative to a `baseline_sigma` single-trial threshold).
pub fn net_depth_gain_mag(
    n_frames: usize,
    n_trials: f64,
    global_alpha: f64,
    baseline_sigma: f64,
) -> f64 {
    let z = look_elsewhere_threshold_sigma(n_trials, global_alpha);
    stack_depth_gain_mag(n_frames) - trials_penalty_mag(z, baseline_sigma)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn normal_quantile_known_points() {
        assert_relative_eq!(normal_quantile(0.975), 1.959964, epsilon = 1e-4);
        assert_relative_eq!(normal_quantile(0.99865), 3.0, epsilon = 2e-3);
        assert_relative_eq!(normal_quantile(0.5), 0.0, epsilon = 1e-9);
    }

    #[test]
    fn billions_of_trials_need_about_six_sigma() {
        // ~10^9 trial orbits at 1% global false-alarm → ~6.5σ per-trial.
        let z = look_elsewhere_threshold_sigma(1.0e9, 0.01);
        assert!((6.0..7.2).contains(&z), "z = {z}");
        // the √(2 ln N) rule is close.
        assert_relative_eq!(asymptotic_threshold_sigma(1.0e9), 6.44, epsilon = 0.1);
    }

    #[test]
    fn trials_penalty_is_cheap_versus_stacking_gain() {
        // The billions-of-trials penalty costs only a few tenths of a mag …
        let z = look_elsewhere_threshold_sigma(1.0e9, 0.01);
        let penalty = trials_penalty_mag(z, 5.0);
        assert!(penalty < 0.5, "penalty = {penalty} mag");
        // … while stacking ~10^5 frames buys ~6 mag, so the net gain stays
        // multi-magnitude despite searching a billion orbits.
        let net = net_depth_gain_mag(100_000, 1.0e9, 0.01, 5.0);
        assert!(net > 5.0, "net gain = {net} mag");
    }

    #[test]
    fn threshold_grows_with_trials_but_only_logarithmically() {
        let z6 = look_elsewhere_threshold_sigma(1.0e6, 0.01);
        let z12 = look_elsewhere_threshold_sigma(1.0e12, 0.01);
        // a million-fold more trials raises the bar by < 2.5σ.
        assert!(z12 > z6 && z12 - z6 < 2.5, "z6={z6}, z12={z12}");
    }
}
