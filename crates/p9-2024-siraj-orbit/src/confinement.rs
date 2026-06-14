//! Mapping the observed ETNO apsidal confinement to a required secular
//! forcing strength S_obs.
//!
//! The clustered ETNOs have longitudes of perihelion ϖ that are concentrated
//! rather than uniform on the circle. We summarize that concentration with the
//! maximum-likelihood von Mises concentration parameter κ̂, estimated from the
//! sample mean resultant length R̄ (computed with the shared p9-core circular
//! statistics). A stronger secular forcing produces a tighter apsidal cluster,
//! so κ̂ maps monotonically to a required forcing strength.
//!
//! We adopt the standard secular-confinement relation in which the equilibrium
//! concentration grows linearly with the forcing strength,
//!
//!   κ = k_cal · S_obs,
//!
//! (a forced apsidal libration of small-amplitude has confinement ∝ the
//! restoring torque, i.e. ∝ S). Inverting gives S_obs = κ̂ / k_cal. The
//! single calibration constant `k_cal` sets the absolute scale; it is fixed
//! once (in `posterior`) so that the maximum-a-posteriori perturber lands on
//! the Siraj et al. (2024) reference orbit, and only enters through the ratio
//! S(m,a_p)/S_obs in the likelihood.

use p9_core::analysis::circular::mean_resultant_length;

/// Maximum-likelihood von Mises concentration κ̂ from a mean resultant length
/// R̄, using the Mardia & Jupp (2000, §5.3.1) rational approximation to the
/// inverse of A(κ) = I₁(κ)/I₀(κ):
///
///   κ̂ = 2R̄ + R̄³ + 5R̄⁵/6                       (R̄ < 0.53)
///   κ̂ = −0.4 + 1.39R̄ + 0.43/(1 − R̄)             (0.53 ≤ R̄ < 0.85)
///   κ̂ = 1/(R̄³ − 4R̄² + 3R̄)                       (R̄ ≥ 0.85)
pub fn kappa_from_resultant(r_bar: f64) -> f64 {
    if r_bar < 0.53 {
        2.0 * r_bar + r_bar.powi(3) + 5.0 * r_bar.powi(5) / 6.0
    } else if r_bar < 0.85 {
        -0.4 + 1.39 * r_bar + 0.43 / (1.0 - r_bar)
    } else {
        1.0 / (r_bar.powi(3) - 4.0 * r_bar * r_bar + 3.0 * r_bar)
    }
}

/// Estimate the von Mises concentration κ̂ directly from a sample of angles
/// (radians), via the shared circular mean resultant length.
pub fn kappa_from_angles(angles: &[f64]) -> f64 {
    kappa_from_resultant(mean_resultant_length(angles))
}

/// Approximate fractional (1σ) uncertainty on κ̂ for a von Mises sample of
/// size `n`. For a concentrated sample the large-κ asymptotic variance of the
/// MLE is Var(κ̂) ≈ 2κ²/n (Fisher 1993, §4.5.5), giving σ_κ/κ ≈ √(2/n). This
/// is the dominant uncertainty band that propagates into S_obs.
pub fn kappa_fractional_uncertainty(n: usize) -> f64 {
    (2.0 / n as f64).sqrt()
}

/// The observed required forcing strength and its (1σ) fractional uncertainty,
/// derived from the ETNO ϖ sample and a calibration constant `k_cal`:
///
///   S_obs = κ̂ / k_cal,    σ_lnS = σ_κ/κ.
///
/// Returns `(s_obs, frac_sigma)`.
pub fn observed_forcing(angles: &[f64], k_cal: f64) -> (f64, f64) {
    let kappa = kappa_from_angles(angles);
    let s_obs = kappa / k_cal;
    let frac = kappa_fractional_uncertainty(angles.len());
    (s_obs, frac)
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::TWO_PI;

    #[test]
    fn test_kappa_monotone_in_resultant() {
        // Tighter clustering (larger R̄) ⇒ larger κ̂, across all three branches.
        let ks: Vec<f64> = [0.2, 0.4, 0.6, 0.8, 0.9, 0.95]
            .iter()
            .map(|&r| kappa_from_resultant(r))
            .collect();
        for w in ks.windows(2) {
            assert!(w[1] > w[0], "kappa not monotone: {ks:?}");
        }
    }

    #[test]
    fn test_kappa_uniform_is_small() {
        // A near-uniform circle has R̄ ≈ 0 ⇒ κ̂ ≈ 0.
        let uniform: Vec<f64> = (0..20).map(|k| k as f64 * TWO_PI / 20.0).collect();
        assert!(kappa_from_angles(&uniform) < 0.2);
    }

    #[test]
    fn test_fractional_uncertainty_shrinks_with_n() {
        assert!(kappa_fractional_uncertainty(40) < kappa_fractional_uncertainty(10));
    }
}
