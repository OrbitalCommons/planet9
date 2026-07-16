//! Seeded synthetic ϖ samples mimicking the stable-ETNO populations of
//! Siraj, Chyba & Tremaine (2026).
//!
//! Two samples are constructed with a fixed [`rand::rngs::StdRng`] seed so the
//! computed significances are reproducible:
//!
//! - [`stable_21`] — n = 21 objects clustered around a mean ϖ with a
//!   concentration tuned so the estimator yields ≈2.7σ.
//! - [`stable_25`] — the same 21 objects plus 4 additional, less-aligned
//!   stable objects (broader scatter, offset mean), yielding ≈1.9σ.
//!
//! These synthetic samples are TUNED (the concentrations were solved from
//! the published 2.7σ/1.9σ outputs), so the headline values they reproduce
//! are circular — they demonstrate only the dilution MECHANISM (same
//! estimator, larger diluted sample ⇒ lower significance). The genuine data
//! test is [`vetted_etno_varpi`]: the real 10-object ϖ sample run through
//! the same estimator (pinned in the tests, and documented in
//! REPRODUCTION_NOTES).
//!
//! Wrapped-normal draws are used (σ ↔ κ via R̄ = exp(−σ²/2)); the von Mises
//! fit then recovers an effective κ from the realized sample.

/// The REAL vetted ETNO ϖ sample from `p9_core::data::etno` (10 objects,
/// Brown 2017 selection) — the workspace's actual data, previously never fed
/// to the estimator (only tuned synthetics were). Not the paper's exact
/// stable-21/25 tables (those element lists are not transcribed here), but a
/// genuine, independent sample the estimator must find clustered.
pub fn vetted_etno_varpi() -> Vec<f64> {
    p9_core::data::etno::longitudes_of_perihelion()
}

use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};
use std::f64::consts::TAU;

/// Mean longitude of perihelion ϖ* of the clustered stable core (rad).
/// Loosely matches the ETNO clustering direction near ϖ ≈ 50°.
const VARPI_STAR: f64 = 0.873; // ~50 deg

/// Seed for the reproducible sample construction.
const SEED: u64 = 2026;

/// Angular scatter of the clustered core (rad). Tuned so the
/// log-likelihood-ratio estimator returns ≈2.7σ on the n = 21 sample.
const CORE_SCATTER: f64 = 1.20;

/// Offset of the 4 less-aligned objects from `VARPI_STAR` (rad, ~103°) and
/// their (broader) scatter. Tuned so the n = 25 significance drops to ≈1.9σ.
const EXTRA_OFFSET: f64 = 1.80;
const EXTRA_SCATTER: f64 = 1.00;

/// The n = 21 stable-ETNO sample: clustered around `VARPI_STAR`.
///
/// The angular scatter (`CORE_SCATTER` rad) is tuned so the
/// log-likelihood-ratio estimator returns a significance ≈2.7σ.
pub fn stable_21() -> Vec<f64> {
    let mut rng = StdRng::seed_from_u64(SEED);
    core_21(&mut rng)
}

/// The n = 25 sample: the 21 clustered objects plus 4 less-aligned stable
/// objects drawn from a broader, offset distribution. Recomputing the same
/// estimator on this diluted sample drops the significance to ≈1.9σ.
pub fn stable_25() -> Vec<f64> {
    // Identical seed/stream so the first 21 draws match `stable_21` exactly;
    // the 4 extras then come from a wider, offset wrapped-normal.
    let mut rng = StdRng::seed_from_u64(SEED);
    let mut sample = core_21(&mut rng);

    // Less-aligned additions: broader scatter, mean pulled ~103° away so they
    // dilute rather than reinforce the clustering.
    let wide = Normal::new(VARPI_STAR + EXTRA_OFFSET, EXTRA_SCATTER).unwrap();
    for _ in 0..4 {
        let x: f64 = wide.sample(&mut rng);
        sample.push(x.rem_euclid(TAU));
    }
    sample
}

/// Draw the 21-object clustered core from the supplied RNG stream.
fn core_21(rng: &mut StdRng) -> Vec<f64> {
    let tight = Normal::new(VARPI_STAR, CORE_SCATTER).unwrap();
    (0..21)
        .map(|_| {
            let x: f64 = tight.sample(rng);
            x.rem_euclid(TAU)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_sizes() {
        assert_eq!(stable_21().len(), 21);
        assert_eq!(stable_25().len(), 25);
    }

    #[test]
    fn stable_25_extends_stable_21() {
        // The first 21 entries of the n = 25 sample are exactly the n = 21
        // sample (same seed and RNG stream).
        let s21 = stable_21();
        let s25 = stable_25();
        for (a, b) in s21.iter().zip(s25.iter()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn real_etno_sample_is_significantly_clustered() {
        // The genuine data test issue #248 asked for: the vetted 10-object
        // ϖ sample (p9-core, Brown 2017 selection) through the SAME
        // estimator gives κ ≈ 1.73 and σ ≈ 2.64 — independently landing at
        // the paper's ~2.7σ scale with real data, unlike the tuned
        // synthetics (which reproduce their targets by construction).
        let v = vetted_etno_varpi();
        assert_eq!(v.len(), 10);
        let fit = crate::estimator::fit(&v);
        let lam = crate::estimator::log_likelihood_ratio(&v, fit.mu, fit.kappa);
        let sig = crate::significance::sigma(lam);
        assert!((1.6..1.9).contains(&fit.kappa), "kappa = {:.3}", fit.kappa);
        assert!((2.3..3.0).contains(&sig), "sigma = {sig:.3}");

        // Small-sample cross-check: the asymptotic Wilks p = exp(−Λ) vs the
        // small-n-corrected Rayleigh series from p9-core on the same sample.
        // At n = 10 they agree to within ~0.3σ; the Wilks form is the more
        // optimistic (uncorrected O(1/n) bias), documented here.
        let p_wilks = crate::significance::lambda_to_p_value(lam);
        let p_rayleigh = p9_core::analysis::circular::rayleigh_p_value(&v);
        let s_wilks = crate::significance::p_value_to_sigma(p_wilks);
        let s_rayleigh = crate::significance::p_value_to_sigma(p_rayleigh);
        assert!(
            (s_wilks - s_rayleigh).abs() < 0.5,
            "Wilks {s_wilks:.2}σ vs small-n Rayleigh {s_rayleigh:.2}σ"
        );
    }
}
