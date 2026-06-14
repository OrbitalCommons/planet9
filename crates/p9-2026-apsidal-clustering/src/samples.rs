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
//! The samples are NOT the answer: the significances are computed from them by
//! the estimator. The point of the construction is that the SAME estimator on
//! the larger, diluted sample returns a LOWER significance.
//!
//! Wrapped-normal draws are used (σ ↔ κ via R̄ = exp(−σ²/2)); the von Mises
//! fit then recovers an effective κ from the realized sample.

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
}
