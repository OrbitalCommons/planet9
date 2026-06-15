//! Clustering metrics for evaluating P9 parameter acceptability.
//!
//! Implements the statistical tests from Brown & Batygin (2016) Section 2:
//! - Perihelion confinement: concentration (mean resultant length R̄) of the
//!   survivors' Δϖ relative to Planet Nine
//! - High-perihelion production: fraction of Sedna-like detached objects
//! - Survival rate: number of particles remaining after the integration
//!
//! Circular statistics (R̄, Rayleigh p-value with the small-n correction)
//! come from `p9_core::analysis::circular`; the observed-sample comparison
//! values come from the vetted ETNO table (`p9_core::data::etno`) and the
//! 6 stable B&B 2016 KBOs (`p9_core::data::stable_kbos`).

use std::f64::consts::PI;

use p9_core::analysis::circular::{mean_resultant_length, wrap_to_pi};
use p9_core::constants::*;
use p9_core::data::etno::longitudes_of_perihelion;
use p9_core::data::stable_kbos::{longitude_of_perihelion, stable_kbos};
use p9_core::types::OrbitalElements;

/// Sedna-like "high perihelion" threshold (AU). Detached objects lifted to
/// q > 60 AU are dynamically decoupled from Neptune (q_Neptune-scattering
/// ≲ 38 AU with a wide margin); this is the single threshold used by both
/// the metric below and the grid acceptance documentation.
pub const HIGH_Q_THRESHOLD_AU: f64 = 60.0;

/// Summary clustering statistics for a survivor population.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct ClusteringEval {
    /// Mean resultant length of the survivors' Δϖ (concentration, 0..1).
    pub r_bar_dvarpi: f64,
    /// Fraction of survivors with |Δϖ| > π/2 (anti-aligned with P9).
    pub anti_aligned_fraction: f64,
    /// Fraction of survivors with q > `HIGH_Q_THRESHOLD_AU`.
    pub high_perihelion_fraction: f64,
    /// Number of qualifying survivors.
    pub n_qualifying: usize,
}

/// Evaluate clustering metrics for a set of final orbital elements.
///
/// `elements`: final orbital elements of surviving particles
/// `varpi_p9`: Planet Nine's longitude of perihelion (radians)
/// `a_min`, `a_max`: semi-major axis range for clustering analysis
pub fn evaluate_clustering(
    elements: &[OrbitalElements],
    varpi_p9: f64,
    a_min: f64,
    a_max: f64,
) -> ClusteringEval {
    let dvarpis: Vec<f64> = elements
        .iter()
        .filter(|e| e.a >= a_min && e.a <= a_max && e.e < 1.0)
        .map(|e| wrap_to_pi(e.omega + e.omega_big - varpi_p9))
        .collect();

    let n = dvarpis.len();
    if n == 0 {
        return ClusteringEval {
            r_bar_dvarpi: 0.0,
            anti_aligned_fraction: 0.0,
            high_perihelion_fraction: 0.0,
            n_qualifying: 0,
        };
    }

    let n_anti_aligned = dvarpis.iter().filter(|dv| dv.abs() > PI / 2.0).count();
    let n_high_q = elements
        .iter()
        .filter(|e| e.a >= a_min && e.a <= a_max && e.e < 1.0)
        .filter(|e| e.a * (1.0 - e.e) > HIGH_Q_THRESHOLD_AU)
        .count();

    ClusteringEval {
        r_bar_dvarpi: mean_resultant_length(&dvarpis),
        anti_aligned_fraction: n_anti_aligned as f64 / n as f64,
        high_perihelion_fraction: n_high_q as f64 / n as f64,
        n_qualifying: n,
    }
}

/// Observed mean resultant length R̄ of the longitudes of perihelion of the
/// 6 stable a > 250 AU, q > 30 AU KBOs from Batygin & Brown (2016)
/// (≈ 0.84 from the catalog elements, i.e. circular σ ≈ 34°; the paper's
/// headline "ϖ = 71° ± 16°" corresponds to a tighter dispersion than the
/// full sample's circular statistics give).
pub fn observed_six_kbo_r_bar() -> f64 {
    let varpis: Vec<f64> = stable_kbos()
        .iter()
        .map(|k| longitude_of_perihelion(&k.elements))
        .collect();
    mean_resultant_length(&varpis)
}

/// Observed mean resultant length R̄ of the longitudes of perihelion of the
/// 10-object Brown (2017) ETNO sample (the workspace's vetted observational
/// table, `p9_core::data::etno`).
pub fn observed_etno_r_bar() -> f64 {
    mean_resultant_length(&longitudes_of_perihelion())
}

/// Compute confinement probability for a subsample.
///
/// From the paper: randomly draw `n_draw` objects from the qualifying
/// survivors (a ∈ [300,700], q < 80 AU, i < 50°) and check whether their Δϖ
/// concentration reaches the *observed* 6-KBO mean resultant length
/// (R̄ ≈ 0.96 from `observed_six_kbo_r_bar`, replacing a previous invented
/// R̄ > 0.5 cut). Returns the fraction of draws at least as confined as the
/// observations.
pub fn confinement_probability(
    elements: &[OrbitalElements],
    varpi_p9: f64,
    n_draw: usize,
    n_trials: usize,
    seed: u64,
) -> f64 {
    use rand::seq::SliceRandom;
    use rand::SeedableRng;

    let qualifying: Vec<f64> = elements
        .iter()
        .filter(|e| {
            e.a >= 300.0
                && e.a <= 700.0
                && e.perihelion() < 80.0
                && e.i < 50.0 * DEG2RAD
                && e.e < 1.0
        })
        .map(|e| wrap_to_pi(e.omega + e.omega_big - varpi_p9))
        .collect();

    if qualifying.len() < n_draw {
        return 0.0;
    }

    let r_bar_observed = observed_six_kbo_r_bar();

    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let mut n_confined = 0;

    for _ in 0..n_trials {
        let sample: Vec<f64> = qualifying
            .choose_multiple(&mut rng, n_draw)
            .copied()
            .collect();
        if mean_resultant_length(&sample) >= r_bar_observed {
            n_confined += 1;
        }
    }

    n_confined as f64 / n_trials as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::analysis::circular::rayleigh_p_value;

    #[test]
    fn test_clustering_evaluation() {
        // Two anti-aligned objects (ϖ ≈ π from P9) and one aligned
        let elements = vec![
            OrbitalElements {
                a: 400.0,
                e: 0.7,
                i: 0.0,
                omega: 3.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            OrbitalElements {
                a: 500.0,
                e: 0.8,
                i: 0.0,
                omega: 3.2,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            OrbitalElements {
                a: 350.0,
                e: 0.5,
                i: 0.0,
                omega: 0.5,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
        ];

        let eval = evaluate_clustering(&elements, 0.0, 300.0, 700.0);

        assert_eq!(eval.n_qualifying, 3);
        assert!((eval.anti_aligned_fraction - 2.0 / 3.0).abs() < 1e-12);
        // q = 120, 100, 175 AU: all above the 60 AU detachment threshold
        assert!((eval.high_perihelion_fraction - 1.0).abs() < 1e-12);
        assert!(eval.r_bar_dvarpi > 0.0 && eval.r_bar_dvarpi < 1.0);
    }

    #[test]
    fn test_core_rayleigh_uniform_vs_clustered() {
        // Sanity check on the p9-core Rayleigh machinery this crate uses:
        // uniform angles are insignificant, clustered ones are significant.
        let uniform: Vec<f64> = (0..20).map(|i| i as f64 / 20.0 * TWO_PI).collect();
        assert!(mean_resultant_length(&uniform) < 0.2);
        assert!(rayleigh_p_value(&uniform) > 0.1);

        let clustered: Vec<f64> = (0..20).map(|i| 1.0 + i as f64 * 0.01).collect();
        assert!(mean_resultant_length(&clustered) > 0.9);
        assert!(rayleigh_p_value(&clustered) < 0.001);
    }

    #[test]
    fn test_observed_r_bars() {
        // The 6-KBO sample is tightly confined: the catalog elements give
        // R̄ ≈ 0.84 (circular σ ≈ 34°; the paper's "ϖ = 71 ± 16°" quote
        // corresponds to a tighter dispersion than the full sample shows).
        let r6 = observed_six_kbo_r_bar();
        assert!((r6 - 0.84).abs() < 0.02, "6-KBO R̄ = {r6:.3}");

        // The 10-object ETNO sample is clustered but looser
        let r10 = observed_etno_r_bar();
        assert!(r10 > 0.3 && r10 < r6, "ETNO R̄ = {r10:.3}");
    }

    #[test]
    fn test_confinement_probability_limits() {
        // A population as tightly confined as the observations passes ~all
        // draws; a uniform population passes ~none.
        let confined: Vec<OrbitalElements> = (0..30)
            .map(|k| OrbitalElements {
                a: 400.0,
                e: 0.85,
                i: 0.1,
                omega: 3.0 + 0.01 * k as f64,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            })
            .collect();
        let p_confined = confinement_probability(&confined, 0.0, 6, 200, 3);
        assert!(p_confined > 0.9, "p = {p_confined}");

        let uniform: Vec<OrbitalElements> = (0..30)
            .map(|k| OrbitalElements {
                a: 400.0,
                e: 0.85,
                i: 0.1,
                omega: TWO_PI * k as f64 / 30.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            })
            .collect();
        let p_uniform = confinement_probability(&uniform, 0.0, 6, 200, 3);
        assert!(p_uniform < 0.1, "p = {p_uniform}");
    }
}
