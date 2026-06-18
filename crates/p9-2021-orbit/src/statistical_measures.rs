//! Clustering significance statistics from Brown & Batygin (2021).
//!
//! Circular statistics come from `p9_core::analysis::circular` (the
//! small-n-corrected Rayleigh machinery); the observed longitudes of
//! perihelion come from the vetted ETNO table in `p9_core::data::etno`
//! (Brown 2017 Table 1 cross-checked against JPL SBDB) instead of an
//! unsourced hard-coded list.
//!
//! In addition to the textbook Rayleigh test, a seeded Monte Carlo
//! significance is provided with either a uniform null or a survey-bias
//! weighted null (the paper's significance is bias-aware: discovery surveys
//! avoid the galactic plane, which the ecliptic crosses near longitudes of
//! ~95 and ~275 degrees, so even unclustered ETNOs would not have uniform
//! observed perihelion longitudes).

use rand::Rng;
use serde::{Deserialize, Serialize};

use p9_core::analysis::circular::{circular_mean, mean_resultant_length, rayleigh_p_value};
use p9_core::constants::{DEG2RAD, TWO_PI};
use p9_core::data::etno::longitudes_of_perihelion;
use p9_core::units::{radians, Angle};

/// Result of the clustering confidence analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClusteringConfidence {
    /// Rayleigh test statistic Z = n * R^2
    pub rayleigh_z: f64,
    /// Mean resultant length R (0 = uniform, 1 = perfectly aligned)
    pub mean_resultant_length: f64,
    /// Mean direction of clustering (radians)
    pub mean_direction: f64,
    /// p-value from the Rayleigh test
    pub p_value: f64,
    /// Confidence level (1 - p_value)
    pub confidence: f64,
    /// Number of objects analyzed
    pub n_objects: usize,
}

impl ClusteringConfidence {
    /// Mean direction of clustering as a dimension-checked [`Angle`].
    pub fn mean_direction_angle(&self) -> Angle {
        radians(self.mean_direction)
    }
}

/// Compute clustering significance using the (analytic, small-n-corrected)
/// Rayleigh test on varpi values.
pub fn compute_clustering_significance(varpi_values: &[f64]) -> ClusteringConfidence {
    let n = varpi_values.len();
    assert!(n > 0, "Need at least one varpi value");

    let r = mean_resultant_length(varpi_values);
    let z = n as f64 * r * r;
    let mean_direction = circular_mean(varpi_values)
        .map(|m| m.rem_euclid(TWO_PI))
        .unwrap_or(0.0);
    let p_value = rayleigh_p_value(varpi_values);

    ClusteringConfidence {
        rayleigh_z: z,
        mean_resultant_length: r,
        mean_direction,
        p_value,
        confidence: 1.0 - p_value,
        n_objects: n,
    }
}

/// Null hypothesis model for the Monte Carlo significance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NullModel {
    /// Isotropic longitudes of perihelion (the null the paper's critics
    /// used; kept as a labeled baseline).
    Uniform,
    /// Survey-bias weighted: perihelion longitudes suppressed near the
    /// galactic-plane crossings of the ecliptic (λ ≈ 95° and 275°), since
    /// high-e ETNOs are discovered near perihelion and surveys avoid the
    /// plane. Simplified single-factor model; see `survey_bias_weight`.
    SurveyBias,
}

/// Relative discovery weight for a perihelion longitude (radians) under the
/// simplified survey-bias model: unit weight away from the galactic plane,
/// suppressed by `1 − depth·exp(−Δλ²/2σ²)` near each crossing
/// (depth = 0.85, σ = 20°).
pub fn survey_bias_weight(varpi: f64) -> f64 {
    let depth = 0.85;
    let sigma = 20.0 * DEG2RAD;
    let mut w = 1.0;
    for crossing_deg in [95.0, 275.0] {
        let mut d = (varpi - crossing_deg * DEG2RAD).rem_euclid(TWO_PI);
        if d > std::f64::consts::PI {
            d -= TWO_PI;
        }
        w -= depth * (-d * d / (2.0 * sigma * sigma)).exp();
    }
    w.max(0.0)
}

/// Seeded Monte Carlo clustering significance.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct McSignificance {
    pub observed_r: f64,
    pub p_value: f64,
    pub confidence: f64,
    pub n_iterations: usize,
}

/// Monte Carlo p-value for the observed mean resultant length under the
/// chosen null: the fraction of synthetic samples (same n) whose R̄ meets or
/// exceeds the observed one.
pub fn clustering_significance_mc(
    varpis: &[f64],
    n_iterations: usize,
    seed: u64,
    null: NullModel,
) -> McSignificance {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let n = varpis.len();
    let observed_r = mean_resultant_length(varpis);

    let draw_angle = |rng: &mut rand::rngs::StdRng| -> f64 {
        match null {
            NullModel::Uniform => rng.gen_range(0.0..TWO_PI),
            NullModel::SurveyBias => loop {
                let lambda = rng.gen_range(0.0..TWO_PI);
                if rng.gen::<f64>() < survey_bias_weight(lambda) {
                    break lambda;
                }
            },
        }
    };

    let mut n_exceed = 0usize;
    let mut sample = vec![0.0; n];
    for _ in 0..n_iterations {
        for s in sample.iter_mut() {
            *s = draw_angle(&mut rng);
        }
        if mean_resultant_length(&sample) >= observed_r {
            n_exceed += 1;
        }
    }

    let p_value = n_exceed as f64 / n_iterations as f64;
    McSignificance {
        observed_r,
        p_value,
        confidence: 1.0 - p_value,
        n_iterations,
    }
}

/// Clustering confidence for the vetted ETNO sample (Brown 2017 Table 1 via
/// `p9_core::data::etno`), analytic Rayleigh version.
///
/// Brown & Batygin (2021) report 99.6% from their 11-object sample with the
/// full bias treatment; this vetted 10-object table with the analytic
/// Rayleigh test lands within a point of that.
pub fn paper_clustering_confidence() -> ClusteringConfidence {
    let varpis = longitudes_of_perihelion();
    compute_clustering_significance(&varpis)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use uom::si::angle::radian;

    #[test]
    fn typed_mean_direction_matches_f64() {
        let result = compute_clustering_significance(&[0.7, 0.8, 0.9, 1.0]);
        assert_relative_eq!(
            result.mean_direction_angle().get::<radian>(),
            result.mean_direction,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_perfectly_clustered() {
        let angles = vec![1.0, 1.0, 1.0, 1.0, 1.0];
        let result = compute_clustering_significance(&angles);
        assert_relative_eq!(result.mean_resultant_length, 1.0, epsilon = 1e-10);
        assert!(result.p_value < 0.01);
        assert!(result.confidence > 0.99);
    }

    /// M8 paper-number regression: confidence pinned against the published
    /// 99.6%, within the tolerance expected from the sample difference
    /// (vetted 10-object Brown 2017 table vs the paper's 11 objects) and
    /// the analytic-vs-MC approximation. Previously this test asserted only
    /// > 0.95.
    #[test]
    fn test_paper_clustering_confidence_pinned() {
        let result = paper_clustering_confidence();
        assert_eq!(result.n_objects, 10);
        assert!(
            (result.confidence - 0.996).abs() < 0.01,
            "confidence {:.4} should be within 0.01 of the paper's 0.996",
            result.confidence
        );
        assert!(
            result.mean_resultant_length > 0.5,
            "R = {:.3} should show clustering",
            result.mean_resultant_length
        );
    }

    /// The seeded uniform-null MC agrees with the analytic Rayleigh p-value.
    #[test]
    fn test_mc_uniform_null_matches_analytic() {
        let varpis = longitudes_of_perihelion();
        let analytic = compute_clustering_significance(&varpis);
        let mc = clustering_significance_mc(&varpis, 50_000, 42, NullModel::Uniform);
        assert!(
            (mc.p_value - analytic.p_value).abs() < 0.005,
            "MC p = {:.4} vs analytic p = {:.4}",
            mc.p_value,
            analytic.p_value
        );
    }

    /// The bias-aware null is more permissive than uniform: apparent
    /// clustering is partly expected from the survey footprint, so the
    /// p-value cannot decrease.
    #[test]
    fn test_bias_aware_null_is_more_conservative() {
        let varpis = longitudes_of_perihelion();
        let uniform = clustering_significance_mc(&varpis, 20_000, 7, NullModel::Uniform);
        let biased = clustering_significance_mc(&varpis, 20_000, 7, NullModel::SurveyBias);
        assert!(
            biased.p_value >= uniform.p_value,
            "bias-aware p {:.4} < uniform p {:.4}",
            biased.p_value,
            uniform.p_value
        );
        // But the observed clustering should still be significant.
        assert!(
            biased.confidence > 0.9,
            "confidence {:.3}",
            biased.confidence
        );
    }

    #[test]
    fn test_survey_bias_weight_shape() {
        // Suppressed at the crossings, near unity at the cluster direction.
        assert!(survey_bias_weight(95.0 * DEG2RAD) < 0.3);
        assert!(survey_bias_weight(275.0 * DEG2RAD) < 0.3);
        assert!(survey_bias_weight(20.0 * DEG2RAD) > 0.8);
    }
}
