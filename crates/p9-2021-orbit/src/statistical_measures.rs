//! Clustering significance statistics from Brown & Batygin (2021).
//!
//! The paper runs 121 N-body simulations with different Planet Nine
//! parameters drawn from the posterior, measuring orbital clustering
//! of distant KBOs. The Rayleigh test is used to assess the significance
//! of clustering in longitude of perihelion (varpi).
//!
//! Key result: 99.6% confidence that distant KBO orbits are clustered.

use serde::{Deserialize, Serialize};

use p9_core::constants::TWO_PI;

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

/// Compute clustering significance using the Rayleigh test on varpi values.
///
/// The Rayleigh test assesses whether a set of angles is uniformly distributed
/// on the circle. The test statistic is Z = n * R^2, where R is the mean
/// resultant length of the unit vectors at the given angles.
///
/// For the paper's 11 ETNOs with a > 250 AU, the observed clustering in
/// longitude of perihelion gives a confidence of 99.6%.
///
/// # Arguments
/// * `varpi_values` - Longitude of perihelion values in radians
pub fn compute_clustering_significance(varpi_values: &[f64]) -> ClusteringConfidence {
    let n = varpi_values.len();
    assert!(n > 0, "Need at least one varpi value");

    let cos_sum: f64 = varpi_values.iter().map(|v| v.cos()).sum();
    let sin_sum: f64 = varpi_values.iter().map(|v| v.sin()).sum();

    let c_bar = cos_sum / n as f64;
    let s_bar = sin_sum / n as f64;

    let r = (c_bar * c_bar + s_bar * s_bar).sqrt();
    let z = n as f64 * r * r;

    let mean_direction = s_bar.atan2(c_bar);
    let mean_direction = if mean_direction < 0.0 {
        mean_direction + TWO_PI
    } else {
        mean_direction
    };

    let p_value = rayleigh_p_value(z, n);

    ClusteringConfidence {
        rayleigh_z: z,
        mean_resultant_length: r,
        mean_direction,
        p_value,
        confidence: 1.0 - p_value,
        n_objects: n,
    }
}

/// Compute the p-value for the Rayleigh test.
///
/// For large n, P(Z > z) ~ exp(-z) with corrections for small samples.
/// Uses the exact formula from Mardia & Jupp (2000):
///   p = exp(-z) * (1 + (2z - z^2) / (4n) - (24z - 132z^2 + 76z^3 - 9z^4) / (288n^2))
fn rayleigh_p_value(z: f64, n: usize) -> f64 {
    let nf = n as f64;

    if z > 700.0 {
        return 0.0;
    }

    let base = (-z).exp();
    let correction1 = (2.0 * z - z * z) / (4.0 * nf);
    let correction2 =
        (24.0 * z - 132.0 * z * z + 76.0 * z * z * z - 9.0 * z.powi(4)) / (288.0 * nf * nf);

    let p = base * (1.0 + correction1 - correction2);
    p.clamp(0.0, 1.0)
}

/// Compute the expected confidence level for the Brown & Batygin (2021) sample.
///
/// The paper reports 99.6% confidence from 11 ETNOs.
/// This returns the reference clustering confidence using the
/// approximate observed varpi values from Table 1 of the paper.
pub fn paper_clustering_confidence() -> ClusteringConfidence {
    let varpi_deg: Vec<f64> = vec![
        68.0, 96.0, 12.0, 318.0, 7.0, 66.0, 61.0, 78.0, 24.0, 354.0, 35.0,
    ];
    let varpi_rad: Vec<f64> = varpi_deg
        .iter()
        .map(|d| d * p9_core::constants::DEG2RAD)
        .collect();

    compute_clustering_significance(&varpi_rad)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_perfectly_clustered() {
        let angles = vec![1.0, 1.0, 1.0, 1.0, 1.0];
        let result = compute_clustering_significance(&angles);
        assert_relative_eq!(result.mean_resultant_length, 1.0, epsilon = 1e-10);
        assert!(result.p_value < 0.01);
        assert!(result.confidence > 0.99);
    }

    #[test]
    fn test_paper_clustering_confidence() {
        let result = paper_clustering_confidence();
        assert_eq!(result.n_objects, 11);
        assert!(
            result.confidence > 0.95,
            "Confidence {:.3} should exceed 95%",
            result.confidence
        );
        assert!(
            result.mean_resultant_length > 0.3,
            "R = {:.3} should show clustering",
            result.mean_resultant_length
        );
    }

    #[test]
    fn test_rayleigh_p_value_properties() {
        let p_small = rayleigh_p_value(0.1, 10);
        let p_large = rayleigh_p_value(5.0, 10);
        assert!(
            p_small > p_large,
            "Larger Z should give smaller p: {:.4} vs {:.4}",
            p_small,
            p_large
        );

        let p_zero = rayleigh_p_value(0.0, 10);
        assert_relative_eq!(p_zero, 1.0, epsilon = 0.05);
    }
}
