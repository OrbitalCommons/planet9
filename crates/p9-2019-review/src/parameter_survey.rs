//! Parameter grid survey with statistical success measures.
//!
//! Implements the statistical measures from the review paper (Figure 16):
//! - a_c: critical semi-major axis (transition from random to confined)
//! - f_varpi: apsidal confinement fraction
//! - eta: forced equilibrium magnitude
//! - mu: forced equilibrium angle
//! - sigma: RMS dispersion in (p, g) space

use p9_core::constants::*;
use p9_core::types::P9Params;

/// Statistical success criteria from the review paper.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SuccessCriteria {
    /// Critical semi-major axis should be in (200, 300) AU
    pub a_c_range: (f64, f64),
    /// Apsidal confinement fraction should be >= 0.8
    pub f_varpi_min: f64,
    /// Forced equilibrium magnitude should be ~0.1
    pub eta_target: f64,
    /// Forced equilibrium angle should be <= 10 degrees
    pub mu_max_deg: f64,
    /// RMS dispersion should be ~0.2
    pub sigma_target: f64,
}

impl SuccessCriteria {
    pub fn paper_defaults() -> Self {
        Self {
            a_c_range: (200.0, 300.0),
            f_varpi_min: 0.80,
            eta_target: 0.1,
            mu_max_deg: 10.0,
            sigma_target: 0.2,
        }
    }
}

/// Result of evaluating a single simulation against success criteria.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SimulationScore {
    pub params: P9ParamsSummary,
    pub a_c: f64,
    pub f_varpi: f64,
    pub eta: f64,
    pub mu_deg: f64,
    pub sigma: f64,
    pub passes: bool,
}

/// Summary of P9 parameters for serialization.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct P9ParamsSummary {
    pub mass_earth: f64,
    pub a: f64,
    pub e: f64,
    pub i_deg: f64,
}

impl From<&P9Params> for P9ParamsSummary {
    fn from(p: &P9Params) -> Self {
        Self {
            mass_earth: p.mass_earth,
            a: p.a,
            e: p.e,
            i_deg: p.i * RAD2DEG,
        }
    }
}

/// Definition of the parameter grid from the review paper.
///
/// Grid: a₉ ∈ (300, 1500) step 100, e₉ ∈ (0.05, 0.95) step 0.1,
///       i₉ ∈ (10, 35) step 5, m₉ ∈ {5, 10, 20} M_Earth.
pub struct ParameterGrid {
    pub a_values: Vec<f64>,
    pub e_values: Vec<f64>,
    pub i_values: Vec<f64>,
    pub mass_values: Vec<f64>,
}

impl ParameterGrid {
    /// Full grid from the paper (1,134 viable simulations).
    pub fn paper_grid() -> Self {
        let a_values: Vec<f64> = (3..=15).map(|i| i as f64 * 100.0).collect();
        let e_values: Vec<f64> = (0..=9).map(|i| 0.05 + i as f64 * 0.1).collect();
        let i_values: Vec<f64> = vec![10.0, 15.0, 20.0, 25.0, 30.0, 35.0];
        let mass_values = vec![5.0, 10.0, 20.0];

        Self {
            a_values,
            e_values,
            i_values,
            mass_values,
        }
    }

    /// Quick grid for testing.
    pub fn quick_grid() -> Self {
        Self {
            a_values: vec![400.0, 600.0, 800.0],
            e_values: vec![0.15, 0.35, 0.55],
            i_values: vec![15.0, 20.0],
            mass_values: vec![5.0, 10.0],
        }
    }

    /// Generate all P9 parameter sets, filtering for viability.
    pub fn generate(&self) -> Vec<P9Params> {
        let mut params = Vec::new();

        for &m in &self.mass_values {
            for &a in &self.a_values {
                for &e in &self.e_values {
                    for &i in &self.i_values {
                        let q = a * (1.0 - e);
                        // Viability: q > 100 AU and q < 500 AU
                        if q > 100.0 && q < 500.0 {
                            params.push(P9Params {
                                mass_earth: m,
                                a,
                                e,
                                i: i * DEG2RAD,
                                omega: 140.0 * DEG2RAD,
                                omega_big: 100.0 * DEG2RAD,
                                mean_anomaly: 0.0,
                            });
                        }
                    }
                }
            }
        }

        params
    }

    /// Count the total number of viable parameter sets.
    pub fn count_viable(&self) -> usize {
        self.generate().len()
    }
}

/// Evaluate a parameter set against success criteria.
///
/// TODO: This requires running the actual simulation. Currently returns
/// analytical estimates based on secular theory scaling relations.
pub fn evaluate_params(params: &P9Params, criteria: &SuccessCriteria) -> SimulationScore {
    // Analytical estimate of critical semi-major axis
    // a_c scales approximately as a₉ * (m₉/M_sun)^(2/7) / (1+e₉)
    let m_ratio = params.mass_earth * EARTH_MASS_SOLAR;
    let a_c = params.a * m_ratio.powf(2.0 / 7.0) / (1.0 + params.e);

    // Apsidal confinement fraction scales with mass and proximity
    let q = params.a * (1.0 - params.e);
    let f_varpi = (params.mass_earth / 10.0 * (400.0 / q).sqrt()).min(1.0);

    // Forced equilibrium magnitude scales with i₉
    let eta = 0.3 * params.i.sin();

    // Forced equilibrium angle
    let mu_deg = (10.0 / params.mass_earth * (params.a / 500.0)).min(30.0);

    // RMS dispersion
    let sigma = 0.15 + 0.1 * params.e;

    let passes = a_c >= criteria.a_c_range.0
        && a_c <= criteria.a_c_range.1
        && f_varpi >= criteria.f_varpi_min
        && mu_deg <= criteria.mu_max_deg;

    SimulationScore {
        params: P9ParamsSummary::from(params),
        a_c,
        f_varpi,
        eta,
        mu_deg,
        sigma,
        passes,
    }
}

/// Run a quick parameter survey and return scores.
pub fn quick_survey() -> Vec<SimulationScore> {
    let grid = ParameterGrid::quick_grid();
    let criteria = SuccessCriteria::paper_defaults();

    grid.generate()
        .iter()
        .map(|p| evaluate_params(p, &criteria))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_paper_grid_size() {
        let grid = ParameterGrid::paper_grid();
        let viable = grid.count_viable();
        // Paper reports ~1,134 viable simulations
        assert!(
            viable > 500 && viable < 2000,
            "Expected ~1134 viable, got {}",
            viable
        );
    }

    #[test]
    fn test_quick_grid() {
        let grid = ParameterGrid::quick_grid();
        let viable = grid.count_viable();
        assert!(viable > 0, "Quick grid should have viable params");
    }

    #[test]
    fn test_evaluate_best_fit() {
        let criteria = SuccessCriteria::paper_defaults();
        let params = P9Params {
            mass_earth: 5.0,
            a: 500.0,
            e: 0.25,
            i: 20.0 * DEG2RAD,
            omega: 140.0 * DEG2RAD,
            omega_big: 100.0 * DEG2RAD,
            mean_anomaly: 0.0,
        };

        let score = evaluate_params(&params, &criteria);
        assert!(score.f_varpi > 0.0);
        assert!(score.eta > 0.0);
    }

    #[test]
    fn test_quick_survey() {
        let scores = quick_survey();
        assert!(!scores.is_empty());
    }
}
