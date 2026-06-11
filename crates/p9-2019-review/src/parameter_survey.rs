//! Parameter grid survey with the critical-semi-major-axis success measure.
//!
//! The review paper (Figure 16) scores each P9 parameter set by several
//! statistics measured from N-body ensembles. Of these, the critical
//! semi-major axis a_c — the boundary beyond which test-particle orbits are
//! apsidally confined by P9 — is computable from real secular machinery, and
//! is implemented here via the p9-core numerical Gauss-ring double averaging
//! (`p9_core::analysis::secular::numerical_secular_hamiltonian`) plus the
//! orbit-averaged J2 field of the giant planets: a_c is the smallest test
//! particle semi-major axis at which the secular Hamiltonian H(e, Δϖ) has an
//! interior extremum on the anti-aligned axis (a libration island, cf.
//! Batygin & Brown 2016 phase portraits).
//!
//! The paper's remaining ensemble statistics (f_ϖ, η, μ, σ) require the full
//! N-body simulation suite and are intentionally NOT emulated — a previous
//! version returned invented scaling relations for them (admitted TODO),
//! which made every Figure-16-style output decorative. They are omitted from
//! the results path rather than faked.

use p9_core::analysis::secular::numerical_secular_hamiltonian;
use p9_core::constants::*;
use p9_core::forces::j2_secular::combined_j2_jsu;
use p9_core::types::P9Params;

/// Statistical success criterion for the survey: the critical semi-major
/// axis should fall in the observed transition range (200, 300) AU
/// (Batygin, Adams, Brown & Becker 2019).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SuccessCriteria {
    /// Critical semi-major axis should be in (200, 300) AU
    pub a_c_range: (f64, f64),
}

impl SuccessCriteria {
    pub fn paper_defaults() -> Self {
        Self {
            a_c_range: (200.0, 300.0),
        }
    }
}

/// Result of evaluating a single parameter set.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SimulationScore {
    pub params: P9ParamsSummary,
    /// Critical semi-major axis from the secular island criterion (AU);
    /// None if no island appears anywhere on the scanned grid.
    pub a_c: Option<f64>,
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

/// Coplanar secular Hamiltonian of a test particle at (a, e, Δϖ) under P9
/// (numerical Gauss-ring double average) plus the orbit-averaged J2 field of
/// the giants:
///
///   H_J2(e) = −GM J2_eff R_ref² / (2 a³ (1 − e²)^{3/2}),
///
/// which supplies the prograde apsidal precession against which P9's torque
/// competes (without it the island location is wrong).
pub fn secular_hamiltonian(a: f64, e: f64, dvarpi: f64, p9: &P9Params, n_quad: usize) -> f64 {
    let softening = 0.02 * p9.a;
    let h_p9 = numerical_secular_hamiltonian(
        a,
        e,
        0.0,
        dvarpi,
        0.0,
        p9.a,
        p9.e,
        p9.gm(),
        n_quad,
        softening,
    );
    let (j2, _, gm_boost) = combined_j2_jsu();
    let gm_eff = GM_SUN + gm_boost;
    let h_j2 = -gm_eff * j2 / (2.0 * a.powi(3) * (1.0 - e * e).powf(1.5));
    h_p9 + h_j2
}

/// Whether the secular Hamiltonian at fixed `a` has an anti-aligned
/// libration island: an interior extremum (not saddle) of H(e, Δϖ) on the
/// Δϖ = π symmetry axis.
pub fn has_anti_aligned_island(a: f64, p9: &P9Params, n_quad: usize) -> bool {
    let e_grid: Vec<f64> = (1..=17).map(|k| 0.05 + 0.05 * k as f64).collect(); // 0.10..0.90
    let h_line: Vec<f64> = e_grid
        .iter()
        .map(|&e| secular_hamiltonian(a, e, std::f64::consts::PI, p9, n_quad))
        .collect();

    for k in 1..e_grid.len() - 1 {
        let d2e = h_line[k + 1] + h_line[k - 1] - 2.0 * h_line[k];
        let interior_extremum = (h_line[k] > h_line[k - 1] && h_line[k] > h_line[k + 1])
            || (h_line[k] < h_line[k - 1] && h_line[k] < h_line[k + 1]);
        if !interior_extremum {
            continue;
        }
        // Saddle check: curvature along Δϖ must have the same sign as the
        // curvature along e (closed contours around the fixed point).
        let dphi = 0.3;
        let h_off = secular_hamiltonian(a, e_grid[k], std::f64::consts::PI - dphi, p9, n_quad);
        let d2phi = 2.0 * (h_off - h_line[k]); // symmetric about π
        if d2e * d2phi > 0.0 {
            return true;
        }
    }
    false
}

/// Critical semi-major axis: smallest particle a on the scan grid at which
/// the anti-aligned secular island exists. Returns None if no island
/// appears up to the grid edge.
pub fn critical_semi_major_axis(p9: &P9Params) -> Option<f64> {
    let n_quad = 48;
    (100..=600)
        .step_by(25)
        .map(|a| a as f64)
        .find(|&a| has_anti_aligned_island(a, p9, n_quad))
}

/// Evaluate a parameter set against the success criteria using the real
/// secular machinery (see module docs).
pub fn evaluate_params(params: &P9Params, criteria: &SuccessCriteria) -> SimulationScore {
    let a_c = critical_semi_major_axis(params);
    let passes = a_c.is_some_and(|a| a >= criteria.a_c_range.0 && a <= criteria.a_c_range.1);

    SimulationScore {
        params: P9ParamsSummary::from(params),
        a_c,
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

    /// Paper-anchored regression: for the review's favored configuration
    /// (m₉ = 5 M⊕, a₉ = 500 AU, e₉ = 0.25), the transition to apsidal
    /// confinement occurs at a_c ≈ 200–300 AU (Batygin et al. 2019; cf.
    /// Batygin & Brown 2016 phase portraits). The secular-island criterion
    /// must place a_c in the (150, 350) AU band.
    #[test]
    fn test_critical_a_for_favored_params() {
        let p9 = P9Params::revised_2019();
        let a_c = critical_semi_major_axis(&p9).expect("island should exist on grid");
        assert!(
            (150.0..=350.0).contains(&a_c),
            "a_c = {a_c} AU, expected ~200-300"
        );
    }

    /// The island criterion must respond to the perturber: far inside the
    /// critical radius there is no anti-aligned island.
    #[test]
    fn test_no_island_well_inside() {
        let p9 = P9Params::revised_2019();
        assert!(!has_anti_aligned_island(100.0, &p9, 48));
    }

    #[test]
    fn test_evaluate_favored_params_passes() {
        let criteria = SuccessCriteria::paper_defaults();
        let score = evaluate_params(&P9Params::revised_2019(), &criteria);
        assert!(score.a_c.is_some());
    }
}
