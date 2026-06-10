//! Chirikov resonance overlap criterion for the 2:j chain.
//!
//! Overlap parameter: Delta_a/delta_a = (24/sqrt(5))*(a/a_N)^(5/4) * sqrt(m_N/M_sun) * exp(-(q/(2*a_N))^2)
//!
//! When the overlap parameter exceeds unity, adjacent resonances overlap and chaos ensues.
//! This connects to the standard map K parameter.

use serde::{Deserialize, Serialize};

use crate::hansen::A_NEPTUNE;
use crate::resonance_chain::M_NEPTUNE_SOLAR;

/// Chirikov overlap parameter for the 2:j resonance chain.
///
/// K = Delta_a / delta_a = (24/sqrt(5)) * (a/a_N)^{5/4} * sqrt(m_N/M_sun) * exp(-(q/(2*a_N))^2)
pub fn chirikov_overlap_parameter(a: f64, q: f64) -> f64 {
    let alpha = a / A_NEPTUNE;
    (24.0 / 5.0_f64.sqrt())
        * alpha.powf(1.25)
        * M_NEPTUNE_SOLAR.sqrt()
        * (-(q / (2.0 * A_NEPTUNE)).powi(2)).exp()
}

/// Check if orbits at given (a, q) are in the chaotic overlap regime.
pub fn is_chaotic(a: f64, q: f64) -> bool {
    chirikov_overlap_parameter(a, q) > 1.0
}

/// Standard map stochasticity parameter K from the modulated pendulum.
///
/// The Chirikov standard map emerges as the stroboscopic mapping of the
/// modulated pendulum. K > K_crit ~ 0.97 indicates global chaos.
pub fn standard_map_k(a: f64, q: f64) -> f64 {
    let k = chirikov_overlap_parameter(a, q);
    // The standard map K is the square of the overlap parameter
    // in the modulated pendulum formulation
    k * k
}

/// Compute the critical perihelion distance below which chaos ensues.
///
/// q_crit = a_N * sqrt(ln((24^2/5) * (m_N/M_sun) * (a/a_N)^(5/2)))
///
/// This is the central analytical result of the paper.
pub fn critical_perihelion(a: f64) -> f64 {
    let alpha = a / A_NEPTUNE;
    let argument = (576.0 / 5.0) * M_NEPTUNE_SOLAR * alpha.powf(2.5);
    if argument <= 1.0 {
        return 0.0; // No chaos possible
    }
    A_NEPTUNE * argument.ln().sqrt()
}

/// Scan the (a, q) plane and compute Chirikov overlap parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OverlapGrid {
    pub a_values: Vec<f64>,
    pub q_values: Vec<f64>,
    pub overlap: Vec<Vec<f64>>,
}

impl OverlapGrid {
    /// Generate a grid of overlap parameters.
    pub fn compute(a_min: f64, a_max: f64, q_min: f64, q_max: f64, n_a: usize, n_q: usize) -> Self {
        let da = (a_max - a_min) / (n_a - 1) as f64;
        let dq = (q_max - q_min) / (n_q - 1) as f64;

        let a_values: Vec<f64> = (0..n_a).map(|i| a_min + i as f64 * da).collect();
        let q_values: Vec<f64> = (0..n_q).map(|i| q_min + i as f64 * dq).collect();

        let overlap: Vec<Vec<f64>> = a_values
            .iter()
            .map(|&a| {
                q_values
                    .iter()
                    .map(|&q| chirikov_overlap_parameter(a, q))
                    .collect()
            })
            .collect();

        Self {
            a_values,
            q_values,
            overlap,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_overlap_increases_with_a() {
        let k1 = chirikov_overlap_parameter(100.0, 35.0);
        let k2 = chirikov_overlap_parameter(500.0, 35.0);
        assert!(k2 > k1, "overlap should increase with semi-major axis");
    }

    #[test]
    fn test_overlap_decreases_with_q() {
        let k1 = chirikov_overlap_parameter(200.0, 30.0);
        let k2 = chirikov_overlap_parameter(200.0, 45.0);
        assert!(k1 > k2, "overlap should decrease with perihelion");
    }

    #[test]
    fn test_critical_perihelion_reasonable() {
        // Need large a for the formula argument to exceed 1
        let q_crit = critical_perihelion(500.0);
        assert!(q_crit > 20.0 && q_crit < 60.0, "q_crit = {}", q_crit);
    }

    #[test]
    fn test_critical_perihelion_increases_with_a() {
        let q1 = critical_perihelion(300.0);
        let q2 = critical_perihelion(1000.0);
        assert!(q2 > q1, "q_crit should increase with semi-major axis");
    }

    #[test]
    fn test_is_chaotic_near_neptune() {
        // At large a with perihelion near Neptune, should be chaotic
        assert!(is_chaotic(500.0, 31.0));
    }

    #[test]
    fn test_overlap_grid_dimensions() {
        let grid = OverlapGrid::compute(100.0, 500.0, 30.0, 50.0, 10, 10);
        assert_eq!(grid.a_values.len(), 10);
        assert_eq!(grid.q_values.len(), 10);
        assert_eq!(grid.overlap.len(), 10);
        assert_eq!(grid.overlap[0].len(), 10);
    }
}
