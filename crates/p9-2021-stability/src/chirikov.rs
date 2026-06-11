//! Chirikov resonance overlap criterion for the Neptune 2:j chain.
//!
//! The overlap parameter and the critical perihelion are sourced from
//! `p9_core::analysis::resonance` — the single, internally consistent
//! implementation (this crate previously used exp(−q²/4a_N²) in the overlap
//! parameter while `critical_perihelion` was the K = 1 root of the
//! exp(−q²/2a_N²) form, so `is_chaotic` and `stability::classify` disagreed
//! about the same orbit by a factor exp(q²/4a_N²) ≈ 1.4 at q = 35 AU).

use serde::{Deserialize, Serialize};

pub use p9_core::analysis::resonance::{
    chirikov_overlap_parameter, critical_perihelion, is_chaotic,
};

/// Standard map stochasticity parameter K from the resonance overlap
/// parameter s = Δω/δω.
///
/// Derivation (Chirikov 1979, §4): the standard map's primary resonances sit
/// at integer frequencies (spacing Δω = 2π in action units of the map) and
/// each has full width 4√K, so s = 4√K/2π = 2√K/π and therefore
///
///   K = (π s / 2)².
///
/// Touching resonances (s = 1) give K = π²/4 ≈ 2.47, well above the global
/// chaos threshold K_c ≈ 0.97; conversely K = K_c corresponds to
/// s = 2√K_c/π ≈ 0.63 — Chirikov's "two-thirds rule".
pub fn standard_map_k(a: f64, q: f64) -> f64 {
    let s = chirikov_overlap_parameter(a, q);
    let half_pi_s = std::f64::consts::FRAC_PI_2 * s;
    half_pi_s * half_pi_s
}

/// Overlap parameter s at which the standard map reaches the global chaos
/// threshold K_c ≈ 0.9716: s = 2√K_c/π ≈ 0.63 (the two-thirds rule).
pub const OVERLAP_GLOBAL_CHAOS: f64 = 0.627_465_891_8;

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

    /// H1 regression: the two public APIs must agree about the chaotic
    /// boundary — K(a, q_crit(a)) = 1 exactly.
    #[test]
    fn test_overlap_and_critical_perihelion_consistent() {
        for &a in &[300.0, 500.0, 1000.0] {
            let q_crit = critical_perihelion(a);
            assert!(q_crit > 0.0);
            let k = chirikov_overlap_parameter(a, q_crit);
            assert!((k - 1.0).abs() < 1e-10, "K(a, q_crit) = {k}");
            assert!(is_chaotic(a, q_crit - 1.0));
            assert!(!is_chaotic(a, q_crit + 1.0));
        }
    }

    /// Paper-number regression: BMN21 Eq. for the critical perihelion,
    /// q_crit = a_N √(ln[(24²/5)(m_N/M_☉)(a/a_N)^{5/2}]), evaluated at
    /// a = 500 AU gives 41.4 AU (cf. their Fig. 4 stability boundary).
    #[test]
    fn test_critical_perihelion_paper_value() {
        let q_crit = critical_perihelion(500.0);
        assert!((q_crit - 41.4).abs() < 0.3, "q_crit(500 AU) = {q_crit:.2}");
    }

    #[test]
    fn test_critical_perihelion_increases_with_a() {
        let q1 = critical_perihelion(300.0);
        let q2 = critical_perihelion(1000.0);
        assert!(q2 > q1, "q_crit should increase with semi-major axis");
    }

    #[test]
    fn test_standard_map_k_at_two_thirds_rule() {
        // K = (π s/2)²: s = 1 → K = π²/4; s = OVERLAP_GLOBAL_CHAOS → K ≈ 0.97.
        let s_to_k = |s: f64| (std::f64::consts::FRAC_PI_2 * s).powi(2);
        assert!((s_to_k(1.0) - std::f64::consts::PI.powi(2) / 4.0).abs() < 1e-12);
        assert!((s_to_k(OVERLAP_GLOBAL_CHAOS) - 0.9716).abs() < 1e-3);
        // And the public function is that map of the overlap parameter.
        let (a, q) = (500.0, 35.0);
        let s = chirikov_overlap_parameter(a, q);
        assert!((standard_map_k(a, q) - s_to_k(s)).abs() < 1e-12);
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
