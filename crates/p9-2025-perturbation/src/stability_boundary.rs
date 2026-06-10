//! Three-regime stability boundary computation.
//!
//! Regime 1 (fully chaotic layer):     q = 13.37 + 0.55 * a^0.64
//! Regime 2 (1:j diffusion barrier):   q = -63.1 + 19.4 * log(a)
//! Regime 3 (fine comb boundary):      q = 10.28 + 6.01 * a^0.34

use serde::{Deserialize, Serialize};

/// Stability boundary regime.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BoundaryRegime {
    /// Inner chaotic layer dominated by strong resonance overlap
    ChaoticLayer,
    /// Intermediate barrier set by 1:j resonance diffusion
    DiffusionBarrier,
    /// Outer fine comb of higher-order resonances
    FineComb,
}

/// Fully chaotic layer boundary.
///
/// q = 13.37 + 0.55 * a^0.64
pub fn chaotic_layer_boundary(a: f64) -> f64 {
    13.37 + 0.55 * a.powf(0.64)
}

/// 1:j resonance diffusion barrier boundary.
///
/// q = -63.1 + 19.4 * ln(a)
pub fn diffusion_barrier_boundary(a: f64) -> f64 {
    -63.1 + 19.4 * a.ln()
}

/// Fine comb boundary from higher-order resonance overlap.
///
/// q = 10.28 + 6.01 * a^0.34
pub fn fine_comb_boundary(a: f64) -> f64 {
    10.28 + 6.01 * a.powf(0.34)
}

/// Determine which regime applies at a given semi-major axis.
///
/// The effective boundary is the maximum of all three regimes.
pub fn effective_boundary(a: f64) -> (f64, BoundaryRegime) {
    let q_chaotic = chaotic_layer_boundary(a);
    let q_diffusion = diffusion_barrier_boundary(a);
    let q_comb = fine_comb_boundary(a);

    if q_chaotic >= q_diffusion && q_chaotic >= q_comb {
        (q_chaotic, BoundaryRegime::ChaoticLayer)
    } else if q_diffusion >= q_comb {
        (q_diffusion, BoundaryRegime::DiffusionBarrier)
    } else {
        (q_comb, BoundaryRegime::FineComb)
    }
}

/// Classify orbit stability using the three-regime boundary.
pub fn is_stable(a: f64, q: f64) -> bool {
    let (q_boundary, _) = effective_boundary(a);
    q > q_boundary
}

/// Compute the stability boundary curve over a range of semi-major axes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryCurve {
    pub a_values: Vec<f64>,
    pub q_chaotic: Vec<f64>,
    pub q_diffusion: Vec<f64>,
    pub q_comb: Vec<f64>,
    pub q_effective: Vec<f64>,
    pub regimes: Vec<BoundaryRegime>,
}

impl BoundaryCurve {
    pub fn compute(a_min: f64, a_max: f64, n_points: usize) -> Self {
        let da = (a_max - a_min) / (n_points - 1).max(1) as f64;
        let mut curve = Self {
            a_values: Vec::with_capacity(n_points),
            q_chaotic: Vec::with_capacity(n_points),
            q_diffusion: Vec::with_capacity(n_points),
            q_comb: Vec::with_capacity(n_points),
            q_effective: Vec::with_capacity(n_points),
            regimes: Vec::with_capacity(n_points),
        };

        for i in 0..n_points {
            let a = a_min + i as f64 * da;
            let (q_eff, regime) = effective_boundary(a);
            curve.a_values.push(a);
            curve.q_chaotic.push(chaotic_layer_boundary(a));
            curve.q_diffusion.push(diffusion_barrier_boundary(a));
            curve.q_comb.push(fine_comb_boundary(a));
            curve.q_effective.push(q_eff);
            curve.regimes.push(regime);
        }

        curve
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chaotic_layer_reasonable() {
        let q = chaotic_layer_boundary(200.0);
        assert!(q > 20.0 && q < 60.0, "q = {}", q);
    }

    #[test]
    fn test_diffusion_barrier_reasonable() {
        let q = diffusion_barrier_boundary(200.0);
        assert!(q > 20.0 && q < 60.0, "q = {}", q);
    }

    #[test]
    fn test_fine_comb_reasonable() {
        let q = fine_comb_boundary(200.0);
        assert!(q > 20.0 && q < 60.0, "q = {}", q);
    }

    #[test]
    fn test_boundaries_increase_with_a() {
        let q1 = chaotic_layer_boundary(100.0);
        let q2 = chaotic_layer_boundary(500.0);
        assert!(q2 > q1, "boundary should increase with a");
    }

    #[test]
    fn test_effective_boundary_is_max() {
        let a = 200.0;
        let (q_eff, _) = effective_boundary(a);
        assert!(q_eff >= chaotic_layer_boundary(a));
        assert!(q_eff >= diffusion_barrier_boundary(a));
        assert!(q_eff >= fine_comb_boundary(a));
    }

    #[test]
    fn test_is_stable_high_q() {
        assert!(is_stable(200.0, 60.0), "q=60 AU should be stable at a=200");
    }

    #[test]
    fn test_is_unstable_low_q() {
        assert!(
            !is_stable(200.0, 25.0),
            "q=25 AU should be unstable at a=200"
        );
    }

    #[test]
    fn test_boundary_curve_dimensions() {
        let curve = BoundaryCurve::compute(100.0, 500.0, 50);
        assert_eq!(curve.a_values.len(), 50);
        assert_eq!(curve.q_effective.len(), 50);
        assert_eq!(curve.regimes.len(), 50);
    }
}
