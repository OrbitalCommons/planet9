//! Three-regime stability boundary computation (Belyakov & Batygin 2025).
//!
//! Regime 1 (fully chaotic layer):     q = 13.37 + 0.55 * a^0.64
//! Regime 2 (1:j diffusion barrier):   q = -63.1 + 19.4 * ln(a)
//! Regime 3 (fine comb):               q = 10.28 + 6.01 * a^0.34
//!
//! Semantics per the paper's figure: the *instability* boundary is set by
//! the fully chaotic layer at small a and by the 1:j diffusion barrier
//! beyond their crossing (a ≈ 90 AU) — NOT by the fine comb. The fine
//! comb of high-order intermediate resonances above the diffusion barrier
//! produces only weak, *bounded* chaos: orbits there are chaotic but do
//! not diffuse to instability, so taking max() of all three curves (as
//! the previous version did) wrongly made the fine comb the binding
//! instability boundary over a ≈ 100–300 AU.
//!
//! Natural-log check for the −63.1 + 19.4·log(a) fit: with log = ln the
//! curve gives q ≈ 40 AU at a = 200 AU and stays positive for a ≳ 26 AU,
//! matching the paper's figure; with log10 it would be negative until
//! a ≈ 1800 AU (q(200) ≈ −18 AU), which is unphysical — the fit uses ln.

use serde::{Deserialize, Serialize};

/// Stability boundary regime.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BoundaryRegime {
    /// Inner chaotic layer dominated by strong resonance overlap
    ChaoticLayer,
    /// Intermediate barrier set by 1:j resonance diffusion
    DiffusionBarrier,
    /// Outer fine comb of higher-order resonances (weak, bounded chaos —
    /// not an instability boundary)
    FineComb,
}

/// Orbit classification against the three regimes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OrbitClass {
    /// Below the instability boundary: macroscopic chaotic diffusion
    Unstable,
    /// Between the instability boundary and the fine comb: chaotic but
    /// bounded (no large-scale diffusion)
    BoundedChaos,
    /// Above the fine comb: regular
    Stable,
}

/// Fully chaotic layer boundary.
///
/// q = 13.37 + 0.55 * a^0.64
pub fn chaotic_layer_boundary(a: f64) -> f64 {
    13.37 + 0.55 * a.powf(0.64)
}

/// 1:j resonance diffusion barrier boundary.
///
/// q = -63.1 + 19.4 * ln(a)   (natural log — see module docs)
pub fn diffusion_barrier_boundary(a: f64) -> f64 {
    -63.1 + 19.4 * a.ln()
}

/// Fine comb boundary from higher-order resonance overlap. Marks the
/// onset of weak *bounded* chaos, not instability.
///
/// q = 10.28 + 6.01 * a^0.34
pub fn fine_comb_boundary(a: f64) -> f64 {
    10.28 + 6.01 * a.powf(0.34)
}

/// The instability boundary at a given semi-major axis, with the regime
/// that sets it: the chaotic layer below the crossing at a ≈ 90 AU, the
/// 1:j diffusion barrier beyond. The fine comb is deliberately excluded
/// (bounded chaos — see module docs).
pub fn effective_boundary(a: f64) -> (f64, BoundaryRegime) {
    let q_chaotic = chaotic_layer_boundary(a);
    let q_diffusion = diffusion_barrier_boundary(a);

    if q_chaotic >= q_diffusion {
        (q_chaotic, BoundaryRegime::ChaoticLayer)
    } else {
        (q_diffusion, BoundaryRegime::DiffusionBarrier)
    }
}

/// Classify orbit stability using the piecewise boundary semantics.
pub fn classify_orbit(a: f64, q: f64) -> OrbitClass {
    let (q_instability, _) = effective_boundary(a);
    if q <= q_instability {
        OrbitClass::Unstable
    } else if q <= fine_comb_boundary(a) {
        OrbitClass::BoundedChaos
    } else {
        OrbitClass::Stable
    }
}

/// Whether an orbit avoids macroscopic instability (bounded chaos in the
/// fine comb counts as surviving).
pub fn is_stable(a: f64, q: f64) -> bool {
    classify_orbit(a, q) != OrbitClass::Unstable
}

/// Compute the stability boundary curves over a range of semi-major axes.
/// All three regime curves are exposed; `q_effective`/`regimes` carry the
/// instability boundary (chaotic layer / diffusion barrier).
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
    fn test_diffusion_barrier_uses_natural_log() {
        // ln: q(200) ≈ 39.7 AU (paper figure); log10 would give −18.5.
        let q = diffusion_barrier_boundary(200.0);
        assert!((q - 39.7).abs() < 0.5, "q = {q}");
        // Positive across the paper's fitted range.
        for a in [50.0, 100.0, 500.0, 1000.0] {
            assert!(diffusion_barrier_boundary(a) > 0.0, "negative at a = {a}");
        }
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
    fn test_regime_handoff_piecewise() {
        // Chaotic layer binds at small a, diffusion barrier beyond ~90 AU.
        let (_, regime_inner) = effective_boundary(60.0);
        assert_eq!(regime_inner, BoundaryRegime::ChaoticLayer);
        let (_, regime_outer) = effective_boundary(200.0);
        assert_eq!(regime_outer, BoundaryRegime::DiffusionBarrier);
    }

    #[test]
    fn test_fine_comb_not_instability_boundary() {
        // At a = 200 AU the fine comb (≈ 46.8 AU) lies above the
        // diffusion barrier (≈ 39.7 AU): with the old max() rule it set
        // the boundary; here q = 43 AU is bounded chaos, not unstable.
        let a = 200.0;
        assert!(fine_comb_boundary(a) > effective_boundary(a).0);
        assert_eq!(classify_orbit(a, 43.0), OrbitClass::BoundedChaos);
        assert!(is_stable(a, 43.0));
        assert_eq!(classify_orbit(a, 50.0), OrbitClass::Stable);
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
        assert_eq!(classify_orbit(200.0, 25.0), OrbitClass::Unstable);
    }

    #[test]
    fn test_boundary_curve_dimensions() {
        let curve = BoundaryCurve::compute(100.0, 500.0, 50);
        assert_eq!(curve.a_values.len(), 50);
        assert_eq!(curve.q_effective.len(), 50);
        assert_eq!(curve.regimes.len(), 50);
        // The effective boundary never exceeds both source curves and the
        // fine comb is exposed separately.
        for i in 0..50 {
            assert!(curve.q_effective[i] <= curve.q_chaotic[i].max(curve.q_diffusion[i]) + 1e-12);
            assert!(curve.q_comb[i] > 0.0);
        }
    }
}
