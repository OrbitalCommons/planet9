//! Mean-motion resonance identification and analysis.
//!
//! Key resonances with Planet Nine: 1:1, 3:2, 2:1, 5:2.
//! The critical transition from secular to resonance-dominated dynamics
//! occurs at P/P_9 > 0.1-0.15. Resonance trapping is required for
//! long-period KBO survival over 4 Gyr.

use p9_core::constants::TWO_PI;
use serde::{Deserialize, Serialize};

/// Description of a mean-motion resonance with Planet Nine.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Resonance {
    /// Test particle integer (p in p:q)
    pub p: i64,
    /// Planet Nine integer (q in p:q)
    pub q: i64,
    /// Nominal semi-major axis for this resonance (AU)
    pub a_nominal: f64,
    /// Resonance width in semi-major axis (AU)
    pub delta_a: f64,
}

impl Resonance {
    /// Compute the nominal semi-major axis for a p:q resonance given a_9.
    pub fn nominal_semimajor(p: i64, q: i64, a9: f64) -> f64 {
        a9 * (q as f64 / p as f64).powf(2.0 / 3.0)
    }

    /// Approximate resonance half-width (AU) from pendulum model.
    ///
    /// delta_a ~ a * (m_9/M_sun * alpha)^(1/2) where alpha = a/a_9.
    pub fn approximate_width(a: f64, m9_solar: f64, a9: f64) -> f64 {
        let alpha = a / a9;
        a * (m9_solar * alpha).sqrt()
    }
}

/// Key resonances identified in the paper for fiducial a_9 = 700 AU.
pub fn key_resonances_700() -> Vec<Resonance> {
    let a9 = 700.0;
    let m9_solar = 10.0 * 3.003e-6;

    let resonance_pairs: Vec<(i64, i64)> =
        vec![(1, 1), (3, 2), (2, 1), (5, 2), (3, 1), (4, 1), (5, 1)];

    resonance_pairs
        .into_iter()
        .map(|(p, q)| {
            let a_nominal = Resonance::nominal_semimajor(p, q, a9);
            let delta_a = Resonance::approximate_width(a_nominal, m9_solar, a9);
            Resonance {
                p,
                q,
                a_nominal,
                delta_a,
            }
        })
        .collect()
}

/// Period ratio P/P_9 for a test particle with semi-major axis a.
pub fn period_ratio(a: f64, a9: f64) -> f64 {
    (a / a9).powf(1.5)
}

/// Critical period ratio threshold below which secular dynamics dominate.
///
/// From the paper: P/P_9 ~ 0.1-0.15 marks the transition.
pub const CRITICAL_PERIOD_RATIO: f64 = 0.125;

/// Classify whether a test particle is in the secular or resonance-dominated regime.
pub fn dynamical_regime(a: f64, a9: f64) -> DynamicalRegime {
    let ratio = period_ratio(a, a9);
    if ratio < 0.1 {
        DynamicalRegime::Secular
    } else if ratio < 0.15 {
        DynamicalRegime::Transitional
    } else {
        DynamicalRegime::ResonanceDominated
    }
}

/// Dynamical regime classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DynamicalRegime {
    Secular,
    Transitional,
    ResonanceDominated,
}

/// Check if a test particle with mean motion n is near a p:q resonance with P9
/// (mean motion n9), within a specified tolerance in resonant angle libration.
pub fn is_near_resonance(n: f64, n9: f64, p: i64, q: i64, tol_fraction: f64) -> bool {
    let expected_ratio = p as f64 / q as f64;
    let actual_ratio = n / n9;
    (actual_ratio - expected_ratio).abs() / expected_ratio < tol_fraction
}

/// Adiabatic invariant J = integral of Lambda' d(phi_res').
///
/// Under the J=0 approximation (negligible libration amplitude),
/// objects are trapped at the exact center of resonance.
/// Returns the resonant semi-major axis deviation from nominal.
pub fn adiabatic_deviation(libration_amplitude: f64, resonance_width: f64) -> f64 {
    libration_amplitude * resonance_width / TWO_PI
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nominal_semimajor_1_1() {
        let a = Resonance::nominal_semimajor(1, 1, 700.0);
        assert!((a - 700.0).abs() < 0.01, "1:1 resonance should be at a_9");
    }

    #[test]
    fn test_nominal_semimajor_2_1() {
        // 2:1 means n_particle = 2 * n_9, so a_particle = a_9 * (1/2)^(2/3)
        let a = Resonance::nominal_semimajor(2, 1, 700.0);
        let expected = 700.0 * (0.5_f64).powf(2.0 / 3.0);
        assert!(
            (a - expected).abs() < 0.1,
            "a = {}, expected {}",
            a,
            expected
        );
        assert!(a < 700.0, "2:1 resonance should be inside P9 orbit");
    }

    #[test]
    fn test_key_resonances_ordered() {
        let res = key_resonances_700();
        assert!(res.len() >= 5, "should have at least 5 key resonances");
        // 1:1 should have largest a_nominal
        assert!(
            res[0].a_nominal > res[1].a_nominal,
            "1:1 should be outermost"
        );
    }

    #[test]
    fn test_period_ratio() {
        let ratio = period_ratio(700.0, 700.0);
        assert!((ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dynamical_regime_secular() {
        // a = 100 AU, a9 = 700 AU => P/P9 = (100/700)^1.5 ~ 0.027
        let regime = dynamical_regime(100.0, 700.0);
        assert_eq!(regime, DynamicalRegime::Secular);
    }

    #[test]
    fn test_dynamical_regime_resonance() {
        // a = 500 AU, a9 = 700 AU => P/P9 = (500/700)^1.5 ~ 0.61
        let regime = dynamical_regime(500.0, 700.0);
        assert_eq!(regime, DynamicalRegime::ResonanceDominated);
    }

    #[test]
    fn test_resonance_width_positive() {
        let dw = Resonance::approximate_width(400.0, 10.0 * 3.003e-6, 700.0);
        assert!(dw > 0.0, "resonance width should be positive");
    }

    #[test]
    fn test_is_near_resonance() {
        let n9 = 1.0;
        let n = 2.0; // exactly 2:1
        assert!(is_near_resonance(n, n9, 2, 1, 0.01));
        assert!(!is_near_resonance(n, n9, 3, 2, 0.01));
    }
}
