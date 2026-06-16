//! Stability classification for scattered disk orbits.
//!
//! The conserved action Psi = L - chi*G/2 approximately conserves perihelion q.
//! Classification is based on the Chirikov overlap parameter s(a, q) from
//! `p9_core::analysis::resonance`, so the boundary width scales with the
//! resonance width instead of an arbitrary fixed ±2 AU band:
//!
//! * s > 1 — adjacent 2:j resonances fully overlap: chaotic;
//! * 0.63 < s < 1 — between Chirikov's two-thirds rule (standard-map
//!   K ≈ 0.97, onset of global chaos) and full overlap: marginal;
//! * s < 0.63 — isolated resonances: stable.

use serde::{Deserialize, Serialize};

use p9_core::analysis::resonance::chirikov_overlap_parameter;
use p9_core::constants::{A_NEPTUNE_AU, GM_SUN, MASS_NEPTUNE_SOLAR};
use p9_core::units::{days, Time};

use crate::chirikov::OVERLAP_GLOBAL_CHAOS;

/// Stability classification of a scattered disk orbit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum StabilityClass {
    /// Orbit is in the fully chaotic stochastic layer
    Chaotic,
    /// Orbit is marginally stable (near the boundary)
    Marginal,
    /// Orbit is stable (perihelion well above critical value)
    Stable,
}

/// Classify the stability of an orbit at given (a, q) from the resonance
/// overlap parameter (see module docs for the band definition).
pub fn classify(a: f64, q: f64) -> StabilityClass {
    let s = chirikov_overlap_parameter(a, q);
    if s > 1.0 {
        StabilityClass::Chaotic
    } else if s > OVERLAP_GLOBAL_CHAOS {
        StabilityClass::Marginal
    } else {
        StabilityClass::Stable
    }
}

/// Lyapunov time estimate for chaotic orbits.
///
/// In the stochastic layer, tau_L ~ orbital period = 2pi * sqrt(a^3 / GM_sun).
/// Units: days (using GM_SUN in AU^3/day^2).
pub fn lyapunov_time_days(a: f64) -> f64 {
    std::f64::consts::TAU * (a * a * a / GM_SUN).sqrt()
}

/// Lyapunov time as a dimension-checked [`Time`] for the semi-major axis `a`
/// (in AU). Typed sibling of [`lyapunov_time_days`].
pub fn lyapunov_time(a: f64) -> Time {
    days(lyapunov_time_days(a))
}

/// Semi-major axis diffusion coefficient (AU^2/yr).
///
/// D_a ~ (8/(5*pi)) * (m_N/M_sun) * sqrt(GM_sun * a_N) * exp(-(q/a_N)^2/2)
///
/// Notably independent of semi-major axis.
pub fn diffusion_coefficient(q: f64) -> f64 {
    let v_n = (GM_SUN * A_NEPTUNE_AU).sqrt(); // ~sqrt(GM * 30) in AU/day
    let d_per_day = (8.0 / (5.0 * std::f64::consts::PI))
        * MASS_NEPTUNE_SOLAR
        * v_n
        * (-(q / A_NEPTUNE_AU).powi(2) / 2.0).exp();
    // Convert from AU^2/day to AU^2/yr
    d_per_day * 365.25
}

/// Timescale for semi-major axis diffusion to change a by delta_a (yr).
///
/// tau ~ delta_a^2 / D_a
pub fn diffusion_timescale(q: f64, delta_a: f64) -> f64 {
    let d = diffusion_coefficient(q);
    if d < 1e-30 {
        return f64::INFINITY;
    }
    delta_a * delta_a / d
}

/// Scan a grid of (a, q) values and classify each orbit.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StabilityMap {
    pub a_values: Vec<f64>,
    pub q_values: Vec<f64>,
    pub classes: Vec<Vec<StabilityClass>>,
}

impl StabilityMap {
    pub fn compute(a_min: f64, a_max: f64, q_min: f64, q_max: f64, n_a: usize, n_q: usize) -> Self {
        let da = (a_max - a_min) / (n_a - 1).max(1) as f64;
        let dq = (q_max - q_min) / (n_q - 1).max(1) as f64;

        let a_values: Vec<f64> = (0..n_a).map(|i| a_min + i as f64 * da).collect();
        let q_values: Vec<f64> = (0..n_q).map(|i| q_min + i as f64 * dq).collect();

        let classes: Vec<Vec<StabilityClass>> = a_values
            .iter()
            .map(|&a| q_values.iter().map(|&q| classify(a, q)).collect())
            .collect();

        Self {
            a_values,
            q_values,
            classes,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::analysis::resonance::{critical_perihelion, is_chaotic};
    use uom::si::time::day;

    #[test]
    fn typed_lyapunov_time_matches_f64() {
        assert_relative_eq!(
            lyapunov_time(500.0).get::<day>(),
            lyapunov_time_days(500.0),
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_chaotic_near_neptune() {
        // Need large a for critical perihelion to be nonzero
        let class = classify(500.0, 31.0);
        assert_eq!(class, StabilityClass::Chaotic);
    }

    #[test]
    fn test_stable_high_perihelion() {
        let class = classify(200.0, 60.0);
        assert_eq!(class, StabilityClass::Stable);
    }

    /// H1 regression: `classify` and `is_chaotic` must agree across the
    /// boundary (they previously used overlap forms differing by
    /// exp(q²/4a_N²)).
    #[test]
    fn test_classify_consistent_with_is_chaotic() {
        for &a in &[300.0, 500.0, 1000.0] {
            let q_crit = critical_perihelion(a);
            assert_eq!(classify(a, q_crit - 0.5), StabilityClass::Chaotic);
            assert!(is_chaotic(a, q_crit - 0.5));
            assert_ne!(classify(a, q_crit + 0.5), StabilityClass::Chaotic);
            assert!(!is_chaotic(a, q_crit + 0.5));
        }
    }

    /// The marginal band scales with the resonance width: its q-extent
    /// (between s = 1 and s = 0.63) is a_N² ln(1/0.63)/q, i.e. wider at
    /// small q and narrower at large q.
    #[test]
    fn test_marginal_band_scales_with_resonance_width() {
        let band = |a: f64| {
            let q_crit = critical_perihelion(a);
            // q where s = OVERLAP_GLOBAL_CHAOS, from s ∝ exp(−q²/2a_N²):
            let a_n2 = A_NEPTUNE_AU * A_NEPTUNE_AU;
            let q_outer = (q_crit * q_crit + 2.0 * a_n2 * (1.0 / OVERLAP_GLOBAL_CHAOS).ln()).sqrt();
            (q_crit, q_outer)
        };
        for &a in &[400.0, 800.0] {
            let (q_crit, q_outer) = band(a);
            assert_eq!(
                classify(a, 0.5 * (q_crit + q_outer)),
                StabilityClass::Marginal
            );
            assert_eq!(classify(a, q_outer + 0.5), StabilityClass::Stable);
        }
        // Band narrows as q_crit grows (larger a).
        let (qc1, qo1) = band(400.0);
        let (qc2, qo2) = band(1000.0);
        assert!(qo2 - qc2 < qo1 - qc1);
    }

    #[test]
    fn test_lyapunov_time_positive() {
        let tau = lyapunov_time_days(200.0);
        assert!(tau > 0.0);
        // ~2800 year orbital period at 200 AU
        let tau_yr = tau / 365.25;
        assert!(tau_yr > 1000.0 && tau_yr < 5000.0, "tau = {} yr", tau_yr);
    }

    #[test]
    fn test_diffusion_coefficient_positive() {
        let d = diffusion_coefficient(35.0);
        assert!(d > 0.0, "D_a should be positive");
    }

    #[test]
    fn test_diffusion_decreases_with_q() {
        let d30 = diffusion_coefficient(30.0);
        let d40 = diffusion_coefficient(40.0);
        assert!(d30 > d40, "diffusion should decrease with perihelion");
    }

    #[test]
    fn test_diffusion_timescale() {
        let tau = diffusion_timescale(35.0, 10.0);
        assert!(tau > 0.0 && tau.is_finite());
    }

    #[test]
    fn test_stability_map_dimensions() {
        let map = StabilityMap::compute(100.0, 500.0, 30.0, 50.0, 5, 5);
        assert_eq!(map.a_values.len(), 5);
        assert_eq!(map.q_values.len(), 5);
        assert_eq!(map.classes.len(), 5);
    }
}
