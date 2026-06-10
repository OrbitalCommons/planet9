//! Stability classification for scattered disk orbits.
//!
//! The conserved action Psi = L - chi*G/2 approximately conserves perihelion q.
//! Objects below q_crit are unstable on timescales comparable to their orbital period.

use serde::{Deserialize, Serialize};

use crate::chirikov::critical_perihelion;
use crate::hansen::A_NEPTUNE;
use crate::resonance_chain::M_NEPTUNE_SOLAR;

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

/// Classify the stability of an orbit at given (a, q).
pub fn classify(a: f64, q: f64) -> StabilityClass {
    let q_crit = critical_perihelion(a);
    if q < q_crit - 2.0 {
        StabilityClass::Chaotic
    } else if q < q_crit + 2.0 {
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
    use p9_core::constants::GM_SUN;
    std::f64::consts::TAU * (a * a * a / GM_SUN).sqrt()
}

/// Semi-major axis diffusion coefficient (AU^2/yr).
///
/// D_a ~ (8/(5*pi)) * (m_N/M_sun) * sqrt(GM_sun * a_N) * exp(-(q/a_N)^2/2)
///
/// Notably independent of semi-major axis.
pub fn diffusion_coefficient(q: f64) -> f64 {
    use p9_core::constants::GM_SUN;
    let v_n = (GM_SUN * A_NEPTUNE).sqrt(); // ~sqrt(GM * 30) in AU/day
    let d_per_day = (8.0 / (5.0 * std::f64::consts::PI))
        * M_NEPTUNE_SOLAR
        * v_n
        * (-(q / A_NEPTUNE).powi(2) / 2.0).exp();
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
