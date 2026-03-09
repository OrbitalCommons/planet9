//! Poincaré action-angle canonical variables for clustering analysis.
//!
//! The paper uses reduced Poincaré variables to properly weight eccentricity
//! and inclination when computing clustering statistics:
//!
//!   Γ = 1 - sqrt(1 - e²),  conjugate angle = -ϖ
//!   Z = sqrt(1 - e²)(1 - cos i),  conjugate angle = -Ω
//!
//! Cartesian projections:
//!   x = sqrt(2Γ) cos(ϖ),  y = sqrt(2Γ) sin(ϖ)   (perihelion clustering)
//!   p = sqrt(2Z) cos(Ω),  q_var = sqrt(2Z) sin(Ω) (pole clustering)

use p9_core::types::OrbitalElements;

/// Poincaré canonical variables for a single KBO.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct PoincareState {
    /// x = sqrt(2Γ) cos(ϖ)
    pub x: f64,
    /// y = sqrt(2Γ) sin(ϖ)
    pub y: f64,
    /// p = sqrt(2Z) cos(Ω)
    pub p: f64,
    /// q = sqrt(2Z) sin(Ω)
    pub q_var: f64,
}

impl PoincareState {
    /// Compute Poincaré variables from orbital elements.
    pub fn from_elements(elem: &OrbitalElements) -> Self {
        let varpi = elem.omega + elem.omega_big;
        let eta = (1.0 - elem.e * elem.e).sqrt();

        // Γ = 1 - sqrt(1 - e²)
        let gamma = 1.0 - eta;
        let sqrt_2gamma = (2.0 * gamma).sqrt();

        // Z = sqrt(1 - e²)(1 - cos i)
        let z = eta * (1.0 - elem.i.cos());
        let sqrt_2z = (2.0 * z).sqrt();

        Self {
            x: sqrt_2gamma * varpi.cos(),
            y: sqrt_2gamma * varpi.sin(),
            p: sqrt_2z * elem.omega_big.cos(),
            q_var: sqrt_2z * elem.omega_big.sin(),
        }
    }
}

/// Compute the mean Poincaré state for a collection of KBOs.
pub fn mean_state(states: &[PoincareState]) -> PoincareState {
    let n = states.len() as f64;
    PoincareState {
        x: states.iter().map(|s| s.x).sum::<f64>() / n,
        y: states.iter().map(|s| s.y).sum::<f64>() / n,
        p: states.iter().map(|s| s.p).sum::<f64>() / n,
        q_var: states.iter().map(|s| s.q_var).sum::<f64>() / n,
    }
}

/// Compute the perihelion clustering distance (distance from origin in x-y plane).
pub fn perihelion_clustering(mean: &PoincareState) -> f64 {
    (mean.x * mean.x + mean.y * mean.y).sqrt()
}

/// Compute the pole clustering distance (distance from origin in p-q plane).
pub fn pole_clustering(mean: &PoincareState) -> f64 {
    (mean.p * mean.p + mean.q_var * mean.q_var).sqrt()
}

/// Compute the 4D clustering distance (combined perihelion + pole).
pub fn combined_clustering(mean: &PoincareState) -> f64 {
    (mean.x * mean.x + mean.y * mean.y + mean.p * mean.p + mean.q_var * mean.q_var).sqrt()
}

/// Mean direction of perihelion clustering (atan2(y, x)).
pub fn mean_varpi_direction(mean: &PoincareState) -> f64 {
    mean.y.atan2(mean.x)
}

/// Mean direction of pole clustering (atan2(q, p)).
pub fn mean_omega_direction(mean: &PoincareState) -> f64 {
    mean.q_var.atan2(mean.p)
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;
    use std::f64::consts::PI;

    #[test]
    fn test_poincare_from_circular() {
        // Circular orbit: e=0 → Γ=0 → x=y=0
        let elem = OrbitalElements {
            a: 300.0,
            e: 0.0,
            i: 10.0 * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        };
        let state = PoincareState::from_elements(&elem);
        assert!(state.x.abs() < 1e-10);
        assert!(state.y.abs() < 1e-10);
    }

    #[test]
    fn test_poincare_from_equatorial() {
        // Equatorial orbit: i=0 → Z=0 → p=q=0
        let elem = OrbitalElements {
            a: 300.0,
            e: 0.8,
            i: 0.0,
            omega: 1.0,
            omega_big: 2.0,
            mean_anomaly: 0.0,
        };
        let state = PoincareState::from_elements(&elem);
        assert!(state.p.abs() < 1e-10);
        assert!(state.q_var.abs() < 1e-10);
    }

    #[test]
    fn test_poincare_eccentric() {
        let elem = OrbitalElements {
            a: 300.0,
            e: 0.85,
            i: 20.0 * DEG2RAD,
            omega: PI / 4.0,
            omega_big: PI / 4.0,
            mean_anomaly: 0.0,
        };
        let state = PoincareState::from_elements(&elem);

        // varpi = pi/4 + pi/4 = pi/2
        // x = sqrt(2Γ) cos(pi/2) ≈ 0
        // y = sqrt(2Γ) sin(pi/2) > 0
        assert!(state.x.abs() < 0.01, "x = {}", state.x);
        assert!(state.y > 0.0, "y should be positive");
    }

    #[test]
    fn test_mean_state() {
        let states = vec![
            PoincareState {
                x: 0.5,
                y: 0.3,
                p: 0.1,
                q_var: 0.2,
            },
            PoincareState {
                x: 0.3,
                y: 0.5,
                p: -0.1,
                q_var: 0.0,
            },
        ];
        let mean = mean_state(&states);
        assert!((mean.x - 0.4).abs() < 1e-10);
        assert!((mean.y - 0.4).abs() < 1e-10);
        assert!((mean.p - 0.0).abs() < 1e-10);
        assert!((mean.q_var - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_clustering_distances() {
        let mean = PoincareState {
            x: 0.3,
            y: 0.4,
            p: 0.0,
            q_var: 0.0,
        };
        assert!((perihelion_clustering(&mean) - 0.5).abs() < 1e-10);
        assert!((pole_clustering(&mean) - 0.0).abs() < 1e-10);
        assert!((combined_clustering(&mean) - 0.5).abs() < 1e-10);
    }
}
