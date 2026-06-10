//! The infinite chain of 2:j resonances with Neptune that shapes the scattered disk.
//!
//! Quadrupole disturbing function (Eq. in the paper):
//! R_q = (Gm_N/4a) * alpha^2 * {X_j^(-3,0) cos(jM) + 3 X_j^(-3,2) cos[jM - 2(M_N - omega)]}
//!
//! Resonance width: Delta_a = 4*a_N*sqrt(2*chi*m_N/(5*M_sun)) * exp(-(q/(2*a_N))^2)

use p9_core::constants::GM_SUN;
use serde::{Deserialize, Serialize};

use crate::hansen::{hansen_x_neg3_2_approx, A_NEPTUNE};

/// Neptune mass in solar masses.
pub const M_NEPTUNE_SOLAR: f64 = 5.15e-5;

/// A single 2:j resonance with Neptune.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NeptuneResonance {
    /// Resonance index j (in n = 2*n_N/j, or equivalently the 2:j resonance)
    pub j: i64,
    /// Nominal semi-major axis (AU)
    pub a_nominal: f64,
    /// Resonance half-width in semi-major axis (AU)
    pub delta_a: f64,
    /// Hansen coefficient X_j^{-3,2}
    pub hansen_coeff: f64,
}

/// Compute the nominal semi-major axis for a 2:j resonance with Neptune.
///
/// From n_particle / n_Neptune = 2/j, so (a_N/a)^{3/2} = 2/j,
/// giving a = a_N * (j/2)^{2/3}.
pub fn resonance_semimajor(j: i64) -> f64 {
    A_NEPTUNE * (j as f64 / 2.0).powf(2.0 / 3.0)
}

/// Resonance half-width (AU) from pendulum model.
///
/// Delta_a = 4*a_N*sqrt(2*chi*m_N/(5*M_sun)) * exp(-(q/(2*a_N))^2)
/// where chi = j/2 for a 2:j resonance.
pub fn resonance_width(j: i64, q: f64) -> f64 {
    let chi = j as f64 / 2.0;
    4.0 * A_NEPTUNE
        * (2.0 * chi * M_NEPTUNE_SOLAR / 5.0).sqrt()
        * (-(q / (2.0 * A_NEPTUNE)).powi(2)).exp()
}

/// Spacing between adjacent 2:j and 2:(j+1) resonances (AU).
pub fn resonance_spacing(j: i64) -> f64 {
    let a_j = resonance_semimajor(j);
    let a_j1 = resonance_semimajor(j + 1);
    (a_j1 - a_j).abs()
}

/// Build the chain of 2:j resonances for j in [j_min, j_max].
pub fn build_resonance_chain(j_min: i64, j_max: i64, q: f64) -> Vec<NeptuneResonance> {
    (j_min..=j_max)
        .map(|j| {
            let a_nominal = resonance_semimajor(j);
            let delta_a = resonance_width(j, q);
            let hansen_coeff = hansen_x_neg3_2_approx(j, q);
            NeptuneResonance {
                j,
                a_nominal,
                delta_a,
                hansen_coeff,
            }
        })
        .collect()
}

/// Isolated resonance pendulum Hamiltonian.
///
/// H_chi = -(3/a_N^2)(chi/2)^{2/3}(Phi_tilde^2/2)
///        - (6*G*m_N/(5*a_N*chi)) * exp(-(q/a_N)^2) * cos(2*phi)
pub fn pendulum_hamiltonian(phi: f64, phi_tilde: f64, j: i64, q: f64) -> f64 {
    let chi = j as f64 / 2.0;
    let kinetic = -3.0 / (A_NEPTUNE * A_NEPTUNE)
        * (chi / 2.0).powf(2.0 / 3.0)
        * (phi_tilde * phi_tilde / 2.0);
    let potential = -(6.0 * GM_SUN * M_NEPTUNE_SOLAR / (5.0 * A_NEPTUNE * chi))
        * (-(q / A_NEPTUNE).powi(2)).exp()
        * (2.0 * phi).cos();
    kinetic + potential
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resonance_semimajor_j2() {
        // 2:2 = 1:1 resonance should be at a_N
        let a = resonance_semimajor(2);
        assert!(
            (a - A_NEPTUNE).abs() < 0.01,
            "2:2 should be at a_N, got {}",
            a
        );
    }

    #[test]
    fn test_resonance_semimajor_increases_with_j() {
        let a10 = resonance_semimajor(10);
        let a20 = resonance_semimajor(20);
        assert!(a20 > a10, "higher j should have larger semi-major axis");
    }

    #[test]
    fn test_resonance_width_positive() {
        let w = resonance_width(20, 35.0);
        assert!(w > 0.0, "resonance width should be positive");
    }

    #[test]
    fn test_resonance_width_decreases_with_q() {
        let w30 = resonance_width(20, 30.0);
        let w40 = resonance_width(20, 40.0);
        assert!(w30 > w40, "width should decrease with perihelion distance");
    }

    #[test]
    fn test_resonance_chain_ordered() {
        let chain = build_resonance_chain(10, 50, 35.0);
        assert_eq!(chain.len(), 41);
        for i in 1..chain.len() {
            assert!(
                chain[i].a_nominal > chain[i - 1].a_nominal,
                "chain should be ordered by a"
            );
        }
    }

    #[test]
    fn test_pendulum_hamiltonian_finite() {
        let h = pendulum_hamiltonian(0.5, 0.1, 20, 35.0);
        assert!(h.is_finite(), "H should be finite, got {}", h);
    }

    #[test]
    fn test_resonance_spacing_positive() {
        let ds = resonance_spacing(20);
        assert!(ds > 0.0, "spacing should be positive");
    }
}
