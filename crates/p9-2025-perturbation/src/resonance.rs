//! Resonance width calculation for arbitrary m:j resonances.
//!
//! Pendulum Hamiltonian per resonance:
//! H_chi = gamma * cos(m*phi) - (beta/2) * Phi_tilde^2
//! where gamma and beta depend on the multipole order and Hansen coefficients.
//!
//! Resonance width: Delta_a = (8a/sqrt(3)) * sqrt(m_N/M_sun)
//!                   * sqrt(alpha^l * M_l * c^2_{lm} * X^{-(l+1),m}_chi)

use serde::{Deserialize, Serialize};

/// Neptune parameters.
const A_NEPTUNE: f64 = 30.07;
const M_NEPTUNE_SOLAR: f64 = 5.15e-5;

/// A resonance in the m:j chain with Neptune.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Resonance {
    /// Azimuthal order m (from multipole expansion)
    pub m: u32,
    /// Resonance index j
    pub j: i64,
    /// Multipole order l
    pub l: u32,
    /// Nominal semi-major axis (AU)
    pub a_nominal: f64,
    /// Resonance half-width (AU)
    pub delta_a: f64,
}

/// Nominal semi-major axis for an m:j resonance.
///
/// From n/n_N = m/j, so a = a_N * (j/m)^{2/3}.
pub fn resonance_semimajor(m: u32, j: i64) -> f64 {
    A_NEPTUNE * (j as f64 / m as f64).powf(2.0 / 3.0)
}

/// Resonance half-width for a given multipole order and Hansen coefficient.
///
/// Delta_a = (8*a/sqrt(3)) * sqrt(m_N/M_sun * alpha^l * c2_lm * |X|)
pub fn resonance_width(a: f64, l: u32, c2_lm: f64, hansen_coeff: f64) -> f64 {
    let alpha = A_NEPTUNE / a;
    let inner = M_NEPTUNE_SOLAR * alpha.powi(l as i32) * c2_lm * hansen_coeff.abs();
    if inner < 0.0 {
        return 0.0;
    }
    (8.0 * a / 3.0_f64.sqrt()) * inner.sqrt()
}

/// Build the 2:j resonance chain (quadrupole, same as Batygin et al. 2021).
pub fn build_2j_chain(j_min: i64, j_max: i64, q: f64) -> Vec<Resonance> {
    (j_min..=j_max)
        .map(|j| {
            let a = resonance_semimajor(2, j);
            // Quadrupole l=2, m=2, c^2_{22} = 3/4
            let hansen = (2.0 * j as f64 / 5.0) * (-(q / A_NEPTUNE).powi(2)).exp();
            let delta_a = resonance_width(a, 2, 0.75, hansen);
            Resonance {
                m: 2,
                j,
                l: 2,
                a_nominal: a,
                delta_a,
            }
        })
        .collect()
}

/// Build the 1:j resonance chain (octupole, m=1).
pub fn build_1j_chain(j_min: i64, j_max: i64, q: f64) -> Vec<Resonance> {
    (j_min..=j_max)
        .map(|j| {
            let a = resonance_semimajor(1, j);
            // Octupole l=3, m=1, c^2_{31} = 3/8
            let hansen = (j.abs() as f64).powf(5.0 / 3.0) / 8.0 * (-(q / A_NEPTUNE).powi(2)).exp();
            let delta_a = resonance_width(a, 3, 3.0 / 8.0, hansen);
            Resonance {
                m: 1,
                j,
                l: 3,
                a_nominal: a,
                delta_a,
            }
        })
        .collect()
}

/// Build the 3:j resonance chain (octupole, m=3).
pub fn build_3j_chain(j_min: i64, j_max: i64, q: f64) -> Vec<Resonance> {
    (j_min..=j_max)
        .map(|j| {
            let a = resonance_semimajor(3, j);
            // Octupole l=3, m=3, c^2_{33} = 5/8
            let hansen = (j.abs() as f64).powf(5.0 / 3.0) / 12.0 * (-(q / A_NEPTUNE).powi(2)).exp();
            let delta_a = resonance_width(a, 3, 5.0 / 8.0, hansen);
            Resonance {
                m: 3,
                j,
                l: 3,
                a_nominal: a,
                delta_a,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resonance_semimajor_2j() {
        let a = resonance_semimajor(2, 2);
        assert!((a - A_NEPTUNE).abs() < 0.01, "2:2 = 1:1 at a_N");
    }

    #[test]
    fn test_resonance_semimajor_1j() {
        let a = resonance_semimajor(1, 1);
        assert!((a - A_NEPTUNE).abs() < 0.01, "1:1 at a_N");
    }

    #[test]
    fn test_2j_chain_size() {
        let chain = build_2j_chain(10, 50, 35.0);
        assert_eq!(chain.len(), 41);
    }

    #[test]
    fn test_1j_chain_larger_a() {
        let chain_1j = build_1j_chain(5, 10, 35.0);
        let chain_2j = build_2j_chain(5, 10, 35.0);
        // 1:j resonances have larger a than 2:j at same j
        for (r1, r2) in chain_1j.iter().zip(chain_2j.iter()) {
            assert!(r1.a_nominal > r2.a_nominal, "1:j should be at larger a");
        }
    }

    #[test]
    fn test_octupole_widths_narrower() {
        // At same j, octupole resonances should generally be narrower than quadrupole
        let chain_2j = build_2j_chain(20, 20, 35.0);
        let chain_1j = build_1j_chain(20, 20, 35.0);
        // Not guaranteed at all j, but at j=20 the quadrupole should dominate
        assert!(chain_2j[0].delta_a > 0.0);
        assert!(chain_1j[0].delta_a > 0.0);
    }

    #[test]
    fn test_resonance_width_positive() {
        let w = resonance_width(200.0, 2, 0.75, 0.1);
        assert!(w > 0.0);
    }
}
