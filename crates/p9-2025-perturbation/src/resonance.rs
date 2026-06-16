//! Resonance width calculation for arbitrary m:j resonances.
//!
//! Pendulum Hamiltonian per resonance:
//! H_chi = gamma * cos(m*phi) - (beta/2) * Phi_tilde^2
//! where gamma and beta depend on the multipole order and Hansen coefficients.
//!
//! Resonance half-width (Belyakov & Batygin 2025, pendulum width):
//!
//!   Delta_a = (8a/sqrt(3)) * sqrt(m_N/M_sun * alpha^l * c^2_{lm} * |X^{-(l+1),m}_j|)
//!
//! Constants come from `p9_core::constants` (no local Neptune
//! redefinitions); the 2:j chain uses the canonical asymptotic
//! `p9_core::analysis::hansen::hansen_x_neg3_2_neptune_chain`, which makes
//! these widths *exactly* consistent with the p9-core Chirikov overlap
//! parameter (see the cross-crate test in `chirikov`); the 1:j, 3:j and
//! 4:j chains evaluate the Hansen coefficients numerically at the
//! resonance geometry e = 1 − q/a.

use p9_core::analysis::hansen::{hansen_coefficient, hansen_x_neg3_2_neptune_chain};
use p9_core::analysis::resonance::resonance_semi_major_axis;
use p9_core::constants::{A_NEPTUNE_AU, MASS_NEPTUNE_SOLAR};
use p9_core::units::{au, Length};
use serde::{Deserialize, Serialize};

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

impl Resonance {
    /// Nominal semi-major axis as a dimension-checked [`Length`] (AU storage).
    pub fn semi_major_axis(&self) -> Length {
        au(self.a_nominal)
    }

    /// Resonance half-width as a dimension-checked [`Length`] (AU storage).
    pub fn half_width(&self) -> Length {
        au(self.delta_a)
    }
}

/// Nominal semi-major axis for an m:j resonance (n/n_N = m/j):
/// a = a_N (j/m)^{2/3}, via the shared p9-core implementation.
pub fn resonance_semimajor(m: u32, j: i64) -> f64 {
    assert!(j > 0, "resonance index must be positive, got {j}");
    resonance_semi_major_axis(j as u32, m, A_NEPTUNE_AU)
}

/// Resonance half-width for a given multipole order and Hansen coefficient.
///
/// Delta_a = (8*a/sqrt(3)) * sqrt(m_N/M_sun * alpha^l * c2_lm * |X|)
pub fn resonance_width(a: f64, l: u32, c2_lm: f64, hansen_coeff: f64) -> f64 {
    let alpha = A_NEPTUNE_AU / a;
    let inner = MASS_NEPTUNE_SOLAR * alpha.powi(l as i32) * c2_lm * hansen_coeff.abs();
    (8.0 * a / 3.0_f64.sqrt()) * inner.sqrt()
}

/// Eccentricity of an orbit with semi-major axis `a` and perihelion `q`,
/// clamped to [0, 1): resonances inside q (a < q) get e = 0, for which
/// the high-j Hansen coefficients vanish identically.
fn eccentricity_at(a: f64, q: f64) -> f64 {
    (1.0 - q / a).max(0.0)
}

/// Numerical-Hansen tolerance for chain construction.
const CHAIN_TOL: f64 = 1.0e-6;

/// Build the 2:j resonance chain (quadrupole, same as Batygin et al. 2021).
///
/// Uses the canonical p9-core asymptotic X_j^{-3,2} ≈ (2j/5)·e^{−(q/a_N)²}
/// — the single sourced form (this fixes the previous three mutually
/// inconsistent exponents across the workspace).
pub fn build_2j_chain(j_min: i64, j_max: i64, q: f64) -> Vec<Resonance> {
    (j_min..=j_max)
        .map(|j| {
            let a = resonance_semimajor(2, j);
            // Quadrupole l=2, m=2, c^2_{22} = 3/4
            let hansen = hansen_x_neg3_2_neptune_chain(j, q);
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

/// Build the 1:j resonance chain (octupole, m=1), with X_j^{-4,1}
/// evaluated numerically at the resonance geometry.
pub fn build_1j_chain(j_min: i64, j_max: i64, q: f64) -> Vec<Resonance> {
    (j_min..=j_max)
        .map(|j| {
            let a = resonance_semimajor(1, j);
            let e = eccentricity_at(a, q);
            // Octupole l=3, m=1, c^2_{31} = 3/8
            let hansen = hansen_coefficient(j, -4, 1, e, CHAIN_TOL);
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

/// Build the 3:j resonance chain (octupole, m=3), with X_j^{-4,3}
/// evaluated numerically at the resonance geometry.
pub fn build_3j_chain(j_min: i64, j_max: i64, q: f64) -> Vec<Resonance> {
    (j_min..=j_max)
        .map(|j| {
            let a = resonance_semimajor(3, j);
            let e = eccentricity_at(a, q);
            // Octupole l=3, m=3, c^2_{33} = 5/8
            let hansen = hansen_coefficient(j, -4, 3, e, CHAIN_TOL);
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

/// Build the 4:j resonance chain (hexadecapole l=4, m=4) promised by the
/// crate docs: c^2_{44} = 35/64 from the P₄ azimuthal decomposition (see
/// `spherical_harmonics::hexadecapole_coefficients`), with X_j^{-5,4}
/// evaluated numerically.
pub fn build_4j_chain(j_min: i64, j_max: i64, q: f64) -> Vec<Resonance> {
    (j_min..=j_max)
        .map(|j| {
            let a = resonance_semimajor(4, j);
            let e = eccentricity_at(a, q);
            let hansen = hansen_coefficient(j, -5, 4, e, CHAIN_TOL);
            let delta_a = resonance_width(a, 4, 35.0 / 64.0, hansen);
            Resonance {
                m: 4,
                j,
                l: 4,
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
        assert!((a - A_NEPTUNE_AU).abs() < 0.01, "2:2 = 1:1 at a_N");
    }

    #[test]
    fn test_resonance_semimajor_1j() {
        let a = resonance_semimajor(1, 1);
        assert!((a - A_NEPTUNE_AU).abs() < 0.01, "1:1 at a_N");
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
    fn test_octupole_narrower_than_quadrupole_at_same_a() {
        // 1:10 and 2:20 sit at the same semi-major axis (a_N·10^{2/3});
        // the octupole width must be smaller than the quadrupole width
        // there (an extra factor of α ≈ 0.215 inside the square root).
        let r_oct = &build_1j_chain(10, 10, 35.0)[0];
        let r_quad = &build_2j_chain(20, 20, 35.0)[0];
        assert!(
            (r_oct.a_nominal - r_quad.a_nominal).abs() < 0.01,
            "chains misaligned: {} vs {}",
            r_oct.a_nominal,
            r_quad.a_nominal
        );
        assert!(r_quad.delta_a > 0.0 && r_oct.delta_a > 0.0);
        // Numerically the ratio is ≈ 1.48 at q = 35 AU.
        assert!(
            r_quad.delta_a > 1.2 * r_oct.delta_a,
            "quadrupole {:.3} AU should exceed octupole {:.3} AU",
            r_quad.delta_a,
            r_oct.delta_a
        );
    }

    #[test]
    fn test_4j_chain_present_and_narrow() {
        // The hexadecapole 4:j chain exists and is narrower than the
        // quadrupole at the same semi-major axis (4:20 vs 2:10 at
        // a_N·5^{2/3}).
        let r_hex = &build_4j_chain(20, 20, 35.0)[0];
        let r_quad = &build_2j_chain(10, 10, 35.0)[0];
        assert!((r_hex.a_nominal - r_quad.a_nominal).abs() < 0.01);
        assert!(r_hex.delta_a > 0.0);
        assert!(
            r_quad.delta_a > r_hex.delta_a,
            "quadrupole {:.3} vs hexadecapole {:.3}",
            r_quad.delta_a,
            r_hex.delta_a
        );
    }

    #[test]
    fn test_width_decreases_with_q() {
        // Raising the perihelion weakens Neptune's kicks at every order.
        for chain in [
            build_2j_chain(20, 20, 30.0),
            build_1j_chain(10, 10, 30.0),
            build_3j_chain(30, 30, 30.0),
        ] {
            assert!(chain[0].delta_a > 0.0);
        }
        assert!(build_2j_chain(20, 20, 30.0)[0].delta_a > build_2j_chain(20, 20, 45.0)[0].delta_a);
        assert!(build_1j_chain(10, 10, 30.0)[0].delta_a > build_1j_chain(10, 10, 45.0)[0].delta_a);
    }

    #[test]
    fn test_resonance_width_positive() {
        let w = resonance_width(200.0, 2, 0.75, 0.1);
        assert!(w > 0.0);
    }

    #[test]
    fn typed_accessors_match_f64_fields() {
        use approx::assert_relative_eq;
        use uom::si::length::astronomical_unit;

        let r = &build_2j_chain(10, 10, 35.0)[0];
        assert_relative_eq!(
            r.semi_major_axis().get::<astronomical_unit>(),
            r.a_nominal,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            r.half_width().get::<astronomical_unit>(),
            r.delta_a,
            epsilon = 1e-12
        );
    }
}
