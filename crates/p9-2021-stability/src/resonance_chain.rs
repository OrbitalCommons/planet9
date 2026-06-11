//! The infinite chain of 2:j resonances with Neptune that shapes the scattered disk.
//!
//! Quadrupole disturbing function (Eq. in the paper):
//! R_q = (Gm_N/4a) * alpha^2 * {X_j^(-3,0) cos(jM) + 3 X_j^(-3,2) cos[jM - 2(M_N - omega)]}
//!
//! Resonance half-width from the pendulum model:
//!   Delta_a = 4 a_N sqrt(2 chi m_N / (5 M_sun)) * exp(-q^2 / (2 a_N^2))
//!
//! The exponent carries q²/(2a_N²): the pendulum width goes as the square
//! root of the resonant potential, which itself carries the Hansen factor
//! exp(−(q/a_N)²). (A previous copy used q²/4a_N², inconsistent with both the
//! potential and the overlap parameter in `p9_core::analysis::resonance`.)

use p9_core::analysis::hansen::hansen_x_neg3_2_neptune_chain;
use p9_core::constants::{A_NEPTUNE_AU, GM_SUN, MASS_NEPTUNE_SOLAR};
use serde::{Deserialize, Serialize};

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
    A_NEPTUNE_AU * (j as f64 / 2.0).powf(2.0 / 3.0)
}

/// Resonance half-width (AU) from the pendulum model.
///
/// Delta_a = 4 a_N sqrt(2 chi m_N / (5 M_sun)) * exp(-q^2/(2 a_N^2))
/// where chi = j/2 for a 2:j resonance.
pub fn resonance_width(j: i64, q: f64) -> f64 {
    let chi = j as f64 / 2.0;
    4.0 * A_NEPTUNE_AU
        * (2.0 * chi * MASS_NEPTUNE_SOLAR / 5.0).sqrt()
        * (-q * q / (2.0 * A_NEPTUNE_AU * A_NEPTUNE_AU)).exp()
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
            let hansen_coeff = hansen_x_neg3_2_neptune_chain(j, q);
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
    let kinetic = -3.0 / (A_NEPTUNE_AU * A_NEPTUNE_AU)
        * (chi / 2.0).powf(2.0 / 3.0)
        * (phi_tilde * phi_tilde / 2.0);
    let potential = -(6.0 * GM_SUN * MASS_NEPTUNE_SOLAR / (5.0 * A_NEPTUNE_AU * chi))
        * (-(q / A_NEPTUNE_AU).powi(2)).exp()
        * (2.0 * phi).cos();
    kinetic + potential
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::analysis::resonance::chirikov_overlap_parameter;

    #[test]
    fn test_resonance_semimajor_j2() {
        // 2:2 = 1:1 resonance should be at a_N
        let a = resonance_semimajor(2);
        assert!(
            (a - A_NEPTUNE_AU).abs() < 0.01,
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

    /// H1 consistency: the chain's width/spacing ratio carries exactly the
    /// same (a, q) dependence — α^{5/4} exp(−q²/2a_N²) — as the canonical
    /// overlap parameter K = (24/√5) α^{5/4} √(m_N/M☉) exp(−q²/2a_N²), up to
    /// a constant factor 1/√2 (single full width over spacing vs the
    /// sum-of-adjacent-half-widths convention in K). The previous /4a_N²
    /// exponent broke this by exp(q²/4a_N²) ≈ 1.4 at q = 35 AU.
    #[test]
    fn test_width_over_spacing_tracks_overlap_parameter() {
        let expected = 1.0 / 2.0_f64.sqrt();
        for &j in &[10_i64, 20, 40, 80] {
            for &q in &[30.0, 35.0, 45.0] {
                let a = resonance_semimajor(j);
                let chain_ratio = resonance_width(j, q) / resonance_spacing(j);
                let k = chirikov_overlap_parameter(a, q);
                let factor = chain_ratio / k;
                assert!(
                    (factor - expected).abs() / expected < 0.05,
                    "j = {j}, q = {q}: width/spacing = {chain_ratio:.4}, K = {k:.4}, factor = {factor:.4}"
                );
            }
        }
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
    fn test_pendulum_hamiltonian_potential_well() {
        // The potential term must create a well: H(phi=0) < H(phi=pi/2) at
        // fixed momentum (cos(2*phi) flips sign).
        let h_center = pendulum_hamiltonian(0.0, 0.1, 20, 35.0);
        let h_separatrix = pendulum_hamiltonian(std::f64::consts::FRAC_PI_2, 0.1, 20, 35.0);
        assert!(h_center < h_separatrix);
    }

    #[test]
    fn test_resonance_spacing_positive() {
        let ds = resonance_spacing(20);
        assert!(ds > 0.0, "spacing should be positive");
    }
}
