//! Multipole expansion of Neptune's disturbing function.
//!
//! R = (Gm_N/a) sum_l sum_m sum_{j',j} c^2_{lm} alpha^l
//!     X^{l,m}_{j'}(e_N) X^{-(l+1),m}_j(e) cos(...)
//!
//! The azimuthal coefficients c^2_{lm} are the coefficients of the
//! coplanar decomposition P_l(cos θ) = Σ_m c²_{lm} cos(mθ):
//!
//!   l=2: P₂ = 1/4 + (3/4)cos2θ          → m=0: 1/4,  m=2: 3/4
//!   l=3: P₃ = (3/8)cosθ + (5/8)cos3θ     → m=1: 3/8,  m=3: 5/8
//!   l=4: P₄ = 9/64 + (20/64)cos2θ + (35/64)cos4θ
//!                                        → m=0: 9/64, m=2: 5/16, m=4: 35/64
//!
//! Quadrupole (l=2): m=0,2 terms yield the 2:j resonance chain.
//! Octupole (l=3): m=1,3 terms yield 1:j and 3:j chains.
//! Hexadecapole (l=4): the m=4 term yields the 4:j chain.

use p9_core::analysis::hansen::hansen_coefficient;
use serde::{Deserialize, Serialize};

/// Multipole coefficients c^2_{lm}.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultipoleCoefficients {
    /// Multipole order l
    pub l: u32,
    /// Azimuthal order m
    pub m: u32,
    /// Coefficient c^2_{lm} of cos(mθ) in the P_l decomposition
    pub coefficient: f64,
}

/// Quadrupole (l=2) coefficients.
pub fn quadrupole_coefficients() -> Vec<MultipoleCoefficients> {
    vec![
        MultipoleCoefficients {
            l: 2,
            m: 0,
            coefficient: 0.25, // (1/4) for l=2, m=0
        },
        MultipoleCoefficients {
            l: 2,
            m: 2,
            coefficient: 0.75, // (3/4) for l=2, m=2
        },
    ]
}

/// Octupole (l=3) coefficients.
pub fn octupole_coefficients() -> Vec<MultipoleCoefficients> {
    vec![
        MultipoleCoefficients {
            l: 3,
            m: 1,
            coefficient: 3.0 / 8.0, // (3/8) for l=3, m=1
        },
        MultipoleCoefficients {
            l: 3,
            m: 3,
            coefficient: 5.0 / 8.0, // (5/8) for l=3, m=3
        },
    ]
}

/// Hexadecapole (l=4) coefficients, from
/// P₄(cos θ) = 9/64 + (20/64)cos 2θ + (35/64)cos 4θ.
pub fn hexadecapole_coefficients() -> Vec<MultipoleCoefficients> {
    vec![
        MultipoleCoefficients {
            l: 4,
            m: 0,
            coefficient: 9.0 / 64.0,
        },
        MultipoleCoefficients {
            l: 4,
            m: 2,
            coefficient: 20.0 / 64.0,
        },
        MultipoleCoefficients {
            l: 4,
            m: 4,
            coefficient: 35.0 / 64.0,
        },
    ]
}

/// Octupole disturbing function contribution.
///
/// R^oct = (Gm_N alpha^3 / 8a) * sum_j [3 X^{-4,1}_j cos(jM - (M_N - omega))
///         + 5 X^{-4,3}_j cos(jM - 3(M_N - omega))]
///
/// `j_values[i]` is the resonance index of `hansen_neg4_1[i]` /
/// `hansen_neg4_3[i]` (the previous version derived j from `enumerate()`,
/// silently wrong for tables not starting at j = 0).
pub fn octupole_disturbing_function(
    alpha: f64,
    gm_n: f64,
    a: f64,
    j_values: &[i64],
    hansen_neg4_1: &[f64],
    hansen_neg4_3: &[f64],
    m_anom: f64,
    m_neptune: f64,
    omega: f64,
) -> f64 {
    assert_eq!(j_values.len(), hansen_neg4_1.len());
    assert_eq!(j_values.len(), hansen_neg4_3.len());

    let prefactor = gm_n * alpha.powi(3) / (8.0 * a);
    let mut sum = 0.0;

    for ((&j, &x1), &x3) in j_values
        .iter()
        .zip(hansen_neg4_1.iter())
        .zip(hansen_neg4_3.iter())
    {
        let jf = j as f64;
        let arg1 = jf * m_anom - (m_neptune - omega);
        let arg3 = jf * m_anom - 3.0 * (m_neptune - omega);
        sum += 3.0 * x1 * arg1.cos() + 5.0 * x3 * arg3.cos();
    }

    prefactor * sum
}

/// Alpha ratio a_N / a for a test particle.
pub fn alpha_ratio(a_neptune: f64, a: f64) -> f64 {
    a_neptune / a
}

/// Secular amplitude envelope of the disturbing function up to multipole
/// order `max_order`, with proper Hansen eccentricity factors:
///
///   R(a, e) = (Gm_N/a) Σ_l α^l Σ_m c²_{lm} |X_0^{-(l+1),m}(e)|
///
/// Neptune is on a near-circular orbit, so its Hansen factor
/// X_{j'}^{l,m}(e_N → 0) = δ_{j'm} contributes exactly 1 at j' = m; the
/// particle's factors X_0^{-(l+1),m}(e) are evaluated with the p9-core
/// convergence-controlled quadrature (the previous version omitted all
/// eccentricity dependence).
pub fn total_disturbing_function(
    a: f64,
    e: f64,
    a_neptune: f64,
    gm_neptune: f64,
    max_order: u32,
) -> f64 {
    let alpha = alpha_ratio(a_neptune, a);
    let mut groups = vec![quadrupole_coefficients()];
    if max_order >= 3 {
        groups.push(octupole_coefficients());
    }
    if max_order >= 4 {
        groups.push(hexadecapole_coefficients());
    }

    let mut total = 0.0;
    for group in &groups {
        for c in group {
            let x = hansen_coefficient(0, -(c.l as i32 + 1), c.m as i32, e, 1e-8);
            total += gm_neptune / a * c.coefficient * alpha.powi(c.l as i32) * x.abs();
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quadrupole_coefficients_count() {
        let coeffs = quadrupole_coefficients();
        assert_eq!(coeffs.len(), 2);
    }

    #[test]
    fn test_octupole_coefficients_count() {
        let coeffs = octupole_coefficients();
        assert_eq!(coeffs.len(), 2);
    }

    #[test]
    fn test_multipole_coefficients_sum_to_one() {
        // P_l(1) = 1, so each azimuthal decomposition sums to 1.
        for group in [
            quadrupole_coefficients(),
            octupole_coefficients(),
            hexadecapole_coefficients(),
        ] {
            let sum: f64 = group.iter().map(|c| c.coefficient).sum();
            assert!((sum - 1.0).abs() < 1e-10, "sum = {}", sum);
        }
    }

    #[test]
    fn test_alpha_ratio_less_than_one() {
        let alpha = alpha_ratio(30.0, 200.0);
        assert!(alpha > 0.0 && alpha < 1.0, "alpha = {}", alpha);
    }

    #[test]
    fn test_total_disturbing_function_positive() {
        let gm_n = p9_core::constants::GM_NEPTUNE;
        let r = total_disturbing_function(200.0, 0.8, 30.07, gm_n, 3);
        assert!(r > 0.0, "R should be positive");
    }

    #[test]
    fn test_octupole_smaller_than_quadrupole() {
        let gm_n = p9_core::constants::GM_NEPTUNE;
        let r2 = total_disturbing_function(200.0, 0.8, 30.07, gm_n, 2);
        let r3 = total_disturbing_function(200.0, 0.8, 30.07, gm_n, 3);
        let r4 = total_disturbing_function(200.0, 0.8, 30.07, gm_n, 4);
        assert!(r3 > r2, "octupole adds positive contribution");
        assert!((r3 - r2) < r2, "octupole should be smaller than quadrupole");
        assert!(r4 > r3 && (r4 - r3) < (r3 - r2), "multipoles decay with l");
    }

    #[test]
    fn test_eccentricity_strengthens_secular_terms() {
        // X_0^{-3,0}(e) = (1−e²)^{-3/2} and the m≠0 secular factors all
        // grow with e: an eccentric orbit feels a stronger envelope.
        let gm_n = p9_core::constants::GM_NEPTUNE;
        let r_circ = total_disturbing_function(200.0, 0.0, 30.07, gm_n, 3);
        let r_ecc = total_disturbing_function(200.0, 0.8, 30.07, gm_n, 3);
        assert!(r_ecc > r_circ, "{r_ecc} should exceed {r_circ}");
        // At e = 0 the m≠0 secular Hansen factors vanish (X_0^{n,m}(0) =
        // δ_{0m}): only the m=0 quadrupole survives.
        let quad_m0 = gm_n / 200.0 * 0.25 * (30.07_f64 / 200.0).powi(2);
        assert!((r_circ - quad_m0).abs() < 1e-3 * quad_m0);
    }

    #[test]
    fn test_octupole_disturbing_function_j_indexing() {
        // A single Hansen entry at j = 7 must produce a j-dependent phase:
        // the function evaluated at m_anom = π/7 with j = 7 sees
        // cos(π − ...), not cos(0 − ...) as the old enumerate() bug did.
        let j_values = [7_i64];
        let x1 = [1.0];
        let x3 = [0.0];
        let r_aligned =
            octupole_disturbing_function(0.15, 1e-8, 200.0, &j_values, &x1, &x3, 0.0, 0.0, 0.0);
        let r_anti = octupole_disturbing_function(
            0.15,
            1e-8,
            200.0,
            &j_values,
            &x1,
            &x3,
            std::f64::consts::PI / 7.0,
            0.0,
            0.0,
        );
        // cos(0) = 1 vs cos(π) = −1 → strict sign flip.
        assert!(r_aligned > 0.0);
        assert!(r_anti < 0.0);
        assert!((r_aligned + r_anti).abs() < 1e-20);
    }
}
