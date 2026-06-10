//! Multipole expansion of Neptune's disturbing function.
//!
//! R = (Gm_N/a) sum_l sum_m sum_{j',j} c^2_{lm} M_l alpha^l
//!     X^{l,m}_{j'}(e_N) X^{-(l+1),m}_j(e) cos(...)
//!
//! Quadrupole (l=2): m=0,2 terms yield the 2:j resonance chain.
//! Octupole (l=3): m=1,3 terms yield 1:j and 3:j chains.

use serde::{Deserialize, Serialize};

/// Multipole coefficients c^2_{lm} and M_l.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultipoleCoefficients {
    /// Multipole order l
    pub l: u32,
    /// Azimuthal order m
    pub m: u32,
    /// Combined coefficient c^2_{lm} * M_l
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

/// Octupole disturbing function contribution.
///
/// R^oct = (Gm_N alpha^3 / 8a) * sum[3 X^{-4,1}_j cos(jM - (M_N - omega))
///         + 5 X^{-4,3}_j cos(jM - 3(M_N - omega))]
pub fn octupole_disturbing_function(
    alpha: f64,
    gm_n: f64,
    a: f64,
    hansen_neg4_1: &[f64],
    hansen_neg4_3: &[f64],
    m_anom: f64,
    m_neptune: f64,
    omega: f64,
) -> f64 {
    let prefactor = gm_n * alpha.powi(3) / (8.0 * a);
    let mut sum = 0.0;

    for (j, (&x1, &x3)) in hansen_neg4_1.iter().zip(hansen_neg4_3.iter()).enumerate() {
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

/// Total disturbing function up to a given multipole order.
pub fn total_disturbing_function(a: f64, a_neptune: f64, gm_neptune: f64, max_order: u32) -> f64 {
    let alpha = alpha_ratio(a_neptune, a);
    let mut total = 0.0;

    // Quadrupole contribution (always included)
    let quad = quadrupole_coefficients();
    for c in &quad {
        total += gm_neptune / a * c.coefficient * alpha.powi(c.l as i32);
    }

    // Octupole contribution (if max_order >= 3)
    if max_order >= 3 {
        let oct = octupole_coefficients();
        for c in &oct {
            total += gm_neptune / a * c.coefficient * alpha.powi(c.l as i32);
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
    fn test_quadrupole_sum_to_one() {
        let coeffs = quadrupole_coefficients();
        let sum: f64 = coeffs.iter().map(|c| c.coefficient).sum();
        assert!((sum - 1.0).abs() < 1e-10, "sum = {}", sum);
    }

    #[test]
    fn test_octupole_sum_to_one() {
        let coeffs = octupole_coefficients();
        let sum: f64 = coeffs.iter().map(|c| c.coefficient).sum();
        assert!((sum - 1.0).abs() < 1e-10, "sum = {}", sum);
    }

    #[test]
    fn test_alpha_ratio_less_than_one() {
        let alpha = alpha_ratio(30.0, 200.0);
        assert!(alpha > 0.0 && alpha < 1.0, "alpha = {}", alpha);
    }

    #[test]
    fn test_total_disturbing_function_positive() {
        let gm_n = 6.836e-12; // GM_Neptune in AU^3/day^2
        let r = total_disturbing_function(200.0, 30.07, gm_n, 3);
        assert!(r > 0.0, "R should be positive");
    }

    #[test]
    fn test_octupole_smaller_than_quadrupole() {
        let gm_n = 6.836e-12;
        let r2 = total_disturbing_function(200.0, 30.07, gm_n, 2);
        let r3 = total_disturbing_function(200.0, 30.07, gm_n, 3);
        assert!(r3 > r2, "octupole adds positive contribution");
        assert!((r3 - r2) < r2, "octupole should be smaller than quadrupole");
    }
}
