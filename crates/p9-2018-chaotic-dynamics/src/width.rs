//! Secular / resonant libration width and its scaling with Planet Nine mass.
//!
//! The pendulum half-width of a mean-motion resonance scales as
//! δa ∝ √(m₉) (the resonant harmonic is linear in the perturber mass, and the
//! pendulum half-width is the square root of the harmonic amplitude). The
//! overlap of these widening islands is what tips the dynamics into chaos, so
//! the libration width is the natural single-resonance precursor of the chaos
//! map in [`crate::chaos`].
//!
//! This module builds the closed-form pendulum half-width directly on
//! `p9_core::analysis::hansen::hansen_coefficient` (the same multipole form used
//! by the 2017-dynamics sibling), so the √m scaling and the eccentricity
//! growth come from real shared physics rather than a hand-tuned curve.

use p9_core::analysis::hansen::hansen_coefficient;
use p9_core::analysis::resonance::resonance_semi_major_axis;
use p9_core::constants::EARTH_MASS_SOLAR;
use p9_core::units::{au, Length};

/// Quadrupole-limit pendulum half-width (AU) of the interior j:1 resonance at
/// particle eccentricity `e`, for a Planet Nine of mass `m9_earth` (Earth
/// masses) and semi-major axis `a9_au`:
///
///   δa/a = √( (16/3) (m₉/M☉) α^{j+1} |X_1^{j,j}(e)| ),   α = a_j/a₉
///
/// with a_j = a₉ (1/j)^{2/3}. The relevant Hansen coefficient X_1^{j,j}(e)
/// (from `p9_core::analysis::hansen`) carries the O(1) eccentricity dependence:
/// it vanishes for a circular orbit and grows strongly toward e → 1, which is
/// what drives the resonance overlap. `e9` is accepted for API symmetry with
/// the chaos map but does not enter the leading multipole width.
pub fn libration_half_width_typed(j: i64, e: f64, m9_earth: f64, a9_au: f64, _e9: f64) -> Length {
    let m9_solar = m9_earth * EARTH_MASS_SOLAR;
    let a_j = resonance_semi_major_axis(1, j as u32, a9_au);
    let alpha = a_j / a9_au;
    let x = hansen_coefficient(1, j as i32, j as i32, e, 1e-8).abs();
    au(a_j * ((16.0 / 3.0) * m9_solar * alpha.powi(j as i32 + 1) * x).sqrt())
}

/// Full libration width as a typed [`Length`] — twice the half-width — of the
/// j:1 resonance (see [`libration_half_width_typed`]).
pub fn libration_width_typed(j: i64, e: f64, m9_earth: f64, a9_au: f64, e9: f64) -> Length {
    2.0 * libration_half_width_typed(j, e, m9_earth, a9_au, e9)
}

/// Libration half-width as a function of mass over a list of Planet Nine
/// masses (Earth masses), all other parameters fixed. Returns parallel
/// (mass, half-width) pairs for trend analysis.
pub fn width_vs_mass(
    masses_earth: &[f64],
    j: i64,
    e: f64,
    a9_au: f64,
    e9: f64,
) -> Vec<(f64, Length)> {
    masses_earth
        .iter()
        .map(|&m| (m, libration_half_width_typed(j, e, m, a9_au, e9)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference;
    use approx::assert_relative_eq;

    const A9: f64 = reference::FIDUCIAL_A9_AU;
    const E9: f64 = reference::FIDUCIAL_E9;

    #[test]
    fn libration_width_grows_monotonically_with_mass() {
        let masses = [1.0, 5.0, 10.0, 20.0, 50.0];
        let widths = width_vs_mass(&masses, 3, 0.7, A9, E9);
        for w in widths.windows(2) {
            let lo = (w[0].1 / au(1.0)).value;
            let hi = (w[1].1 / au(1.0)).value;
            assert!(
                hi > lo,
                "width must grow with mass: {} M⊕ -> {} AU, {} M⊕ -> {} AU",
                w[0].0,
                lo,
                w[1].0,
                hi
            );
        }
    }

    #[test]
    fn libration_half_width_follows_sqrt_mass_scaling() {
        // δa ∝ √m₉: doubling the mass multiplies the half-width by √2.
        let w1 = libration_half_width_typed(3, 0.7, 10.0, A9, E9);
        let w2 = libration_half_width_typed(3, 0.7, 40.0, A9, E9);
        // 4x mass -> 2x width.
        assert_relative_eq!((w2 / w1).value, 2.0, max_relative = 1e-6);
    }

    #[test]
    fn typed_width_matches_inline_formula() {
        use p9_core::analysis::hansen::hansen_coefficient;
        use p9_core::analysis::resonance::resonance_semi_major_axis;
        use p9_core::constants::EARTH_MASS_SOLAR;
        // Reproduce the closed-form half-width inline and compare to the typed
        // accessor (in AU).
        let (j, e, m9_earth) = (3i64, 0.7, 10.0);
        let m9_solar = m9_earth * EARTH_MASS_SOLAR;
        let a_j = resonance_semi_major_axis(1, j as u32, A9);
        let alpha = a_j / A9;
        let x = hansen_coefficient(1, j as i32, j as i32, e, 1e-8).abs();
        let expected = a_j * ((16.0 / 3.0) * m9_solar * alpha.powi(j as i32 + 1) * x).sqrt();
        let half = libration_half_width_typed(j, e, m9_earth, A9, E9);
        assert_relative_eq!((half / au(1.0)).value, expected, max_relative = 1e-12);
        // Full width is exactly twice the half-width.
        let full = libration_width_typed(j, e, m9_earth, A9, E9);
        assert_relative_eq!((full / au(1.0)).value, 2.0 * expected, max_relative = 1e-12);
    }

    #[test]
    fn full_width_is_twice_half_width() {
        let h = libration_half_width_typed(4, 0.6, 10.0, A9, E9);
        let f = libration_width_typed(4, 0.6, 10.0, A9, E9);
        assert_relative_eq!(
            (f / au(1.0)).value,
            2.0 * (h / au(1.0)).value,
            max_relative = 1e-12
        );
    }

    #[test]
    fn width_is_finite_and_positive() {
        for &j in &[2i64, 3, 5, 8] {
            for &e in &[0.1, 0.5, 0.85] {
                let w = (libration_half_width_typed(j, e, 10.0, A9, E9) / au(1.0)).value;
                assert!(w.is_finite() && w > 0.0, "j={j} e={e} w={w}");
            }
        }
    }
}
