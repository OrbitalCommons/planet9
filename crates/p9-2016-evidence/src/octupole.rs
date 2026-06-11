//! Secular Hamiltonian beyond quadrupole order.
//!
//! Two fidelity levels, dispatched on the semi-major axis ratio α = a/a_p:
//!
//! 1. **Analytic hierarchical expansion** (`coplanar_hierarchical_hamiltonian`):
//!    the doubly-averaged quadrupole + octupole disturbing function for an
//!    interior test particle and an exterior eccentric perturber,
//!
//!      H = −(gm_p/a_p)·[ (α²/4)(1 + 3e²/2)/(1−e_p²)^{3/2}
//!          − (15/16) α³ e e_p (1 + 3e²/4) cos Δϖ / (1−e_p²)^{5/2} ]
//!
//!    (standard coplanar octupole form, e.g. Lee & Peale 2003 Eq. 5;
//!    Naoz 2016 test-particle limit). Only the cos Δϖ harmonic appears at
//!    octupole order — a previous implementation carried unsourced
//!    coefficients and a cos 3Δϖ term that match no standard form, plus a
//!    bolted-on cos(i) "3D" factor; both are removed. Valid only in the
//!    hierarchical regime (see `hierarchical_expansion_valid`).
//!
//! 2. **Numerical Gauss-ring double averaging** via
//!    `p9_core::analysis::secular::numerical_secular_hamiltonian`: the exact
//!    doubly-averaged 1/Δ. The multipole expansion diverges for the
//!    α ≈ 0.3–0.8 of every use in this crate (test particles cross P9's
//!    orbit); Batygin & Brown (2016) used numerical ring averaging for
//!    exactly this reason, and the phase portraits here do the same.
//!
//! Reference: Batygin & Brown (2016) Section 3; Lee & Peale (2003);
//! Naoz (2016)

use p9_core::analysis::secular;

/// Octupole coupling strength epsilon_oct = (a/a_p) * e_p / (1 - e_p^2).
/// When this is large, the octupole term significantly modifies the dynamics.
pub fn octupole_epsilon(a: f64, a_p: f64, e_p: f64) -> f64 {
    (a / a_p) * e_p / (1.0 - e_p * e_p)
}

/// Largest α = a/a_p for which the truncated quadrupole+octupole expansion
/// is used. Beyond this the multipole series converges too slowly (or not at
/// all, once orbits cross) and the numerical average is required.
pub const HIERARCHICAL_ALPHA_MAX: f64 = 0.1;

/// Whether the analytic hierarchical expansion is trustworthy: requires both
/// α = a/a_p below `HIERARCHICAL_ALPHA_MAX` and a non-crossing geometry
/// (particle aphelion inside the perturber perihelion).
pub fn hierarchical_expansion_valid(a: f64, e: f64, a_p: f64, e_p: f64) -> bool {
    (a / a_p) < HIERARCHICAL_ALPHA_MAX && a * (1.0 + e) < a_p * (1.0 - e_p)
}

/// Coplanar secular Hamiltonian to octupole order (analytic, hierarchical
/// regime only — see `hierarchical_expansion_valid` and the module docs for
/// the source of the coefficients).
///
/// `a`: test particle semi-major axis (AU)
/// `e`: test particle eccentricity
/// `delta_varpi`: longitude of perihelion difference (radians)
/// `a_p`: perturber semi-major axis (AU)
/// `e_p`: perturber eccentricity
/// `gm_p`: perturber GM (AU^3/day^2)
pub fn coplanar_hierarchical_hamiltonian(
    a: f64,
    e: f64,
    delta_varpi: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
) -> f64 {
    let alpha = a / a_p;
    let e2 = e * e;

    // Quadrupole part (axisymmetric in Δϖ at this order)
    let h_quad = secular::coplanar_quadrupole(a, e, a_p, e_p, gm_p);

    // Octupole: H_oct = +(15/16)(gm_p/a_p) α³ e e_p (1 + 3e²/4) cosΔϖ / (1−e_p²)^{5/2}
    let h_oct = (15.0 / 16.0) * (gm_p / a_p) * alpha.powi(3) * e * e_p * (1.0 + 0.75 * e2)
        / (1.0 - e_p * e_p).powf(2.5)
        * delta_varpi.cos();

    h_quad + h_oct
}

/// Quadrature nodes per anomaly for the numerical double average.
const N_QUAD: usize = 64;

/// Softening for crossing geometries, as a fraction of a_p (matches
/// p9-core's portrait default).
const SOFTENING_FRAC: f64 = 0.01;

/// Coplanar secular Hamiltonian, automatically choosing the analytic
/// hierarchical expansion (small, non-crossing α) or the numerically
/// double-averaged exact interaction (everything else — the regime of every
/// phase portrait in this crate).
pub fn coplanar_secular_hamiltonian(
    a: f64,
    e: f64,
    delta_varpi: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
) -> f64 {
    if hierarchical_expansion_valid(a, e, a_p, e_p) {
        coplanar_hierarchical_hamiltonian(a, e, delta_varpi, a_p, e_p, gm_p)
    } else {
        secular::numerical_secular_hamiltonian(
            a,
            e,
            0.0,
            delta_varpi,
            0.0,
            a_p,
            e_p,
            gm_p,
            N_QUAD,
            SOFTENING_FRAC * a_p,
        )
    }
}

/// Generate a coplanar phase-space portrait of the secular Hamiltonian using
/// the numerical Gauss-ring double average (the α here is far outside the
/// hierarchical regime).
///
/// Returns (e_vals, dvarpi_vals, portrait) where portrait[i][j] = H(e_i, dvarpi_j).
pub fn octupole_phase_portrait(
    a: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
    n_e: usize,
    n_dvarpi: usize,
) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    secular::phase_portrait(a, a_p, e_p, gm_p, n_e, n_dvarpi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::types::P9Params;
    use std::f64::consts::PI;

    #[test]
    fn test_octupole_reduces_to_quadrupole_at_small_alpha() {
        // When a << a_p, the octupole correction is O(α·e_p) of the
        // quadrupole and the analytic H approaches the quadrupole-only H.
        let p9 = P9Params::nominal_2016();
        let gm_p = p9.gm();

        let a = 50.0; // alpha = 50/700 = 0.071
        let e = 0.5;
        let dv = 1.0;

        let h_oct = coplanar_hierarchical_hamiltonian(a, e, dv, p9.a, p9.e, gm_p);
        let h_quad = secular::coplanar_quadrupole(a, e, p9.a, p9.e, gm_p);

        let rel_diff = ((h_oct - h_quad) / h_quad).abs();
        assert!(
            rel_diff < 0.1,
            "Octupole differs from quadrupole by {:.1}% at small alpha",
            rel_diff * 100.0
        );
    }

    #[test]
    fn test_analytic_octupole_matches_numerical_at_small_alpha() {
        // The Δϖ-dependent part of the exact double average is dominated by
        // the octupole at small α. The aligned/anti-aligned difference
        //   ΔH = H(Δϖ=0) − H(Δϖ=π) = 2·(15/16)(gm_p/a_p) α³ e e_p (1+3e²/4)/(1−e_p²)^{5/2}
        // must match the numerical average (which carries all higher
        // multipoles; the next odd term is O(α²) ≈ 0.6% here).
        let (a, e, a_p, e_p) = (55.0, 0.4, 700.0, 0.6);
        let gm_p = P9Params::nominal_2016().gm();
        assert!(hierarchical_expansion_valid(a, e, a_p, e_p));

        let d_analytic = coplanar_hierarchical_hamiltonian(a, e, 0.0, a_p, e_p, gm_p)
            - coplanar_hierarchical_hamiltonian(a, e, PI, a_p, e_p, gm_p);
        let num = |dv: f64| {
            secular::numerical_secular_hamiltonian(a, e, 0.0, dv, 0.0, a_p, e_p, gm_p, 128, 0.0)
        };
        let d_numeric = num(0.0) - num(PI);

        let rel = ((d_analytic - d_numeric) / d_numeric).abs();
        assert!(
            rel < 0.05,
            "analytic octupole ΔH = {d_analytic:.3e} vs numerical {d_numeric:.3e} (rel {rel:.2e})"
        );
    }

    #[test]
    fn test_octupole_breaks_aligned_antialigned_symmetry() {
        // At quadrupole order H(e, 0) = H(e, π). The octupole breaks this by
        // an amount with a definite physical scale: |ΔH| of order
        // ε_oct·e·|H_quad| (not merely nonzero at machine precision).
        let p9 = P9Params::nominal_2016();
        let gm_p = p9.gm();

        let a = 350.0; // alpha = 0.5: numerical-average regime
        let e = 0.7;

        let h_aligned = coplanar_secular_hamiltonian(a, e, 0.0, p9.a, p9.e, gm_p);
        let h_anti = coplanar_secular_hamiltonian(a, e, PI, p9.a, p9.e, gm_p);
        let h_quad = secular::coplanar_quadrupole(a, e, p9.a, p9.e, gm_p);

        let scale = octupole_epsilon(a, p9.a, p9.e) * e * h_quad.abs();
        let dh = (h_aligned - h_anti).abs();
        assert!(
            dh > 0.01 * scale,
            "Δϖ asymmetry |ΔH| = {dh:.3e} below the octupole scale {scale:.3e}"
        );
    }

    #[test]
    fn test_dispatch_uses_numerical_outside_hierarchical_regime() {
        // For crossing geometry the dispatcher must agree with the numerical
        // average, not the (divergent) expansion.
        let p9 = P9Params::nominal_2016();
        let gm_p = p9.gm();
        let (a, e) = (350.0, 0.7);
        assert!(!hierarchical_expansion_valid(a, e, p9.a, p9.e));

        let h = coplanar_secular_hamiltonian(a, e, 1.0, p9.a, p9.e, gm_p);
        let h_num = secular::numerical_secular_hamiltonian(
            a,
            e,
            0.0,
            1.0,
            0.0,
            p9.a,
            p9.e,
            gm_p,
            N_QUAD,
            SOFTENING_FRAC * p9.a,
        );
        assert_eq!(h, h_num);
    }

    #[test]
    fn test_phase_portrait_at_six_semimajor_axes() {
        // Paper uses a = 50, 150, 250, 350, 450, 550 AU
        let p9 = P9Params::nominal_2016();
        let gm_p = p9.gm();

        let test_axes = [50.0, 150.0, 250.0, 350.0, 450.0, 550.0];

        for &a in &test_axes {
            let (e_vals, dv_vals, portrait) = octupole_phase_portrait(a, p9.a, p9.e, gm_p, 8, 12);

            assert_eq!(e_vals.len(), 8);
            assert_eq!(dv_vals.len(), 12);
            assert_eq!(portrait.len(), 8);
            assert_eq!(portrait[0].len(), 12);

            // Hamiltonian should be finite and negative (bound interaction)
            for row in &portrait {
                for &h in row {
                    assert!(h.is_finite() && h < 0.0, "H invalid at a = {} AU", a);
                }
            }
        }
    }
}
