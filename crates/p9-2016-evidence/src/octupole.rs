//! Octupole-order secular Hamiltonian.
//!
//! Extends p9-core's quadrupole Hamiltonian (Eq. 1 of Batygin & Brown 2016)
//! with the octupole correction (Eq. 5), which becomes important when the
//! semi-major axis ratio alpha = a/a_p is not small.
//!
//! The octupole Hamiltonian captures asymmetry between aligned and anti-aligned
//! libration islands, breaking the Δϖ → -Δϖ symmetry present at quadrupole order.
//!
//! Reference: Batygin & Brown (2016), Equations 1-5; Mardling (2013)

use std::f64::consts::PI;

/// Octupole coupling coefficient epsilon_oct = (a/a_p) * e_p / (1 - e_p^2).
/// When this is large, the octupole term significantly modifies the dynamics.
pub fn octupole_epsilon(a: f64, a_p: f64, e_p: f64) -> f64 {
    (a / a_p) * e_p / (1.0 - e_p * e_p)
}

/// Coplanar secular Hamiltonian to octupole order.
///
/// H = H_quad + H_oct
///
/// H_quad = -C_quad * [2 + 3e^2 - 5e^2 cos^2(Δϖ)]
///
/// H_oct = C_oct * e * [A * cos(Δϖ) + B * cos(3Δϖ)]
///
/// where the octupole terms (from Mardling 2013, adapted by Batygin):
///   A = -35/8 * e^2 + 25/8
///   B = -35/8 * e^2 + 15/8
///
/// `a`: test particle semi-major axis (AU)
/// `e`: test particle eccentricity
/// `delta_varpi`: longitude of perihelion difference (radians)
/// `a_p`: perturber semi-major axis (AU)
/// `e_p`: perturber eccentricity
/// `gm_p`: perturber GM (AU^3/day^2)
pub fn coplanar_octupole_hamiltonian(
    a: f64,
    e: f64,
    delta_varpi: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
) -> f64 {
    let alpha = a / a_p;
    let e2 = e * e;

    // Quadrupole coupling constant
    let c_quad = gm_p * alpha * alpha / (4.0 * a_p * (1.0 - e_p * e_p).powf(1.5));

    // Quadrupole Hamiltonian (coplanar)
    let h_quad = -c_quad * (2.0 + 3.0 * e2 - 5.0 * e2 * delta_varpi.cos().powi(2));

    // Octupole coupling constant
    let eps_oct = octupole_epsilon(a, a_p, e_p);
    let c_oct = c_quad * eps_oct * 15.0 / 4.0;

    // Octupole terms
    let cos_dv = delta_varpi.cos();
    let cos_3dv = (3.0 * delta_varpi).cos();

    let a_coeff = -35.0 / 8.0 * e2 + 25.0 / 8.0;
    let b_coeff = -35.0 / 8.0 * e2 + 15.0 / 8.0;

    let h_oct = c_oct * e * (a_coeff * cos_dv + b_coeff * cos_3dv);

    h_quad + h_oct
}

/// Full 3D secular Hamiltonian to octupole order.
///
/// Includes inclination-dependent terms from the Mardling (2013) expansion.
/// For the inclined case, the quadrupole part uses the full p9-core formula,
/// and the octupole adds inclination-modulated corrections.
pub fn full_octupole_hamiltonian(
    a: f64,
    e: f64,
    delta_varpi: f64,
    i: f64,
    delta_omega: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
) -> f64 {
    // Quadrupole part from p9-core
    let h_quad = p9_core::analysis::secular::quadrupole_hamiltonian(
        a,
        e,
        delta_varpi,
        i,
        delta_omega,
        a_p,
        e_p,
        gm_p,
    );

    let alpha = a / a_p;
    let e2 = e * e;

    // Octupole coupling
    let c_quad = gm_p * alpha * alpha / (4.0 * a_p * (1.0 - e_p * e_p).powf(1.5));
    let eps_oct = octupole_epsilon(a, a_p, e_p);
    let c_oct = c_quad * eps_oct * 15.0 / 4.0;

    // Coplanar octupole terms (simplified — the full 3D octupole has
    // many inclination terms; for now we use the dominant coplanar part)
    // TODO: Implement full 3D octupole from Mardling (2013) Eq. A11-A14
    let cos_i = i.cos();

    let cos_dv = delta_varpi.cos();
    let cos_3dv = (3.0 * delta_varpi).cos();

    let a_coeff = -35.0 / 8.0 * e2 + 25.0 / 8.0;
    let b_coeff = -35.0 / 8.0 * e2 + 15.0 / 8.0;

    // Inclination modulation: in the coplanar limit (i=0), cos_i = 1
    let h_oct = c_oct * e * cos_i * (a_coeff * cos_dv + b_coeff * cos_3dv);

    h_quad + h_oct
}

/// Generate phase-space portrait with octupole Hamiltonian.
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
    let e_vals: Vec<f64> = (0..n_e)
        .map(|i| (i as f64 + 0.5) / n_e as f64 * 0.95)
        .collect();

    let dvarpi_vals: Vec<f64> = (0..n_dvarpi)
        .map(|j| (j as f64) / n_dvarpi as f64 * 2.0 * PI - PI)
        .collect();

    let mut portrait = vec![vec![0.0; n_dvarpi]; n_e];

    for (i, &e) in e_vals.iter().enumerate() {
        for (j, &dv) in dvarpi_vals.iter().enumerate() {
            portrait[i][j] = coplanar_octupole_hamiltonian(a, e, dv, a_p, e_p, gm_p);
        }
    }

    (e_vals, dvarpi_vals, portrait)
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::types::P9Params;

    #[test]
    fn test_octupole_reduces_to_quadrupole_at_small_alpha() {
        // When a << a_p, octupole epsilon -> 0, so octupole should match quadrupole
        let p9 = P9Params::nominal_2016();
        let gm_p = p9.gm();

        let a = 50.0; // Very small alpha = 50/700 = 0.071
        let e = 0.5;
        let dv = 1.0;

        let h_oct = coplanar_octupole_hamiltonian(a, e, dv, p9.a, p9.e, gm_p);
        let h_quad = p9_core::analysis::secular::coplanar_quadrupole(a, e, dv, p9.a, p9.e, gm_p);

        // Should be very close since alpha is small
        let rel_diff = ((h_oct - h_quad) / h_quad).abs();
        assert!(
            rel_diff < 0.1,
            "Octupole differs from quadrupole by {:.1}% at small alpha",
            rel_diff * 100.0
        );
    }

    #[test]
    fn test_octupole_breaks_aligned_antialigned_symmetry() {
        // At quadrupole order, H(e, 0) = H(e, π) because cos²(0) = cos²(π) = 1.
        // The octupole breaks this: cos(0) = 1 ≠ cos(π) = -1.
        let p9 = P9Params::nominal_2016();
        let gm_p = p9.gm();

        let a = 350.0; // Larger alpha where octupole matters
        let e = 0.7;

        let h_aligned = coplanar_octupole_hamiltonian(a, e, 0.0, p9.a, p9.e, gm_p);
        let h_anti = coplanar_octupole_hamiltonian(a, e, PI, p9.a, p9.e, gm_p);

        // These should NOT be equal (unlike pure quadrupole where cos²(0) = cos²(π))
        assert!(
            (h_aligned - h_anti).abs() > 1e-20,
            "Octupole should break aligned/anti-aligned symmetry: H(0) = {}, H(π) = {}",
            h_aligned,
            h_anti
        );
    }

    #[test]
    fn test_phase_portrait_at_six_semimajor_axes() {
        // Paper uses a = 50, 150, 250, 350, 450, 550 AU
        let p9 = P9Params::nominal_2016();
        let gm_p = p9.gm();

        let test_axes = [50.0, 150.0, 250.0, 350.0, 450.0, 550.0];

        for &a in &test_axes {
            let (e_vals, dv_vals, portrait) = octupole_phase_portrait(a, p9.a, p9.e, gm_p, 20, 40);

            assert_eq!(e_vals.len(), 20);
            assert_eq!(dv_vals.len(), 40);
            assert_eq!(portrait.len(), 20);
            assert_eq!(portrait[0].len(), 40);

            // Hamiltonian should be finite everywhere
            for row in &portrait {
                for &h in row {
                    assert!(h.is_finite(), "H not finite at a = {} AU", a);
                }
            }
        }
    }
}
