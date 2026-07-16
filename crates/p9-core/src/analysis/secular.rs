//! Secular perturbation theory.
//!
//! Two levels of fidelity:
//!
//! 1. **Analytic quadrupole** (`quadrupole_hamiltonian`): the doubly-averaged
//!    quadrupole for an inner test particle and an outer eccentric perturber.
//!    This is *axisymmetric in the apsidal angle*: it depends on (e, i, ω)
//!    only — at quadrupole order there is no Δϖ dependence and therefore no
//!    aligned/anti-aligned structure. (The previous implementation carried an
//!    unphysical −5e²cos²Δϖ term.) Valid only for α = a/a_p ≪ 1.
//!
//! 2. **Numerical Gauss-ring double averaging**
//!    (`numerical_secular_hamiltonian`): the exact 1/Δ averaged over both
//!    mean anomalies. The multipole expansion diverges for the α ≈ 0.3–0.9
//!    relevant to the Planet Nine problem; Batygin & Brown (2016) used
//!    numerical ring averaging for exactly this reason. This captures the
//!    octupole-and-higher Δϖ dependence (the anti-aligned libration islands)
//!    automatically. For orbit-crossing geometries the doubly averaged
//!    integral is singular and a softening length must be supplied.
//!
//! Reference: Kozai (1962); Naoz (2016) §3; Batygin & Brown (2016)

use std::f64::consts::PI;

use crate::constants::{GM_SUN, G_AU3_MSUN_DAY2, YEAR_DAYS};

/// Radians per arcsecond (1 arcsec = π/648000 rad).
const ARCSEC_PER_RAD: f64 = 648_000.0 / PI;

/// Days per Julian century.
const DAYS_PER_CENTURY: f64 = 100.0 * YEAR_DAYS;

/// Secular precession rate of an inner planet's longitude of perihelion induced
/// by a coplanar external perturber, in **arcseconds per century**.
///
/// This is the doubly orbit-averaged quadrupole apsidal rate (Lagrange's
/// planetary equation for ϖ in the coplanar limit applied to
/// [`coplanar_quadrupole`]), whose closed form is
///
/// ```text
///   d(ϖ)/dt = (3/4) n (m_p/M_sun) (a/a_p)^3 * G(e, e_p),
///   G(e, e_p) = sqrt(1 - e^2) / (1 - e_p^2)^{3/2},
/// ```
///
/// with `n = sqrt(GM_sun / a^3)` the planet's Keplerian mean motion. The rate
/// scales linearly in the perturber mass and as the inverse cube of its
/// distance (the tidal-field gradient ∝ GM_p / a_p^3), and is independent of
/// the planet's own eccentricity to leading order.
///
/// * `a_au`: inner planet semi-major axis (AU)
/// * `e`: inner planet eccentricity
/// * `mass_p_solar`: perturber mass in solar masses
/// * `a_p_au`: perturber semi-major axis (AU)
/// * `e_p`: perturber eccentricity
pub fn perturber_perihelion_precession(
    a_au: f64,
    e: f64,
    mass_p_solar: f64,
    a_p_au: f64,
    e_p: f64,
) -> f64 {
    let n = (GM_SUN / a_au.powi(3)).sqrt(); // rad/day
    let alpha = a_au / a_p_au;
    let geometry = (1.0 - e * e).sqrt() / (1.0 - e_p * e_p).powf(1.5);

    // rate in rad/day
    let rate_rad_per_day = 0.75 * n * mass_p_solar * alpha * alpha * alpha * geometry;

    // rad/day -> arcsec/century
    rate_rad_per_day * DAYS_PER_CENTURY * ARCSEC_PER_RAD
}

/// Quadrupole coupling constant C = gm_p α² / (4 a_p (1 − e_p²)^{3/2}).
fn coupling(a: f64, a_p: f64, e_p: f64, gm_p: f64) -> f64 {
    let alpha = a / a_p;
    gm_p * alpha * alpha / (4.0 * a_p * (1.0 - e_p * e_p).powf(1.5))
}

/// Doubly-averaged quadrupole Hamiltonian for an inner test particle
/// perturbed by an outer body (Kozai-Lidov form):
///
///   H_quad = −(C/4) · [ (2 + 3e²)(3cos²i − 1) + 15 e² sin²i cos 2ω ]
///
/// with C = gm_p α²/(4 a_p (1−e_p²)^{3/2}). Angles are measured in the
/// perturber's orbital plane: `i` is the mutual inclination and `omega` the
/// particle's argument of perihelion from the common node.
///
/// Note the Hamiltonian is independent of the nodal and apsidal longitudes —
/// at quadrupole order the averaged perturber is an axisymmetric ring, so no
/// Δϖ-dependent (aligned/anti-aligned) dynamics exists at this order. Use
/// `numerical_secular_hamiltonian` for apsidal structure.
///
/// `a`: test particle semi-major axis (AU)
/// `e`: test particle eccentricity
/// `i`: mutual inclination (radians)
/// `omega`: argument of perihelion from the common node (radians)
/// `a_p`: perturber semi-major axis (AU)
/// `e_p`: perturber eccentricity
/// `gm_p`: perturber GM (AU^3/day^2)
pub fn quadrupole_hamiltonian(
    a: f64,
    e: f64,
    i: f64,
    omega: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
) -> f64 {
    let c = coupling(a, a_p, e_p, gm_p);
    let e2 = e * e;
    let cos_i2 = i.cos().powi(2);
    let sin_i2 = 1.0 - cos_i2;

    -(c / 4.0)
        * ((2.0 + 3.0 * e2) * (3.0 * cos_i2 - 1.0) + 15.0 * e2 * sin_i2 * (2.0 * omega).cos())
}

/// Coplanar quadrupole Hamiltonian (i = 0):
///   H = −C (2 + 3e²) / 2.
///
/// A function of e alone — the coplanar quadrupole has *no* libration islands;
/// those appear only at octupole order (∝ e e_p α³) and beyond, which the
/// numerical averaging includes.
pub fn coplanar_quadrupole(a: f64, e: f64, a_p: f64, e_p: f64, gm_p: f64) -> f64 {
    quadrupole_hamiltonian(a, e, 0.0, 0.0, a_p, e_p, gm_p)
}

/// Position on an orbit at eccentric anomaly `ea`, rotated into the reference
/// frame by (omega, i, omega_big). Returns the 3D position (units of `a`).
fn orbit_position(a: f64, e: f64, i: f64, omega: f64, omega_big: f64, ea: f64) -> [f64; 3] {
    let x_orb = a * (ea.cos() - e);
    let y_orb = a * (1.0 - e * e).sqrt() * ea.sin();

    let (sw, cw) = omega.sin_cos();
    let (si, ci) = i.sin_cos();
    let (so, co) = omega_big.sin_cos();

    let x = (co * cw - so * sw * ci) * x_orb + (-co * sw - so * cw * ci) * y_orb;
    let y = (so * cw + co * sw * ci) * x_orb + (-so * sw + co * cw * ci) * y_orb;
    let z = sw * si * x_orb + cw * si * y_orb;

    [x, y, z]
}

/// Doubly-averaged secular Hamiltonian from numerical Gauss-ring averaging of
/// the exact interaction:
///
///   H = −gm_p ⟨⟨ 1/Δ ⟩⟩
///     = −gm_p/(4π²) ∮∮ dM dM_p / sqrt(Δ² + b²)
///
/// (The indirect term −gm_p r·r_p/r_p³ double-averages to zero because the
/// orbit average of the Keplerian acceleration vanishes.)
///
/// The integral is evaluated on uniform eccentric-anomaly grids with the
/// Jacobians dM = (1 − e cos E) dE, which is spectrally accurate for smooth
/// (non-crossing) geometries. For crossing orbits the unsoftened average
/// diverges; pass a softening length `softening_au` > 0 (a fraction of a_p)
/// to regularize — the resulting portraits follow Batygin & Brown (2016).
///
/// The perturber's apsidal line defines the x-axis and its orbit plane the
/// reference plane, so the particle's `omega`/`omega_big` are relative to the
/// perturber (coplanar: Δϖ = omega + omega_big).
///
/// `n_quad` is the number of quadrature nodes per anomaly (n_quad² total
/// evaluations); 64–128 suffices for portrait work.
#[allow(clippy::too_many_arguments)]
pub fn numerical_secular_hamiltonian(
    a: f64,
    e: f64,
    i: f64,
    omega: f64,
    omega_big: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
    n_quad: usize,
    softening_au: f64,
) -> f64 {
    let b2 = softening_au * softening_au;
    let de = 2.0 * PI / n_quad as f64;

    // Precompute perturber ring samples (perturber defines the frame:
    // i_p = 0, omega_p = 0, Omega_p = 0).
    let mut p_pos = Vec::with_capacity(n_quad);
    let mut p_weight = Vec::with_capacity(n_quad);
    for kp in 0..n_quad {
        let ea_p = (kp as f64 + 0.5) * de;
        p_pos.push(orbit_position(a_p, e_p, 0.0, 0.0, 0.0, ea_p));
        p_weight.push(1.0 - e_p * ea_p.cos());
    }

    let mut sum = 0.0;
    for k in 0..n_quad {
        let ea = (k as f64 + 0.5) * de;
        let w = 1.0 - e * ea.cos();
        let r = orbit_position(a, e, i, omega, omega_big, ea);

        for kp in 0..n_quad {
            let rp = &p_pos[kp];
            let dx = r[0] - rp[0];
            let dy = r[1] - rp[1];
            let dz = r[2] - rp[2];
            let delta = (dx * dx + dy * dy + dz * dz + b2).sqrt();
            sum += w * p_weight[kp] / delta;
        }
    }

    // (1/2π)² ∮∮ ... dE dE_p with the dM Jacobians folded into the weights.
    -gm_p * sum * (de * de) / (4.0 * PI * PI)
}

/// Generate a phase-space portrait grid: evaluate the numerically averaged
/// Hamiltonian at a grid of (e, Δϖ) values for a fixed semi-major axis,
/// coplanar with the perturber.
///
/// Uses Gauss-ring averaging of the exact 1/Δ (the multipole expansion does
/// not converge for the relevant α), with a softening of 0.01·a_p to
/// regularize crossing geometries.
///
/// Returns (e_vals, dvarpi_vals, portrait) where portrait[i][j] = H(e_i, Δϖ_j).
pub fn phase_portrait(
    a: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
    n_e: usize,
    n_dvarpi: usize,
) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let softening = 0.01 * a_p;
    let n_quad = 64;

    let e_vals: Vec<f64> = (0..n_e)
        .map(|i| (i as f64 + 0.5) / n_e as f64 * 0.95)
        .collect();

    let dvarpi_vals: Vec<f64> = (0..n_dvarpi)
        .map(|j| (j as f64) / n_dvarpi as f64 * 2.0 * PI - PI)
        .collect();

    let mut portrait = vec![vec![0.0; n_dvarpi]; n_e];

    for (i, &e) in e_vals.iter().enumerate() {
        for (j, &dv) in dvarpi_vals.iter().enumerate() {
            portrait[i][j] = numerical_secular_hamiltonian(
                a, e, 0.0, dv, 0.0, a_p, e_p, gm_p, n_quad, softening,
            );
        }
    }

    (e_vals, dvarpi_vals, portrait)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GM_SUN;

    const GM_P: f64 = 10.0 * 3.003489e-6 * GM_SUN; // 10 Earth masses

    #[test]
    fn test_precession_scales_linearly_with_mass() {
        // Two masses, same a, e and perturber geometry: rate scales linearly.
        let r1 = perturber_perihelion_precession(9.5371, 0.0, 3.0 * 3.003489e-6, 400.0, 0.0);
        let r2 = perturber_perihelion_precession(9.5371, 0.0, 9.0 * 3.003489e-6, 400.0, 0.0);
        assert!((r2 / r1 - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_precession_scales_as_inverse_cube_of_distance() {
        // Two perturber distances, same mass: ratio follows (a1/a2)^3.
        let m = 5.0 * 3.003489e-6;
        let r_near = perturber_perihelion_precession(9.5371, 0.0, m, 200.0, 0.0);
        let r_far = perturber_perihelion_precession(9.5371, 0.0, m, 400.0, 0.0);
        let cube = (400.0_f64 / 200.0).powi(3); // = 8
        assert!((r_near / r_far - cube).abs() < 1e-6 * cube);
    }

    #[test]
    fn test_eccentric_perturber_enhances_precession() {
        // (1 - e_p^2)^{-3/2} enhancement from perturber eccentricity.
        let m = 5.0 * 3.003489e-6;
        let r_c = perturber_perihelion_precession(9.5371, 0.0, m, 300.0, 0.0);
        let r_e = perturber_perihelion_precession(9.5371, 0.0, m, 300.0, 0.5);
        let factor = (1.0 - 0.5_f64 * 0.5).powf(-1.5);
        assert!((r_e / r_c - factor).abs() < 1e-12 * factor);
    }

    #[test]
    fn test_quadrupole_has_no_apsidal_dependence() {
        // The coplanar quadrupole must be a function of e only.
        let h1 = coplanar_quadrupole(100.0, 0.5, 700.0, 0.6, GM_P);
        let h2 = coplanar_quadrupole(100.0, 0.5, 700.0, 0.6, GM_P);
        assert_eq!(h1, h2);
        // And independent of omega at i = 0:
        let ha = quadrupole_hamiltonian(100.0, 0.5, 0.0, 0.3, 700.0, 0.6, GM_P);
        let hb = quadrupole_hamiltonian(100.0, 0.5, 0.0, 2.1, 700.0, 0.6, GM_P);
        assert!((ha - hb).abs() < 1e-30);
    }

    #[test]
    fn test_numerical_matches_quadrupole_at_small_alpha() {
        // Small α, circular perturber (octupole vanishes): the e-dependence
        // of the numerical average must match the quadrupole formula. The
        // constant monopole term −gm_p/a_p·(...) cancels in the difference.
        let a = 30.0;
        let a_p = 700.0;
        let e_p = 0.0;
        let n = 128;

        let num =
            |e: f64| numerical_secular_hamiltonian(a, e, 0.0, 0.0, 0.0, a_p, e_p, GM_P, n, 0.0);
        let quad = |e: f64| coplanar_quadrupole(a, e, a_p, e_p, GM_P);

        let d_num = num(0.6) - num(0.1);
        let d_quad = quad(0.6) - quad(0.1);
        let rel = ((d_num - d_quad) / d_quad).abs();
        // Hexadecapole correction is O(α²) ≈ 0.2%
        assert!(rel < 0.02, "numerical vs quadrupole ΔH mismatch: {rel:.3e}");
    }

    #[test]
    fn test_numerical_axisymmetric_for_circular_perturber() {
        // e_p = 0: the averaged perturber is an axisymmetric ring, so H must
        // not depend on Δϖ even at large α.
        let h0 =
            numerical_secular_hamiltonian(300.0, 0.4, 0.0, 0.0, 0.0, 700.0, 0.0, GM_P, 96, 0.0);
        let h1 = numerical_secular_hamiltonian(
            300.0,
            0.4,
            0.0,
            PI / 2.0,
            0.0,
            700.0,
            0.0,
            GM_P,
            96,
            0.0,
        );
        let rel = ((h0 - h1) / h0).abs();
        assert!(
            rel < 1e-10,
            "circular perturber should be axisymmetric: {rel:.3e}"
        );
    }

    #[test]
    fn test_numerical_apsidal_structure_with_eccentric_perturber() {
        // Eccentric perturber at moderate α: the exact average distinguishes
        // aligned from anti-aligned apsides (octupole+), which the quadrupole
        // cannot.
        let h_aligned =
            numerical_secular_hamiltonian(250.0, 0.3, 0.0, 0.0, 0.0, 700.0, 0.6, GM_P, 128, 0.0);
        let h_anti =
            numerical_secular_hamiltonian(250.0, 0.3, 0.0, PI, 0.0, 700.0, 0.6, GM_P, 128, 0.0);
        let rel = ((h_aligned - h_anti) / h_aligned).abs();
        assert!(rel > 1e-6, "expected Δϖ asymmetry, got {rel:.3e}");
    }

    #[test]
    fn test_numerical_quadrature_converged() {
        // n vs 2n agreement on a non-crossing geometry.
        let h64 =
            numerical_secular_hamiltonian(200.0, 0.5, 0.2, 1.0, 0.4, 700.0, 0.5, GM_P, 64, 0.0);
        let h128 =
            numerical_secular_hamiltonian(200.0, 0.5, 0.2, 1.0, 0.4, 700.0, 0.5, GM_P, 128, 0.0);
        let rel = ((h64 - h128) / h128).abs();
        assert!(rel < 1e-8, "quadrature not converged: {rel:.3e}");
    }

    #[test]
    fn test_phase_portrait_shape_and_finite() {
        let (e_vals, dv_vals, portrait) = phase_portrait(400.0, 700.0, 0.6, GM_P, 8, 12);
        assert_eq!(e_vals.len(), 8);
        assert_eq!(dv_vals.len(), 12);
        assert_eq!(portrait.len(), 8);
        assert_eq!(portrait[0].len(), 12);
        for row in &portrait {
            for &h in row {
                assert!(h.is_finite() && h < 0.0);
            }
        }
    }
}

/// Quadrupole secular coupling constant of two coplanar wires of masses
/// `m1`, `m2` (M☉) at semi-major axes `a_inner` < `a_outer` (AU), with outer
/// eccentricity factor `epsilon_outer` = √(1 − e²) (energy units,
/// M☉·AU²/day²):
///
///   C = G m₁ m₂ / (4 a_outer) · (a_inner/a_outer)² / ε_outer³
///
/// For a circular outer orbit ε = 1.
pub fn quadrupole_coupling_constant(
    m1: f64,
    m2: f64,
    a_inner: f64,
    a_outer: f64,
    epsilon_outer: f64,
) -> f64 {
    let ratio = a_inner / a_outer;
    G_AU3_MSUN_DAY2 * m1 * m2 / (4.0 * a_outer) * ratio * ratio / epsilon_outer.powi(3)
}

/// Laplace coefficient b_s^{(j)}(α), optionally softened:
///
///   b_s^{(j)}(α) = (1/π) ∫₀^{2π} cos(jψ) (1 − 2α cos ψ + α² + ε²)^{-s} dψ
///
/// with `softening` ε added in quadrature to the denominator distance
/// (ε = 0 recovers the classical coefficient, valid for α bounded away
/// from 1; a nonzero ε represents finite ring thickness / disk scale
/// height). Midpoint quadrature with 2048 nodes — spectrally converged for
/// the smooth kernel. For α → 0, b_{3/2}^{(1)} → 3α.
pub fn laplace_coefficient(s: f64, j: i32, alpha: f64, softening: f64) -> f64 {
    let n = 2048;
    let dpsi = 2.0 * PI / n as f64;
    let eps2 = softening * softening;
    let mut sum = 0.0;
    for k in 0..n {
        let psi = (k as f64 + 0.5) * dpsi;
        let denom = 1.0 - 2.0 * alpha * psi.cos() + alpha * alpha + eps2;
        sum += (j as f64 * psi).cos() / denom.powf(s);
    }
    sum * dpsi / PI
}

#[cfg(test)]
mod laplace_tests {
    use super::*;

    #[test]
    fn laplace_small_alpha_series() {
        // Series expansions (Murray & Dermott Eq. 6.69-6.70), small alpha:
        //   b_{1/2}^{(0)} ~ 2 + (1/2) alpha^2, b_{1/2}^{(1)} ~ alpha + (3/8) alpha^3,
        //   b_{3/2}^{(1)} ~ 3 alpha, b_{3/2}^{(2)} ~ (15/4) alpha^2.
        let alpha = 0.05;
        assert!(
            (laplace_coefficient(0.5, 0, alpha, 0.0) - (2.0 + 0.5 * alpha * alpha)).abs() < 1e-4
        );
        assert!(
            (laplace_coefficient(0.5, 1, alpha, 0.0) - (alpha + 0.375 * alpha.powi(3))).abs()
                < 1e-4
        );
        let a2 = 0.04;
        assert!((laplace_coefficient(1.5, 1, a2, 0.0) - 3.0 * a2).abs() < 5e-4);
        assert!((laplace_coefficient(1.5, 2, a2, 0.0) - 3.75 * a2 * a2).abs() < 1e-4);
    }

    #[test]
    fn laplace_softening_tames_divergence() {
        // Near alpha = 1 the bare b_{3/2} blows up; softening keeps it finite
        // and below the bare value.
        let bare = laplace_coefficient(1.5, 1, 0.98, 0.0);
        let soft = laplace_coefficient(1.5, 1, 0.98, 0.05);
        assert!(soft.is_finite() && soft < bare, "{soft} vs {bare}");
    }
}
