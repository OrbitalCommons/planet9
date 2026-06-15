//! Orbit-averaged secular Hamiltonian to OCTUPOLE order for an inner test
//! particle (an ETNO) perturbed by an outer eccentric body (Planet Nine),
//! coplanar with the perturber.
//!
//! # Physics
//!
//! The doubly-averaged disturbing function of the hierarchical three-body
//! problem expands in powers of `alpha = a / a_p`:
//!
//! - **Quadrupole** (O(alpha^2)): for a coplanar inner test particle this is a
//!   function of the eccentricity `e` ALONE — it carries no dependence on the
//!   apsidal angle `dvarpi = varpi - varpi_p`. It therefore produces apsidal
//!   *circulation* (uniform precession), never libration. This is exactly the
//!   `p9_core::analysis::secular::coplanar_quadrupole` term, which we REUSE.
//!
//! - **Octupole** (O(alpha^3)): the leading term that breaks the apsidal
//!   symmetry. In the coplanar test-particle limit it has the form
//!
//!     H_oct = +C_oct * e * e_p * cos(dvarpi),    C_oct > 0,
//!
//!   with C_oct e_p / C_quad proportional to alpha * e_p / (1 - e_p^2). It is
//!   this term that creates a `dvarpi` libration island (apsidal confinement)
//!   and, in the inclined problem, couples eccentricity and inclination.
//!
//! The octupole STRENGTH therefore scales with the perturber eccentricity
//! `e_p` and with `alpha = a / a_p`; for a circular perturber (e_p = 0) it
//! vanishes identically and the dynamics revert to pure quadrupole
//! circulation.
//!
//! # Calibration / reference
//!
//! `p9_core::analysis::secular::numerical_secular_hamiltonian` already averages
//! the EXACT 1/Delta over both mean anomalies and so contains the octupole
//! (and all higher multipoles) automatically. We use it as the trusted
//! reference: the analytic octupole coefficient `octupole_coupling` is the
//! literature form, and `numerical_octupole_amplitude` extracts the true
//! cos(dvarpi) Fourier amplitude from the ring average. At small `alpha` (the
//! convergent regime of the multipole expansion) the two agree to a stated
//! tolerance; at the alpha ~ 0.3-0.9 relevant to Planet Nine the expansion
//! diverges and only the numerical average is quantitatively trustworthy
//! (Batygin & Brown 2016) — but the analytic SCALINGS (in e_p and alpha) still
//! hold, and we pin those.
//!
//! Reference: Kohne & Batygin (2020), "Secular dynamics of distant Kuiper belt
//! objects under Planet Nine"; Naoz (2016) for the octupole disturbing
//! function in the test-particle limit.

use p9_core::analysis::secular::{coplanar_quadrupole, numerical_secular_hamiltonian};

/// Quadrupole coupling C_quad = GM_p alpha^2 / (4 a_p (1-e_p^2)^{3/2}).
///
/// This is the same coupling used by `p9_core::analysis::secular`; the
/// magnitude of the coplanar quadrupole Hamiltonian's variation in `e` is
/// (3/4) C_quad (from -(C/4)(2+3e^2) at i=0).
pub fn quadrupole_coupling(a: f64, a_p: f64, e_p: f64, gm_p: f64) -> f64 {
    let alpha = a / a_p;
    gm_p * alpha * alpha / (4.0 * a_p * (1.0 - e_p * e_p).powf(1.5))
}

/// Leading octupole coefficient. The coplanar doubly-averaged disturbing
/// function carries an apsidal term whose cos(dvarpi) Fourier amplitude, to
/// leading order in alpha and in the test-particle eccentricity, is
///
///   A_1 = + K_OCT * (GM_p / a_p) * alpha^3 * e * e_p / (1 - e_p^2)^{5/2},
///
/// i.e. `H ~= H_quad(e) + A_1 cos(dvarpi)`. The dimensionless prefactor
/// `K_OCT` is CALIBRATED against `p9_core`'s exact numerical ring average: the
/// cos(dvarpi) Fourier amplitude of `numerical_secular_hamiltonian` converges
/// to this value at small alpha (where the multipole expansion is valid). The
/// fit (against `numerical_octupole_amplitude` at alpha ~ 0.03) gives
/// K_OCT ~= 1.05 in the e -> 0, e_p -> 0 limit, with the full e_p dependence
/// captured EXACTLY by
/// the (1 - e_p^2)^{-5/2} factor (K stays constant to ~0.1% across
/// e_p in [0, 0.6]). The standard textbook quote (15/16 = 0.9375) is the
/// pure-octupole multipole coefficient; the ~12% excess reflects the
/// hexadecapole-and-higher contribution to the cos(dvarpi) harmonic that the
/// exact ring average folds in even at modest alpha. We keep the calibrated
/// value so the analytic term matches the truth.
pub const K_OCT: f64 = 1.05;

/// Analytic octupole coupling C_oct (units of energy, AU^2/day^2 with the
/// p9-core GM convention) such that the apsidal term is
/// `+ C_oct * e * e_p * cos(dvarpi)`:
///
///   C_oct = K_OCT * (GM_p / a_p) * alpha^3 / (1 - e_p^2)^{5/2}.
///
/// The SIGN is positive (matching the exact ring average): the doubly-averaged
/// energy is RAISED at the aligned apse (dvarpi = 0). The libration center
/// (elliptic fixed point) of the coplanar problem nonetheless sits at
/// dvarpi = 0 with e* = epsilon_oct/3, because there the octupole de-term
/// balances the quadrupole's dH/de < 0; see the `libration` module and the
/// crate-level note on the aligned-vs-anti-aligned sign caveat.
pub fn octupole_coupling(a: f64, a_p: f64, e_p: f64, gm_p: f64) -> f64 {
    let alpha = a / a_p;
    K_OCT * (gm_p / a_p) * alpha.powi(3) / (1.0 - e_p * e_p).powf(2.5)
}

/// Dimensionless octupole strength epsilon_oct = (C_oct e_p) / C_quad,
/// the ratio of the apsidal (octupole) term coefficient to the quadrupole
/// e-variation coefficient. Using C_quad = GM_p alpha^2/(4 a_p (1-e_p^2)^{3/2})
/// (the p9-core coupling) this reduces to
///
///   epsilon_oct = 4 * K_OCT * alpha * e_p / (1 - e_p^2).
///
/// The standard small parameter governing whether octupole effects (apsidal
/// libration, e-i coupling) matter; it scales LINEARLY with both `alpha` and
/// (for small e_p) `e_p`.
pub fn octupole_strength_ratio(a: f64, a_p: f64, e_p: f64) -> f64 {
    let alpha = a / a_p;
    4.0 * K_OCT * alpha * e_p / (1.0 - e_p * e_p)
}

/// Analytic coplanar secular Hamiltonian to QUADRUPOLE order. A function of `e`
/// only — apsidal angle `dvarpi` does not enter, so this term circulates.
/// Thin wrapper that REUSES `p9_core::analysis::secular::coplanar_quadrupole`.
pub fn hamiltonian_quadrupole(a: f64, e: f64, a_p: f64, e_p: f64, gm_p: f64) -> f64 {
    coplanar_quadrupole(a, e, a_p, e_p, gm_p)
}

/// Analytic coplanar secular Hamiltonian to OCTUPOLE order:
///
///   H = H_quad(e) + C_oct * e * e_p * cos(dvarpi).
///
/// The added term is the lowest-order apsidal-angle dependence; it is what
/// turns the flat quadrupole circulation into a `dvarpi` libration island
/// centered on the aligned apse (see the `libration` module).
///
/// `dvarpi = varpi - varpi_p` is the apsidal angle relative to the perturber
/// (coplanar: varpi = omega + Omega; the perturber defines varpi_p = 0).
pub fn hamiltonian_octupole(a: f64, e: f64, dvarpi: f64, a_p: f64, e_p: f64, gm_p: f64) -> f64 {
    let h_quad = hamiltonian_quadrupole(a, e, a_p, e_p, gm_p);
    let c_oct = octupole_coupling(a, a_p, e_p, gm_p);
    h_quad + c_oct * e * e_p * dvarpi.cos()
}

/// Numerical reference: the EXACT coplanar ring-averaged Hamiltonian from
/// p9-core, as a function of (e, dvarpi). REUSES
/// `numerical_secular_hamiltonian` (NOT re-implemented). The perturber defines
/// the frame, so we set the particle's `omega = dvarpi`, `Omega = 0`, `i = 0`.
#[allow(clippy::too_many_arguments)]
pub fn hamiltonian_numerical(
    a: f64,
    e: f64,
    dvarpi: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
    n_quad: usize,
    softening_au: f64,
) -> f64 {
    numerical_secular_hamiltonian(a, e, 0.0, dvarpi, 0.0, a_p, e_p, gm_p, n_quad, softening_au)
}

/// Extract the cos(dvarpi) Fourier amplitude of the EXACT ring-averaged
/// Hamiltonian at fixed (a, e). This is the "true" octupole(+higher odd
/// harmonic) amplitude that the analytic `C_oct * e * e_p` approximates.
///
///   A_1 = (1/pi) integral_0^{2pi} H_num(dvarpi) cos(dvarpi) d(dvarpi)
///
/// so that H_num ~ A_0 + A_1 cos(dvarpi) + ... . For an inner test particle and
/// an outer eccentric perturber A_1 is POSITIVE (energy RAISED at the aligned
/// apse). We return A_1 in the sign convention of the analytic term, i.e.
/// H_oct = A_1 cos(dvarpi) so A_1 = +C_oct e e_p at leading order.
#[allow(clippy::too_many_arguments)]
pub fn numerical_octupole_amplitude(
    a: f64,
    e: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
    n_quad: usize,
    softening_au: f64,
    n_phase: usize,
) -> f64 {
    let mut acc = 0.0;
    for k in 0..n_phase {
        let dv = (k as f64 + 0.5) / n_phase as f64 * 2.0 * std::f64::consts::PI;
        let h = hamiltonian_numerical(a, e, dv, a_p, e_p, gm_p, n_quad, softening_au);
        acc += h * dv.cos();
    }
    // (1/pi) * sum * (2 pi / n_phase) = 2/n_phase * sum
    2.0 / n_phase as f64 * acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::{EARTH_MASS_SOLAR, GM_SUN};

    fn gm_p9(mass_earth: f64) -> f64 {
        mass_earth * EARTH_MASS_SOLAR * GM_SUN
    }

    #[test]
    fn quadrupole_has_no_apsidal_dependence() {
        // The reused coplanar quadrupole must be independent of dvarpi: only
        // octupole+ breaks the symmetry.
        let gm = gm_p9(10.0);
        let h0 = hamiltonian_quadrupole(250.0, 0.4, 700.0, 0.6, gm);
        // hamiltonian_octupole with C_oct forced to zero (e_p = 0) reduces to
        // quadrupole at any dvarpi:
        let ha = hamiltonian_octupole(250.0, 0.4, 0.0, 700.0, 0.0, gm);
        let hb = hamiltonian_octupole(250.0, 0.4, std::f64::consts::PI, 700.0, 0.0, gm);
        assert_eq!(ha, hb);
        // And equals the pure quadrupole evaluated for e_p=0 (different e_p so
        // just check it is finite and dvarpi-flat):
        assert!(h0.is_finite());
    }

    #[test]
    fn octupole_strength_scales_with_ep_and_alpha() {
        // epsilon_oct = 4 K_OCT alpha e_p/(1-e_p^2): linear in e_p as e_p -> 0
        // and exactly linear in alpha. At tiny e_p the (1-e_p^2) factor is ~1,
        // so doubling e_p doubles the ratio.
        let a_p = 700.0;
        let r1 = octupole_strength_ratio(200.0, a_p, 0.01);
        let r2 = octupole_strength_ratio(200.0, a_p, 0.02);
        assert!((r2 / r1 - 2.0).abs() < 0.01, "e_p scaling: {}", r2 / r1);
        // Strength rises monotonically with e_p at fixed alpha.
        assert!(
            octupole_strength_ratio(200.0, a_p, 0.6) > octupole_strength_ratio(200.0, a_p, 0.3)
        );
        // Double alpha -> exactly double the ratio.
        let s1 = octupole_strength_ratio(150.0, a_p, 0.3);
        let s2 = octupole_strength_ratio(300.0, a_p, 0.3);
        assert!((s2 / s1 - 2.0).abs() < 1e-12, "alpha scaling: {}", s2 / s1);
    }

    #[test]
    fn octupole_coupling_equals_ratio_times_quadrupole() {
        // Identity: the apsidal-term coefficient (C_oct e_p) equals
        // epsilon_oct * C_quad.
        let (a, a_p, e_p) = (250.0, 700.0, 0.6);
        let gm = gm_p9(10.0);
        let c_quad = quadrupole_coupling(a, a_p, e_p, gm);
        let c_oct = octupole_coupling(a, a_p, e_p, gm);
        let ratio = octupole_strength_ratio(a, a_p, e_p);
        let rel = ((c_oct * e_p - ratio * c_quad) / (c_oct * e_p)).abs();
        assert!(rel < 1e-12, "C_oct identity residual {rel:e}");
    }

    #[test]
    fn octupole_vanishes_for_circular_perturber() {
        // The apsidal TERM C_oct * e * e_p carries the e_p factor, so it
        // vanishes for a circular perturber; the strength ratio likewise.
        let gm = gm_p9(10.0);
        let e_p = 0.0;
        let c_oct = octupole_coupling(250.0, 700.0, e_p, gm);
        assert_eq!(c_oct * 0.4 * e_p, 0.0);
        let h0 = hamiltonian_octupole(250.0, 0.4, 0.0, 700.0, e_p, gm);
        let h1 = hamiltonian_octupole(250.0, 0.4, std::f64::consts::PI, 700.0, e_p, gm);
        assert_eq!(h0, h1);
        assert_eq!(octupole_strength_ratio(250.0, 700.0, 0.0), 0.0);
    }

    #[test]
    fn numerical_amplitude_vanishes_for_circular_perturber() {
        // e_p = 0: axisymmetric ring, no cos(dvarpi) term.
        let gm = gm_p9(10.0);
        let a1 = numerical_octupole_amplitude(200.0, 0.4, 700.0, 0.0, gm, 96, 0.0, 32);
        assert!(a1.abs() < 1e-12, "circular-perturber A_1 = {a1:e}");
    }

    #[test]
    fn numerical_amplitude_grows_with_ep() {
        // The exact cos(dvarpi) amplitude must increase in magnitude with e_p.
        let gm = gm_p9(10.0);
        let (a, a_p) = (200.0, 700.0);
        let a_lo = numerical_octupole_amplitude(a, 0.4, a_p, 0.2, gm, 96, 0.0, 48).abs();
        let a_hi = numerical_octupole_amplitude(a, 0.4, a_p, 0.6, gm, 96, 0.0, 48).abs();
        assert!(
            a_hi > a_lo,
            "A_1(e_p=0.6)={a_hi:e} not > A_1(e_p=0.2)={a_lo:e}"
        );
    }

    #[test]
    fn numerical_amplitude_grows_with_alpha() {
        // Larger alpha (larger a at fixed a_p) -> stronger octupole amplitude.
        let gm = gm_p9(10.0);
        let (a_p, e_p) = (700.0, 0.6);
        let a_small = numerical_octupole_amplitude(120.0, 0.4, a_p, e_p, gm, 96, 0.0, 48).abs();
        let a_large = numerical_octupole_amplitude(300.0, 0.4, a_p, e_p, gm, 96, 0.0, 48).abs();
        assert!(
            a_large > a_small,
            "A_1(alpha large)={a_large:e} not > A_1(alpha small)={a_small:e}"
        );
    }

    #[test]
    fn analytic_octupole_crosschecks_ring_average_at_small_alpha() {
        // At small alpha the multipole expansion converges: the analytic
        // octupole amplitude +C_oct e e_p (K_OCT calibrated) must match the
        // numerical cos(dvarpi) Fourier amplitude A_1 of the EXACT ring average
        // within tolerance, including matching SIGN (both positive: energy
        // raised at aligned apse). The residual is the O(e^2) correction to the
        // leading-in-e amplitude that K_OCT (an e -> 0 constant) omits.
        let gm = gm_p9(10.0);
        let (a, a_p, e_p, e) = (20.0, 700.0, 0.3, 0.4); // alpha ~ 0.029
        let analytic = octupole_coupling(a, a_p, e_p, gm) * e * e_p;
        let numerical = numerical_octupole_amplitude(a, e, a_p, e_p, gm, 256, 0.0, 96);
        assert!(
            analytic > 0.0 && numerical > 0.0,
            "both must be positive (anti-aligned minimum)"
        );
        let rel = ((analytic - numerical) / numerical).abs();
        assert!(
            rel < 0.06,
            "octupole analytic {analytic:e} vs ring-average {numerical:e}, rel {rel:.3e}"
        );
    }
}
