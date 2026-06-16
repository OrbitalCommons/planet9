//! Doubly-averaged secular Hamiltonian for the P9–KBO problem.
//!
//! Structure (Batygin & Morbidelli 2017, Eq. 1):
//!   H = H_J2 (giant-planet precession) + H_frame (rotating frame)
//!     + H_quad (standard Kozai-Lidov quadrupole, axisymmetric in Δϖ)
//!     + H_oct (octupole, the leading Δϖ-dependent term, ∝ e e₉)
//!
//! The quadrupole term is the standard doubly-averaged Legendre form from
//! `p9_core::analysis::secular` — it carries *no* Δϖ dependence (the
//! orbit-averaged perturber is an axisymmetric ring at this order). Apsidal
//! coupling first appears at octupole order with amplitude
//! ∝ (m₉/M)(α³/a₉)·e·e₉/(1−e₉²)^{5/2}·cos Δϖ, which vanishes for a circular
//! Planet Nine — the previous implementation carried a Δϖ term at quadrupole
//! scaling independent of e₉, predicting apsidal confinement even for e₉ = 0,
//! which is impossible.
//!
//! Giant planets are modeled as an effective J2 moment (Eq. 3).
//!
//! Caveat: the multipole expansion converges only for α = a/a₉ ≲ 0.3; for the
//! orbit-approaching geometries of the real problem use
//! `p9_core::analysis::secular::numerical_secular_hamiltonian`.

use p9_core::analysis::secular::quadrupole_hamiltonian;
use p9_core::constants::{
    EARTH_MASS_SOLAR, GM_SUN, MASS_JUPITER_SOLAR, MASS_NEPTUNE_SOLAR, MASS_SATURN_SOLAR,
    MASS_URANUS_SOLAR, TWO_PI, YEAR_DAYS,
};
use p9_core::units::{
    au, days, gm_from_au3_day2, radians, solar_masses, AngularVelocity, GravitationalParameter,
    Length, Mass,
};
use serde::{Deserialize, Serialize};

/// Parameters for the doubly-averaged secular Hamiltonian.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SecularHamiltonianParams {
    /// Planet Nine semi-major axis (AU)
    pub a9: f64,
    /// Planet Nine eccentricity
    pub e9: f64,
    /// Planet Nine mass in solar masses
    pub m9_solar: f64,
    /// Effective J2 coefficient from giant planets (AU^2)
    pub j2_eff: f64,
    /// Planet Nine apsidal precession rate (rad/day; see `precession_rate_j2`
    /// for the matching unit convention)
    pub precession_rate_9: f64,
}

impl SecularHamiltonianParams {
    /// Fiducial parameters from the paper: a_9=700, e_9=0.6, m_9=10 ME.
    pub fn default_paper() -> Self {
        Self {
            a9: 700.0,
            e9: 0.6,
            m9_solar: 10.0 * EARTH_MASS_SOLAR,
            j2_eff: compute_j2_effective(),
            precession_rate_9: 0.0,
        }
    }

    /// Alternative parameter set: a_9=600, e_9=0.5, m_9=10 ME.
    pub fn alternative_600() -> Self {
        Self {
            a9: 600.0,
            e9: 0.5,
            m9_solar: 10.0 * EARTH_MASS_SOLAR,
            j2_eff: compute_j2_effective(),
            precession_rate_9: 0.0,
        }
    }

    /// Planet Nine GM in AU³/day².
    pub fn gm9(&self) -> f64 {
        self.m9_solar * GM_SUN
    }

    /// Planet Nine semi-major axis as a dimension-checked [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a9)
    }

    /// Planet Nine mass as a dimension-checked [`Mass`] (stored in solar masses).
    pub fn mass(&self) -> Mass {
        solar_masses(self.m9_solar)
    }

    /// Planet Nine gravitational parameter as a dimension-checked
    /// [`GravitationalParameter`] (AU³/day² source).
    pub fn gm9_typed(&self) -> GravitationalParameter {
        gm_from_au3_day2(self.gm9())
    }

    /// Planet Nine apsidal precession rate as a dimension-checked
    /// [`AngularVelocity`] (the stored value is in rad/day).
    pub fn precession_rate_9_typed(&self) -> AngularVelocity {
        (radians(self.precession_rate_9) / days(1.0)).into()
    }
}

/// Compute the effective J2 coefficient from the four giant planets (Eq. 3).
///
/// J2_eff = (1/2) sum_j (m_j/M_sun) * a_j^2
/// This captures the quadrupolar potential of the giant planet system.
pub fn compute_j2_effective() -> f64 {
    let giants: [(f64, f64); 4] = [
        (MASS_JUPITER_SOLAR, 5.203),
        (MASS_SATURN_SOLAR, 9.537),
        (MASS_URANUS_SOLAR, 19.189),
        (MASS_NEPTUNE_SOLAR, 30.070),
    ];

    giants.iter().map(|(m, a)| 0.5 * m * a * a).sum()
}

/// Perihelion precession rate from the giant planet J2 potential (Eq. 3).
///
/// dot_varpi_J2 = (3/2) * n * J2_eff / (a^2 * (1-e^2)^2)
///
/// Returns **rad/day** (the mean motion n carries AU³/day² units from
/// GM_SUN). Use `precession_rate_j2_per_year` for rad/yr. A previous version
/// documented rad/yr while returning rad/day.
pub fn precession_rate_j2(a: f64, e: f64, j2_eff: f64) -> f64 {
    let n = (GM_SUN / (a * a * a)).sqrt(); // mean motion (rad/day)
    let eta_sq = (1.0 - e * e) * (1.0 - e * e);
    1.5 * n * j2_eff / (a * a * eta_sq)
}

/// Perihelion precession rate from the giant-planet J2 potential in rad/yr.
pub fn precession_rate_j2_per_year(a: f64, e: f64, j2_eff: f64) -> f64 {
    precession_rate_j2(a, e, j2_eff) * YEAR_DAYS
}

/// Octupole (leading Δϖ-dependent) term of the doubly-averaged P9–KBO
/// interaction for a coplanar-to-low-inclination configuration:
///
///   H_oct = −(15/16) (gm₉ α³ / a₉) · e e₉ (1 + 3e²/4) / (1−e₉²)^{5/2} · cos Δϖ
///
/// (Ford, Kozinsky & Rasio 2000, coplanar test-particle limit). The
/// amplitude is ∝ e·e₉ and vanishes as e₉ → 0: a circular Planet Nine cannot
/// produce apsidal confinement.
pub fn octupole_term(a: f64, e: f64, delta_varpi: f64, params: &SecularHamiltonianParams) -> f64 {
    let alpha = a / params.a9;
    let eta9_sq = 1.0 - params.e9 * params.e9;
    -(15.0 / 16.0)
        * (params.gm9() * alpha.powi(3) / params.a9)
        * e
        * params.e9
        * (1.0 + 0.75 * e * e)
        / eta9_sq.powf(2.5)
        * delta_varpi.cos()
}

/// Secular Hamiltonian value at given test particle orbital elements
/// (AU²/day² energy units).
///
/// H = H_J2(precession) + H_frame(rotating) + H_quad + H_oct
///
/// Inputs are test particle elements: a (AU), e, i (rad), omega (rad,
/// argument of perihelion from the common node with P9's orbit plane),
/// delta_varpi = varpi - varpi_9 (rad).
pub fn secular_hamiltonian(
    a: f64,
    e: f64,
    i: f64,
    omega: f64,
    delta_varpi: f64,
    params: &SecularHamiltonianParams,
) -> f64 {
    let eta = (1.0 - e * e).sqrt();

    // Term 1: Giant-planet J2 precession
    let h_j2 = -1.5 * GM_SUN * params.j2_eff / (a * a * a * eta * eta * eta);

    // Term 2: Rotating frame contribution (linear in the action √(GMa)·η)
    let h_frame = -params.precession_rate_9 * (GM_SUN * a).sqrt() * eta;

    // Term 3: standard doubly-averaged quadrupole (axisymmetric in Δϖ).
    let h_quad = quadrupole_hamiltonian(a, e, i, omega, params.a9, params.e9, params.gm9());

    // Term 4: octupole — the leading Δϖ-dependent coupling, ∝ e e₉.
    let h_oct = octupole_term(a, e, delta_varpi, params);

    h_j2 + h_frame + h_quad + h_oct
}

/// High-inclination secular resonance angle (theta).
///
/// theta = 2*Omega - varpi - varpi_9
/// Conjugate action: Theta = sqrt(1-e^2)/2 * (1 - cos(i))
pub fn high_inclination_angle(omega_big: f64, varpi: f64, varpi_9: f64) -> f64 {
    let theta = 2.0 * omega_big - varpi - varpi_9;
    ((theta % TWO_PI) + TWO_PI) % TWO_PI
}

/// Conjugate action for the high-inclination secular resonance.
pub fn high_inclination_action(e: f64, i: f64) -> f64 {
    (1.0 - e * e).sqrt() / 2.0 * (1.0 - i.cos())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::constants::DEG2RAD;

    #[test]
    fn typed_secular_param_accessors_match_f64() {
        let p = SecularHamiltonianParams::default_paper();
        assert_relative_eq!(
            (p.semi_major_axis() / au(1.0)).value,
            p.a9,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            (p.mass() / solar_masses(1.0)).value,
            p.m9_solar,
            max_relative = 1e-9
        );
        assert_relative_eq!(
            (p.gm9_typed() / gm_from_au3_day2(1.0)).value,
            p.gm9(),
            max_relative = 1e-9
        );
    }

    #[test]
    fn test_j2_effective_positive() {
        let j2 = compute_j2_effective();
        assert!(j2 > 0.0, "J2_eff = {}", j2);
        // Should be dominated by Jupiter: ~0.5 * 9.5e-4 * 27 ~ 0.013
        assert!(
            j2 > 0.01 && j2 < 0.1,
            "J2_eff = {} out of expected range",
            j2
        );
    }

    #[test]
    fn test_precession_rate_units() {
        let j2 = compute_j2_effective();
        let per_day = precession_rate_j2(300.0, 0.7, j2);
        let per_year = precession_rate_j2_per_year(300.0, 0.7, j2);
        assert!(per_day > 0.0);
        assert!((per_year / per_day - YEAR_DAYS).abs() < 1e-12);
        // Sanity scale: ϖ precession at a = 300 AU is ~ tens of arcsec/cycle —
        // far below one rad/yr.
        assert!(per_year < 1e-3, "rate = {per_year} rad/yr");
    }

    #[test]
    fn test_fiducial_params() {
        let params = SecularHamiltonianParams::default_paper();
        assert!((params.a9 - 700.0).abs() < 0.01);
        assert!((params.e9 - 0.6).abs() < 0.01);
        assert!((params.m9_solar - 10.0 * EARTH_MASS_SOLAR).abs() < 1e-12);
    }

    #[test]
    fn test_secular_hamiltonian_finite() {
        let params = SecularHamiltonianParams::default_paper();
        let h = secular_hamiltonian(300.0, 0.7, 20.0 * DEG2RAD, 0.0, 0.5, &params);
        assert!(h.is_finite(), "Hamiltonian should be finite, got {}", h);
    }

    #[test]
    fn test_apsidal_coupling_vanishes_for_circular_p9() {
        // The defining physics fix (H2): the Δϖ-dependent part of H must
        // vanish as e₉ → 0 — a circular P9 averages to an axisymmetric ring
        // and cannot confine apsides.
        let mut params = SecularHamiltonianParams::default_paper();
        params.e9 = 0.0;
        let (a, e, i, omega) = (300.0, 0.7, 10.0 * DEG2RAD, 0.4);
        let h0 = secular_hamiltonian(a, e, i, omega, 0.0, &params);
        let hpi = secular_hamiltonian(a, e, i, omega, std::f64::consts::PI, &params);
        assert!(
            (h0 - hpi).abs() < 1e-30,
            "Δϖ dependence must vanish at e9 = 0: ΔH = {:.3e}",
            h0 - hpi
        );

        // And it must scale linearly with e₉ for small e₉.
        let dv_amplitude = |e9: f64| {
            let mut p = SecularHamiltonianParams::default_paper();
            p.e9 = e9;
            let h0 = secular_hamiltonian(a, e, i, omega, 0.0, &p);
            let hpi = secular_hamiltonian(a, e, i, omega, std::f64::consts::PI, &p);
            (h0 - hpi).abs()
        };
        let a1 = dv_amplitude(0.1);
        let a2 = dv_amplitude(0.2);
        let ratio = a2 / a1;
        // (1−e9²)^{-5/2} grows slightly faster than linear; allow 2.0–2.2.
        assert!(
            (1.9..2.3).contains(&ratio),
            "octupole amplitude should scale ~linearly in e9, ratio = {ratio:.3}"
        );
    }

    #[test]
    fn test_quadrupole_part_axisymmetric_in_delta_varpi() {
        // With the octupole removed (e9 = 0 handles that), the remaining H is
        // Δϖ-independent; with e9 > 0 the *only* Δϖ dependence is ∝ cos Δϖ.
        let params = SecularHamiltonianParams::default_paper();
        let (a, e, i, omega) = (250.0, 0.6, 5.0 * DEG2RAD, 1.1);
        let h0 = secular_hamiltonian(a, e, i, omega, 0.0, &params);
        let h_pi2 = secular_hamiltonian(a, e, i, omega, std::f64::consts::FRAC_PI_2, &params);
        let h_pi = secular_hamiltonian(a, e, i, omega, std::f64::consts::PI, &params);
        // cos(π/2) = 0 → H(π/2) must be the mean of H(0) and H(π).
        assert!(
            (h_pi2 - 0.5 * (h0 + h_pi)).abs() < 1e-25 * h0.abs().max(1.0),
            "Δϖ dependence should be pure cos Δϖ"
        );
    }

    #[test]
    fn test_high_inclination_action_zero_at_zero_i() {
        let theta = high_inclination_action(0.5, 0.0);
        assert!(
            theta.abs() < 1e-10,
            "action at i=0 should be ~0, got {}",
            theta
        );
    }

    #[test]
    fn test_high_inclination_action_positive() {
        let theta = high_inclination_action(0.5, 30.0 * DEG2RAD);
        assert!(theta > 0.0, "action should be positive at i=30deg");
    }
}
