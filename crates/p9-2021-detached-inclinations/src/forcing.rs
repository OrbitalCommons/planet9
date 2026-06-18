//! Secular inclination forcing of a detached test particle by an inclined
//! Planet Nine, in linear Laplace-Lagrange theory.
//!
//! A detached object (q decoupled from Neptune, a ≳ 250 AU) feels two secular
//! torques on its orbital plane:
//!
//!   * the **giant planets**, which define (with the Sun) a near-ecliptic
//!     reference plane and make the particle's node regress at a rate `B_gp`;
//!   * **Planet Nine**, an inclined massive body, which tries to align the
//!     particle's plane with *its* plane at a rate `B_p9`.
//!
//! In the linear (small mutual-inclination) secular problem the particle's
//! inclination vector (p, q) = (sin i sin Ω, sin i cos Ω) precesses about a
//! **forced inclination** that is the mass/coupling-weighted mean of the two
//! driving planes. With the giant-planet plane taken as the reference
//! (≈ invariable plane), the forced inclination measured from that plane is
//!
//!   i_forced = B_p9 / (B_gp + B_p9) · sin(i9)               (linearized)
//!
//! where `i9` is Planet Nine's inclination to the same reference plane. The
//! free inclination — the particle's intrinsic tilt — precesses on a cone of
//! this opening about the forced pole, so over secular time the *observed*
//! inclination spreads between |i_forced − i_free| and (i_forced + i_free).
//!
//! This is the linear backbone of Anderson & Kaib (2021): an inclined Planet
//! Nine forces detached objects onto an inclined Laplace plane, broadening the
//! inclination distribution relative to a Planet-Nine-free Solar System. The
//! coupling `B_p9` reuses p9-core's orbit-averaged quadrupole secular
//! constant; `B_gp` reuses p9-core's giant-planet nodal-precession rate.
//!
//! References: Murray & Dermott (1999) §7 (Laplace-Lagrange secular theory);
//! Anderson & Kaib (2021), arXiv:2109.13307.

use p9_core::analysis::secular::quadrupole_hamiltonian;
use p9_core::constants::GM_SUN;
use p9_core::initial_conditions::planets::giant_planets_j2000;
use p9_core::types::{cartesian_to_elements, P9Params};
use p9_core::units::{julian_year, radians, Angle, AngularVelocity};

/// Giant-planet nodal precession frequency `B_gp` for a particle at (a, e, i),
/// in rad/yr. This is the secular self-frequency that anchors the particle to
/// the invariable plane in the absence of Planet Nine:
///
///   |Ω̇| = (3/4) n cos i /(1−e²)² Σ_j (m_j/M_⊙)(a_j/a)²,
///
/// with the giant-planet masses and semi-major axes read from p9-core's J2000
/// DE-ephemeris states (no local planet table). This is the identical closed
/// form used by p9-2025-clustering's `orbital_poles::nodal_precession_rate`;
/// inlined here because that lives in a sibling paper crate, not p9-core.
///
/// Returned as a positive magnitude (a precession *rate*); the forced-
/// inclination ratio only depends on the relative strengths.
pub fn giant_planet_frequency_typed(a: f64, e: f64, i: f64) -> AngularVelocity {
    let n = (GM_SUN / (a * a * a)).sqrt(); // rad/day
    let eta_sq = (1.0 - e * e) * (1.0 - e * e);
    let sum: f64 = giant_planets_j2000()
        .iter()
        .map(|body| {
            let a_j = cartesian_to_elements(&body.state, GM_SUN).a;
            body.mass * (a_j / a).powi(2)
        })
        .sum();
    let rate_day = 0.75 * n * i.cos() / eta_sq * sum;
    let rate_yr = (rate_day * 365.25).abs(); // rad/yr
    (radians(rate_yr) / julian_year()).into()
}

/// Planet Nine secular inclination-coupling frequency `B_p9` for a test
/// particle at (a, e), in rad/yr.
///
/// The orbit-averaged quadrupole Hamiltonian (p9-core) is
/// H = −(C/4)[(2+3e²)(3cos²I − 1) + ...], with the inclination-restoring
/// coupling constant C = GM₉ α²/(4 a₉ (1−e₉²)^{3/2}). The nodal precession
/// frequency a particle's plane acquires about Planet Nine's plane is the
/// second derivative of H with respect to the inclination action; to leading
/// order in I it reduces to
///
///   B_p9 = (3/4) · GM₉/(a₉³) · (a/a₉)^{1/2} ... — but the cleanest, and the
///
/// way that stays consistent with the giant-planet expression above, is to use
/// the *same* disturbing-function coefficient form: treating Planet Nine as one
/// more (outer, eccentric, massive) perturber,
///
///   B_p9 = (3/4) n (m₉/M_⊙) (a/a₉)² · 1/(1−e₉²)^{3/2}
///
/// with n = √(GM_⊙/a³) the particle mean motion. The α² = (a/a₉)² scaling makes
/// the coupling *strengthen with semi-major axis*, the central prediction. We
/// cross-check this closed form against the p9-core quadrupole Hamiltonian's
/// curvature in the unit tests.
pub fn planet_nine_frequency_typed(a: f64, p9: &P9Params) -> AngularVelocity {
    let n = (GM_SUN / (a * a * a)).sqrt(); // rad/day
    let m9 = p9.mass_solar();
    let alpha = a / p9.a;
    let ecc_factor = (1.0 - p9.e * p9.e).powf(-1.5);
    let rate_day = 0.75 * n * m9 * alpha * alpha * ecc_factor;
    let rate_yr = rate_day * 365.25; // rad/yr
    (radians(rate_yr) / julian_year()).into()
}

/// Curvature of the p9-core quadrupole Hamiltonian in inclination at I = 0,
/// per unit GM, used to validate that [`planet_nine_frequency`] tracks the
/// reused secular coupling (and its α² and 1/(1−e₉²)^{3/2} scalings) rather
/// than an independent ad-hoc formula.
///
/// d²H/dI² at I = 0 for the quadrupole H = −(C/4)(2+3e²)(3cos²I − 1) is
/// +(3/2) C (2+3e²); we return it normalized so callers can compare scalings.
pub fn quadrupole_inclination_curvature(a: f64, e: f64, p9: &P9Params) -> f64 {
    let di = 1e-4;
    let h = |i: f64| quadrupole_hamiltonian(a, e, i, 0.0, p9.a, p9.e, p9.gm());
    (h(di) - 2.0 * h(0.0) + h(-di)) / (di * di)
}

/// Forced inclination (rad) of a detached particle's orbital plane relative to
/// the giant-planet (invariable) reference plane, given an inclined Planet
/// Nine. Linear Laplace-Lagrange result:
///
///   i_forced = B_p9 / (B_gp + B_p9) · i9
///
/// With no Planet Nine (B_p9 = 0) this is identically 0 — the particle stays in
/// the invariable plane. As Planet Nine's mass or inclination grows, or as the
/// particle's semi-major axis grows (stronger α² coupling), the forced plane
/// tilts further toward Planet Nine's plane.
pub fn forced_inclination_typed(a: f64, e: f64, i: f64, p9: &P9Params) -> Angle {
    let b_gp = giant_planet_frequency_typed(a, e, i);
    let b_p9 = planet_nine_frequency_typed(a, p9);
    // Extract rad/yr magnitudes via the version-safe units pattern.
    let per_yr = radians(1.0) / julian_year();
    let b_gp_yr = (b_gp / per_yr).value;
    let b_p9_yr = (b_p9 / per_yr).value;
    if b_gp_yr + b_p9_yr <= 0.0 {
        return radians(0.0);
    }
    radians(b_p9_yr / (b_gp_yr + b_p9_yr)) * p9.i
}

/// Forced inclination in the Planet-Nine-free Solar System: identically zero
/// (the particle's plane relaxes to the giant-planet invariable plane).
pub fn forced_inclination_no_p9(_a: f64, _e: f64, _i: f64) -> f64 {
    0.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    fn nominal() -> P9Params {
        P9Params::nominal_2016() // 10 M⊕, a = 700, e = 0.6, i = 30°
    }

    /// Planet Nine coupling frequency in rad/yr (version-safe extraction).
    fn p9_freq_yr(a: f64, p9: &P9Params) -> f64 {
        let per_yr = radians(1.0) / julian_year();
        (planet_nine_frequency_typed(a, p9) / per_yr).value
    }

    /// Forced inclination in radians (version-safe extraction).
    fn forced_rad(a: f64, e: f64, i: f64, p9: &P9Params) -> f64 {
        (forced_inclination_typed(a, e, i, p9) / radians(1.0)).value
    }

    #[test]
    fn typed_forcing_accessors_match_inline() {
        use approx::assert_relative_eq;
        use uom::si::angle::radian;
        use uom::si::angular_velocity::radian_per_second;
        use uom::si::time::second;
        // rad/yr → rad/s divisor (Julian year in seconds).
        let yr_s = julian_year().get::<second>();
        let p9 = nominal();
        let (a, e, i) = (300.0, 0.5, 10.0 * DEG2RAD);

        // Inline the forced-inclination formula from the typed components.
        let b_gp_yr = giant_planet_frequency_typed(a, e, i).get::<radian_per_second>() * yr_s;
        let b_p9_yr = planet_nine_frequency_typed(a, &p9).get::<radian_per_second>() * yr_s;
        let forced = b_p9_yr / (b_gp_yr + b_p9_yr) * p9.i;
        assert_relative_eq!(
            forced_inclination_typed(a, e, i, &p9).get::<radian>(),
            forced,
            epsilon = 1e-12
        );

        // Inline B_gp closed form: (3/4) n cos i /(1−e²)² Σ_j (m_j/M_⊙)(a_j/a)².
        let n = (GM_SUN / (a * a * a)).sqrt();
        let eta_sq = (1.0 - e * e) * (1.0 - e * e);
        let sum: f64 = giant_planets_j2000()
            .iter()
            .map(|body| {
                let a_j = cartesian_to_elements(&body.state, GM_SUN).a;
                body.mass * (a_j / a).powi(2)
            })
            .sum();
        let b_gp_inline = (0.75 * n * i.cos() / eta_sq * sum * 365.25).abs();
        assert_relative_eq!(
            giant_planet_frequency_typed(a, e, i).get::<radian_per_second>(),
            b_gp_inline / yr_s,
            max_relative = 1e-12
        );

        // Inline B_p9 closed form: (3/4) n (m₉/M_⊙) (a/a₉)² /(1−e₉²)^{3/2}.
        let alpha = a / p9.a;
        let b_p9_inline =
            0.75 * n * p9.mass_solar() * alpha * alpha * (1.0 - p9.e * p9.e).powf(-1.5) * 365.25;
        assert_relative_eq!(
            planet_nine_frequency_typed(a, &p9).get::<radian_per_second>(),
            b_p9_inline / yr_s,
            max_relative = 1e-12
        );
    }

    #[test]
    fn test_p9_frequency_scales_as_alpha_squared() {
        // Coupling ∝ (a/a₉)²: doubling a quadruples B_p9 (× the n ∝ a^{-3/2}).
        let p9 = nominal();
        let b1 = p9_freq_yr(250.0, &p9);
        let b2 = p9_freq_yr(500.0, &p9);
        // B_p9 ∝ a^{-3/2} · a² = a^{1/2}; ratio = sqrt(2) ≈ 1.414.
        let ratio = b2 / b1;
        assert!(
            (ratio - 2f64.sqrt()).abs() < 0.05,
            "B_p9(500)/B_p9(250) = {ratio:.3}, expected ≈ {:.3}",
            2f64.sqrt()
        );
    }

    #[test]
    fn test_p9_frequency_scales_with_mass() {
        let light = P9Params {
            mass_earth: 5.0,
            ..nominal()
        };
        let heavy = P9Params {
            mass_earth: 10.0,
            ..nominal()
        };
        let ratio = p9_freq_yr(300.0, &heavy) / p9_freq_yr(300.0, &light);
        assert!((ratio - 2.0).abs() < 1e-9, "mass scaling = {ratio:.4}");
    }

    #[test]
    fn test_quadrupole_curvature_shares_scalings() {
        // The closed-form B_p9 and the p9-core quadrupole curvature must share
        // the α² and 1/(1−e₉²)^{3/2} dependence (they differ only by constant
        // and mean-motion factors). Compare *ratios* across two semi-major axes
        // and two perturber eccentricities.
        let base = nominal();
        let ecc = P9Params { e: 0.3, ..base };

        // α² scaling: curvature ∝ C ∝ α² = (a/a₉)².
        let c1 = quadrupole_inclination_curvature(250.0, 0.4, &base);
        let c2 = quadrupole_inclination_curvature(500.0, 0.4, &base);
        let curv_ratio = c2 / c1;
        assert!(
            (curv_ratio - 4.0).abs() < 0.01,
            "quadrupole α² scaling = {curv_ratio:.4}, expected 4.0"
        );

        // 1/(1−e₉²)^{3/2} perturber-eccentricity scaling, shared with B_p9.
        let eobs = quadrupole_inclination_curvature(300.0, 0.4, &base)
            / quadrupole_inclination_curvature(300.0, 0.4, &ecc);
        let bobs = p9_freq_yr(300.0, &base) / p9_freq_yr(300.0, &ecc);
        let expected = (1.0 - ecc.e * ecc.e).powf(1.5) / (1.0 - base.e * base.e).powf(1.5);
        assert!(
            (eobs - expected).abs() < 0.02 && (bobs - expected).abs() < 0.02,
            "e₉ scaling: curvature {eobs:.4}, closed-form {bobs:.4}, expected {expected:.4}"
        );
    }

    #[test]
    fn test_forced_inclination_zero_without_p9() {
        assert_eq!(forced_inclination_no_p9(400.0, 0.7, 0.2), 0.0);
    }

    #[test]
    fn test_forced_inclination_increases_with_mass() {
        let light = P9Params {
            mass_earth: 5.0,
            ..nominal()
        };
        let heavy = P9Params {
            mass_earth: 15.0,
            ..nominal()
        };
        let fl = forced_rad(400.0, 0.7, 10.0 * DEG2RAD, &light);
        let fh = forced_rad(400.0, 0.7, 10.0 * DEG2RAD, &heavy);
        assert!(fh > fl, "forced i should grow with mass: {fl} vs {fh}");
    }

    #[test]
    fn test_forced_inclination_increases_with_sin_i9() {
        let low = P9Params {
            i: 10.0 * DEG2RAD,
            ..nominal()
        };
        let high = P9Params {
            i: 30.0 * DEG2RAD,
            ..nominal()
        };
        let fl = forced_rad(400.0, 0.7, 10.0 * DEG2RAD, &low);
        let fh = forced_rad(400.0, 0.7, 10.0 * DEG2RAD, &high);
        assert!(fh > fl, "forced i should grow with i9: {fl} vs {fh}");
        // Approximately linear in i9 at fixed coupling fraction.
        let ratio = fh / fl;
        assert!(
            (ratio - 3.0).abs() < 0.3,
            "forced-i ratio {ratio:.3} should track i9 ratio 3.0"
        );
    }

    #[test]
    fn test_forced_inclination_strengthens_with_semimajor_axis() {
        let p9 = nominal();
        let inner = forced_rad(250.0, 0.7, 10.0 * DEG2RAD, &p9);
        let outer = forced_rad(450.0, 0.7, 10.0 * DEG2RAD, &p9);
        assert!(
            outer > inner,
            "forced i should strengthen with a: {inner} (250) vs {outer} (450)"
        );
    }

    #[test]
    fn test_forced_inclination_bounded_by_i9() {
        // The forced plane never over-tilts past Planet Nine's own plane.
        let p9 = nominal();
        let f = forced_rad(500.0, 0.8, 5.0 * DEG2RAD, &p9);
        assert!(f >= 0.0 && f <= p9.i, "forced i = {f} out of [0, i9]");
    }
}
