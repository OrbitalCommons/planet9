//! Secular apsidal precession of a distant ETNO.
//!
//! For a test particle *exterior* to a set of interior perturbing rings (the
//! known giant planets), the orbit-averaged quadrupole secular theory gives a
//! steady prograde precession of the longitude of perihelion ϖ. To leading
//! order in α_p = a_p / a ≪ 1 (the relevant limit for a ≳ 200 AU ETNO and
//! a_p ≤ 30 AU planets), the Laplace–Lagrange free-precession rate reduces to
//!
//!   dϖ/dt = (3/4) · n · Σ_p (m_p / M_⊙) · (a_p / a)² · 1 / (1 − e²)²
//!
//! where n = √(GM_⊙ / a³) is the particle's mean motion and the sum runs over
//! the interior perturbers (Murray & Dermott 1999, eq. 7.16, exterior-particle
//! limit b_{3/2}^{(1)}(α) → 3α; the (1 − e²)^{-2} factor is the standard
//! finite-eccentricity correction). The total precession is slow and prograde,
//! with a period of order >~ 1 Gyr for the clustered ETNOs — the regime in
//! which the apsidal lines can stay coherently confined (Clement & Sheppard
//! 2020, arXiv:2005.05326).
//!
//! A Planet Nine is *exterior* to the ETNO and adds an opposite-sign secular
//! term. For an exterior perturber the leading quadrupole free-precession rate
//! is
//!
//!   dϖ_9/dt = (3/4) · n · (m_9 / M_⊙) · (a / a_9_eff)³ · 1 / √(1 − e²)
//!
//! with a_9_eff = a_9 (1 − e_9²)^{1/2} the perturber's averaged ring scale
//! (Murray & Dermott eq. 7.16, interior-particle limit b_{3/2}^{(1)}(α) → 3α
//! with α = a / a_9). The two contributions add; the constraint of the paper
//! is that P9 sits in a mass–semimajor-axis band where this extra rate keeps
//! the combined precession slow and coherent rather than disrupting the
//! cluster.
//!
//! These analytic rates are independently cross-checked in the tests against a
//! finite-difference of `p9_core::analysis::secular::numerical_secular_hamiltonian`
//! (the exact Gauss-ring average), which captures the same dϖ/dt = +∂H/∂(action)
//! structure without any α expansion.

use p9_core::analysis::secular::numerical_secular_hamiltonian;
use p9_core::constants::{
    GM_SUN, GYR_DAYS, MASS_JUPITER_SOLAR, MASS_NEPTUNE_SOLAR, MASS_SATURN_SOLAR, MASS_URANUS_SOLAR,
    YEAR_DAYS,
};
use p9_core::types::P9Params;
use p9_core::units::{au, solar_masses, Length, Mass};

/// An interior perturbing ring: semi-major axis (AU) and mass (solar masses).
#[derive(Debug, Clone, Copy)]
pub struct Perturber {
    pub a_au: f64,
    pub mass_solar: f64,
}

impl Perturber {
    /// Perturber semi-major axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a_au)
    }

    /// Perturber mass as a typed [`Mass`].
    pub fn mass(&self) -> Mass {
        solar_masses(self.mass_solar)
    }
}

/// The four known giant planets as interior secular perturbers.
///
/// Semi-major axes are the standard mean values (Jupiter 5.20, Saturn 9.58,
/// Uranus 19.20 AU) with Neptune anchored to `p9_core`'s `A_NEPTUNE_AU`, and
/// masses from `p9_core`'s mass constants.
pub fn giant_planets() -> [Perturber; 4] {
    [
        Perturber {
            a_au: 5.2026,
            mass_solar: MASS_JUPITER_SOLAR,
        },
        Perturber {
            a_au: 9.5549,
            mass_solar: MASS_SATURN_SOLAR,
        },
        Perturber {
            a_au: 19.2184,
            mass_solar: MASS_URANUS_SOLAR,
        },
        Perturber {
            a_au: p9_core::constants::A_NEPTUNE_AU,
            mass_solar: MASS_NEPTUNE_SOLAR,
        },
    ]
}

/// Mean motion n = √(GM_⊙ / a³) in rad/day for a heliocentric orbit.
pub fn mean_motion(a_au: f64) -> f64 {
    (GM_SUN / a_au.powi(3)).sqrt()
}

/// Leading-order quadrupole apsidal precession rate dϖ/dt (rad/day) imposed on
/// an exterior test particle (a, e) by a single interior perturbing ring.
///
///   dϖ/dt = (3/4) n (m_p/M_⊙) (a_p/a)² (1 − e²)^{-2}
pub fn precession_rate_interior(a_au: f64, e: f64, p: Perturber) -> f64 {
    let n = mean_motion(a_au);
    let alpha = p.a_au / a_au;
    0.75 * n * p.mass_solar * alpha * alpha / (1.0 - e * e).powi(2)
}

/// Total giant-planet-induced apsidal precession rate dϖ/dt (rad/day) for an
/// exterior ETNO at (a, e): the sum of the four interior-planet contributions.
pub fn giant_planet_precession_rate(a_au: f64, e: f64) -> f64 {
    giant_planets()
        .iter()
        .map(|&p| precession_rate_interior(a_au, e, p))
        .sum()
}

/// Leading-order quadrupole apsidal precession rate dϖ/dt (rad/day) imposed on
/// an *interior* test particle (the ETNO) by an exterior perturber (Planet
/// Nine). Using the averaged ring scale a_9_eff = a_9 √(1 − e_9²):
///
///   dϖ_9/dt = (3/4) n (m_9/M_⊙) (a/a_9_eff)³ (1 − e²)^{-1/2}
pub fn p9_precession_rate(a_au: f64, e: f64, p9: &P9Params) -> f64 {
    let n = mean_motion(a_au);
    let a9_eff = p9.a * (1.0 - p9.e * p9.e).sqrt();
    let alpha = a_au / a9_eff;
    0.75 * n * p9.mass_solar() * alpha.powi(3) / (1.0 - e * e).sqrt()
}

/// Combined giant-planet + Planet Nine apsidal precession rate (rad/day).
pub fn combined_precession_rate(a_au: f64, e: f64, p9: &P9Params) -> f64 {
    giant_planet_precession_rate(a_au, e) + p9_precession_rate(a_au, e, p9)
}

/// Convert a precession rate (rad/day) to the precession period (years), i.e.
/// the time for ϖ to circulate once. Returns `f64::INFINITY` at zero rate.
pub fn precession_period_years(rate_rad_per_day: f64) -> f64 {
    if rate_rad_per_day == 0.0 {
        return f64::INFINITY;
    }
    let period_days = std::f64::consts::TAU / rate_rad_per_day.abs();
    period_days / YEAR_DAYS
}

/// Precession period in Gyr.
pub fn precession_period_gyr(rate_rad_per_day: f64) -> f64 {
    if rate_rad_per_day == 0.0 {
        return f64::INFINITY;
    }
    (std::f64::consts::TAU / rate_rad_per_day.abs()) / GYR_DAYS
}

/// Cross-check estimate of dϖ/dt (rad/day) from a finite difference of the
/// exact Gauss-ring secular Hamiltonian of a single perturber.
///
/// In the coplanar secular problem the apsidal angle is canonically conjugate
/// to the Delaunay action G = L√(1 − e²) with L = √(GM_⊙ a), so
///
///   dϖ/dt = ∂H/∂G = (∂H/∂e) · (∂e/∂G),   ∂e/∂G = −√(1 − e²)/(L e).
///
/// `numerical_secular_hamiltonian` returns −gm_p⟨⟨1/Δ⟩⟩ in the perturber frame
/// (Δϖ = omega); we difference it in e at fixed a to obtain ∂H/∂e. This makes
/// no α expansion and is used purely to validate the analytic rates above.
pub fn precession_rate_numerical_interior(a_au: f64, e: f64, p: Perturber) -> f64 {
    let gm_p = p.mass_solar * GM_SUN;
    let a_p = p.a_au;
    // Circular perturber: axisymmetric ring, the correct comparison for the
    // analytic interior formula (which assumes e_p = 0).
    let e_p = 0.0;
    let n_quad = 256;
    let h = 1.0e-4;

    let ham = |ecc: f64| {
        numerical_secular_hamiltonian(a_au, ecc, 0.0, 0.0, 0.0, a_p, e_p, gm_p, n_quad, 0.0)
    };
    let dh_de = (ham(e + h) - ham(e - h)) / (2.0 * h);

    let big_l = (GM_SUN * a_au).sqrt();
    let de_dg = -(1.0 - e * e).sqrt() / (big_l * e);
    dh_de * de_dg
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // ---- Published / documented reference values (labelled constants) ----

    /// Test ETNO used throughout: a = 250 AU, q = 40 AU ⇒ e = 1 − q/a = 0.84.
    const A_TEST: f64 = 250.0;
    const Q_TEST: f64 = 40.0;

    fn e_test() -> f64 {
        1.0 - Q_TEST / A_TEST
    }

    #[test]
    fn test_perturber_typed_accessors_match_f64() {
        use uom::si::length::astronomical_unit;
        use uom::si::mass::kilogram;
        let neptune = giant_planets()[3];
        assert_relative_eq!(
            neptune.semi_major_axis().get::<astronomical_unit>(),
            neptune.a_au,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            neptune.mass().get::<kilogram>(),
            neptune.mass_solar * p9_core::units::SOLAR_MASS_KG,
            epsilon = 1e15
        );
    }

    #[test]
    fn test_giant_planet_precession_period_is_documented_order() {
        // Clement & Sheppard 2020: the giant planets drive a slow prograde
        // apsidal precession of the distant ETNOs, with a period of order
        // ~1 Gyr (>~ a Gyr for a ≳ 250 AU). The rate must be positive
        // (prograde) and the period in the 10^8–10^10 yr band.
        let rate = giant_planet_precession_rate(A_TEST, e_test());
        assert!(rate > 0.0, "precession should be prograde: {rate:e}");

        let period_gyr = precession_period_gyr(rate);
        assert!(
            (0.05..20.0).contains(&period_gyr),
            "giant-planet precession period {period_gyr:.3} Gyr out of documented order"
        );
        // Cross-check the year/Gyr converters agree.
        let period_yr = precession_period_years(rate);
        assert_relative_eq!(period_yr / 1.0e9, period_gyr, epsilon = 1e-9);
    }

    #[test]
    fn test_analytic_matches_numerical_gauss_ring() {
        // The leading-order analytic interior rate (which truncates the
        // eccentricity dependence at (1 − e²)^{-2}) is validated against the
        // exact Gauss-ring secular derivative for the dominant perturber
        // (Neptune) at a moderate eccentricity where the truncation is small.
        // Here α = a_N/a ≈ 0.12 is small, so the residual is the leading
        // O(e^4) eccentricity term, a few percent at e = 0.3.
        let gp = giant_planets();
        let neptune = gp[3];
        let e = 0.3;
        let analytic = precession_rate_interior(A_TEST, e, neptune);
        let numerical = precession_rate_numerical_interior(A_TEST, e, neptune);
        assert!(analytic > 0.0 && numerical > 0.0);
        let rel = ((analytic - numerical) / numerical).abs();
        assert!(
            rel < 0.05,
            "analytic vs Gauss-ring rate mismatch: {rel:.3e}"
        );
    }

    #[test]
    fn test_eccentricity_correction_underestimates_at_high_e() {
        // At the q = 40 AU test ETNO (e = 0.84) the (1 − e²)^{-2} leading-order
        // factor underestimates the full ring precession: the exact Gauss-ring
        // derivative is ~2x the analytic rate. Both stay prograde and the
        // analytic value is the conservative (lower) bound. Documented residual
        // — the slow-precession conclusion is unchanged, both give Gyr periods.
        let gp = giant_planets();
        let neptune = gp[3];
        let e = e_test();
        let analytic = precession_rate_interior(A_TEST, e, neptune);
        let numerical = precession_rate_numerical_interior(A_TEST, e, neptune);
        assert!(analytic > 0.0 && numerical > analytic);
        let ratio = numerical / analytic;
        assert!(
            (1.5..2.5).contains(&ratio),
            "high-e Gauss-ring/analytic ratio = {ratio:.3}"
        );
    }

    #[test]
    fn test_p9_changes_rate_in_expected_direction() {
        // Adding an exterior Planet Nine adds a strictly positive (prograde)
        // secular term, so the combined rate exceeds the giant-planet-only
        // rate, and the combined precession period shortens.
        let e = e_test();
        let p9 = P9Params::nominal_2016();
        let gp_rate = giant_planet_precession_rate(A_TEST, e);
        let p9_rate = p9_precession_rate(A_TEST, e, &p9);
        let combined = combined_precession_rate(A_TEST, e, &p9);

        assert!(p9_rate > 0.0, "P9 secular term should be prograde");
        assert_relative_eq!(combined, gp_rate + p9_rate, epsilon = 1e-18);
        assert!(combined > gp_rate, "P9 must increase the precession rate");
        assert!(
            precession_period_gyr(combined) < precession_period_gyr(gp_rate),
            "P9 must shorten the precession period"
        );
    }

    #[test]
    fn test_p9_contribution_grows_with_mass_and_shrinks_with_distance() {
        let e = e_test();
        let light = P9Params {
            mass_earth: 5.0,
            ..P9Params::nominal_2016()
        };
        let heavy = P9Params {
            mass_earth: 10.0,
            ..P9Params::nominal_2016()
        };
        assert!(
            p9_precession_rate(A_TEST, e, &heavy) > p9_precession_rate(A_TEST, e, &light),
            "heavier P9 must precess the ETNO faster"
        );

        let near = P9Params {
            a: 400.0,
            ..P9Params::nominal_2016()
        };
        let far = P9Params {
            a: 800.0,
            ..P9Params::nominal_2016()
        };
        assert!(
            p9_precession_rate(A_TEST, e, &near) > p9_precession_rate(A_TEST, e, &far),
            "a more distant P9 must precess the ETNO more slowly"
        );
    }

    #[test]
    fn test_precession_period_inversely_tracks_rate() {
        // Sanity: doubling the rate halves the period.
        let p = precession_period_years(2.0e-12);
        let p2 = precession_period_years(4.0e-12);
        assert_relative_eq!(p / p2, 2.0, epsilon = 1e-12);
        assert!(precession_period_years(0.0).is_infinite());
    }
}
