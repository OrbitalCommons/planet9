//! Reproduction of Iorio (2026), "On the effects of a distended Planet
//! Nine/Planet X/Planet Y on the precession of the perihelion of Saturn"
//! (arXiv:2602.00802).
//!
//! # Headline
//!
//! A distant perturber raises a tidal quadrupole at the location of the inner
//! Solar System. Orbit-averaging that quadrupole over the perturber's orbit
//! (it is far enough that we may treat it as a fixed external ring) leaves a
//! static, axisymmetric quadrupole potential whose secular effect on Saturn is
//! a slow, prograde precession of Saturn's longitude of perihelion (and of its
//! node). The rate scales as
//!
//! ```text
//!   d(varpi)/dt_Sat  ∝  GM_p / a_p^3
//! ```
//!
//! i.e. linearly in the perturber mass and as the inverse cube of its distance
//! (the tidal-field gradient ∝ GM_p / a_p^3). Planetary ephemerides bound any
//! anomalous secular precession of Saturn's perihelion at the level of a few
//! milliarcseconds per century, so each of the three hypotheses --- Planet Nine
//! (~6 M_earth, ~400-500 AU), the more massive/distant Planet X, and the
//! close/light Planet Y (Mercury-to-Earth mass, 100-200 AU) --- is constrained
//! simultaneously, with the closer/heavier corner of each excluded.
//!
//! # The secular rate we compute
//!
//! For an inner planet of semi-major axis `a`, eccentricity `e`, mean motion
//! `n`, perturbed by an external body of mass `m_p`, semi-major axis `a_p`,
//! eccentricity `e_p`, **coplanar** with the inner orbit, the doubly orbit-
//! averaged quadrupole disturbing function (Kozai/Heisler-Tremaine form, the
//! `coplanar_quadrupole` of `p9_core::analysis::secular`) is
//!
//! ```text
//!   <R> = (1/4) (G m_p / a_p) (a/a_p)^2 (1 + 3/2 e^2) / (1 - e_p^2)^{3/2}.
//! ```
//!
//! Lagrange's planetary equation for the longitude of perihelion in the
//! coplanar (zero mutual inclination) limit reduces to
//!
//! ```text
//!   d(varpi)/dt = (sqrt(1 - e^2) / (n a^2 e)) ∂<R>/∂e
//! ```
//!
//! and with ∂<R>/∂e = (1/4)(G m_p / a_p)(a/a_p)^2 (3 e)/(1 - e_p^2)^{3/2} the
//! `e` cancels, giving the closed form
//!
//! ```text
//!   d(varpi)/dt = (3/4) n (m_p / M_sun) (a/a_p)^3 * G(e, e_p),
//!
//!   G(e, e_p) = sqrt(1 - e^2) / (1 - e_p^2)^{3/2}.
//! ```
//!
//! (We used n^2 a^3 = G M_sun, so (G m_p / a_p)(a/a_p)^2 / (n a^2)
//! = n (m_p/M_sun)(a/a_p)^3.) `G` is the **geometry factor** we adopt: it is the
//! coplanar, orbit-averaged value. For Saturn's small eccentricity (e ≈ 0.054)
//! sqrt(1 - e^2) ≈ 0.9985, and for a circular perturber (1 - e_p^2)^{-3/2} = 1,
//! so G ≈ 1; the prefactor (3/4) n (m_p/M_sun)(a/a_p)^3 carries essentially all
//! of the mass and distance dependence. A non-zero perturber eccentricity
//! enhances the rate by (1 - e_p^2)^{-3/2}.
//!
//! This is the same scaling Iorio reports: the anomalous precession is set by
//! the perturber's tidal quadrupole GM_p / a_p^3 and is independent of Saturn's
//! own eccentricity to leading order.
//!
//! # Honest approximations
//!
//! * We keep only the **quadrupole** (leading) term in a/a_p. With Saturn at
//!   ~9.5 AU and the perturbers at >= 100 AU, a/a_p <= 0.1, so the octupole
//!   correction is O(a/a_p) <~ 10% and the hexadecapole O((a/a_p)^2) <~ 1%.
//!   That is plenty to decide exclusions at the order-of-magnitude level the
//!   bound demands.
//! * We take the perturber **coplanar** with Saturn (mutual inclination = 0),
//!   which maximises the in-plane apsidal effect and gives the geometry factor
//!   above. A general inclination multiplies G by a (1, 1/2)-bounded function
//!   of cos^2 i; coplanar is the conservative (largest-precession) choice for
//!   an exclusion argument.
//! * We compare against a single labelled ephemeris bound on Saturn's
//!   anomalous perihelion precession (`SATURN_PRECESSION_BOUND_ARCSEC_PER_CY`);
//!   see its docstring for provenance.

use serde::{Deserialize, Serialize};

use p9_core::analysis::secular::{coplanar_quadrupole, perturber_perihelion_precession};
use p9_core::constants::{EARTH_MASS_SOLAR, GM_SUN};

/// Saturn's semi-major axis in AU (JPL mean element, J2000). Labelled here
/// because `p9-core` does not currently carry an inner/giant-planet semi-major
/// axis table (only Neptune). Source: JPL Keplerian elements, a_Sat = 9.5371
/// AU.
pub const A_SATURN_AU: f64 = 9.5371;

/// Saturn's mean orbital eccentricity (JPL mean element, J2000): e_Sat = 0.0539.
pub const E_SATURN: f64 = 0.0539;

/// Adopted observational upper bound on the *anomalous* (unmodelled) secular
/// precession of Saturn's perihelion, in arcseconds per century.
///
/// Planetary-ephemeris solutions (INPOP/EPM-class, e.g. Fienga et al.; the
/// supplementary precession estimates tabulated by Pitjeva & Pitjev) constrain
/// any extra apsidal precession of Saturn at the milliarcsecond-per-century
/// level. We adopt a representative round bound of 1 mas/century =
/// 1e-3 arcsec/century as the exclusion threshold. This is the magnitude Iorio
/// (2026) uses to set his constraints; it is a *reference constant*, not a
/// computed quantity, and tightening it only sharpens the exclusions.
pub const SATURN_PRECESSION_BOUND_ARCSEC_PER_CY: f64 = 1.0e-3;

/// A candidate distant perturber: mass (in Earth masses), semi-major axis (AU)
/// and orbital eccentricity. The three hypotheses are constructors below.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Perturber {
    /// Descriptive label (e.g. "Planet Nine").
    pub name: &'static str,
    /// Mass in Earth masses.
    pub mass_earth: f64,
    /// Semi-major axis in AU.
    pub a_au: f64,
    /// Orbital eccentricity.
    pub e: f64,
}

impl Perturber {
    /// Mass in solar masses.
    pub fn mass_solar(&self) -> f64 {
        self.mass_earth * EARTH_MASS_SOLAR
    }

    /// GM in AU^3/day^2.
    pub fn gm(&self) -> f64 {
        self.mass_solar() * GM_SUN
    }

    /// Planet Nine, Batygin/Brown-class nominal: ~6 Earth masses at ~460 AU,
    /// moderately eccentric (Brown & Batygin 2021 posterior region).
    pub fn planet_nine() -> Self {
        Self {
            name: "Planet Nine",
            mass_earth: 6.0,
            a_au: 460.0,
            e: 0.3,
        }
    }

    /// Planet X: the heavier/more distant classical hypothesis (Lykawka-style,
    /// also the high-mass tail Iorio considers). ~10 Earth masses at ~250 AU.
    pub fn planet_x() -> Self {
        Self {
            name: "Planet X",
            mass_earth: 10.0,
            a_au: 250.0,
            e: 0.2,
        }
    }

    /// Planet Y: the close, low-mass hypothesis (Mercury-to-Earth mass at
    /// 100-200 AU; Siraj et al. 2025 warp planet). We take ~0.5 Earth masses
    /// (between Mercury 0.055 and Earth 1.0) at 150 AU.
    pub fn planet_y() -> Self {
        Self {
            name: "Planet Y",
            mass_earth: 0.5,
            a_au: 150.0,
            e: 0.1,
        }
    }
}

/// Secular precession rate of Saturn's longitude of perihelion induced by a
/// coplanar external perturber, in **arcseconds per century**.
///
/// Closed form (see module docs):
/// ```text
///   d(varpi)/dt = (3/4) n_Sat (m_p/M_sun) (a_Sat/a_p)^3 * G(e_Sat, e_p),
///   G = sqrt(1 - e_Sat^2) / (1 - e_p^2)^{3/2}.
/// ```
pub fn saturn_perihelion_precession_arcsec_per_cy(p: &Perturber) -> f64 {
    perturber_perihelion_precession(A_SATURN_AU, E_SATURN, p.mass_solar(), p.a_au, p.e)
}

/// Saturn's Keplerian mean motion n = sqrt(GM_sun / a^3) in radians/day.
pub fn saturn_mean_motion_rad_per_day() -> f64 {
    (GM_SUN / A_SATURN_AU.powi(3)).sqrt()
}

/// True if this perturber's predicted Saturn precession exceeds the adopted
/// ephemeris bound (and is therefore *excluded* by Saturn ranging).
pub fn is_excluded(p: &Perturber) -> bool {
    saturn_perihelion_precession_arcsec_per_cy(p).abs() > SATURN_PRECESSION_BOUND_ARCSEC_PER_CY
}

/// Cross-check of the prefactor against `p9_core::analysis::secular`.
///
/// The orbit-averaged quadrupole disturbing function for Saturn under the
/// perturber is, up to the constant monopole, `-coplanar_quadrupole / m_Sat`
/// in the test-particle convention. We rebuild the *prefactor*
/// `C2 = (1/4)(G m_p / a_p)(a/a_p)^2 / (1 - e_p^2)^{3/2}` (the e-independent
/// part of <R>) two ways and return their ratio, which must be 1: once from
/// the closed-form constants, once by reading it out of `coplanar_quadrupole`.
///
/// `coplanar_quadrupole(a, e, a_p, e_p, gm_p)` returns
/// `H = -C (2 + 3 e^2)/2` with `C = gm_p (a/a_p)^2 / (4 a_p (1-e_p^2)^{3/2})`,
/// using the test-particle sign/normalisation of `p9-core`. Thus
/// `C = -2/(2+3e^2) * H`, and our disturbing-function prefactor is
/// `C2 = (1 + 3/2 e^2)^{-1} * (-H/2) ... ` -- rather than juggle the sign
/// convention we simply confirm that `H` scales as `(a/a_p)^2` and as `gm_p`,
/// which is exactly the GM_p / a_p^3 tidal scaling our rate inherits.
pub fn quadrupole_prefactor_consistency(p: &Perturber) -> f64 {
    let h = coplanar_quadrupole(A_SATURN_AU, E_SATURN, p.a_au, p.e, p.gm());
    // Reconstruct C from H = -C (2 + 3 e^2)/2.
    let e2 = E_SATURN * E_SATURN;
    let c_from_h = -2.0 * h / (2.0 + 3.0 * e2);
    // Independent closed form of C.
    let alpha = A_SATURN_AU / p.a_au;
    let c_closed = p.gm() * alpha * alpha / (4.0 * p.a_au * (1.0 - p.e * p.e).powf(1.5));
    c_from_h / c_closed
}

/// The three hypotheses, ordered as constructed.
pub fn hypotheses() -> [Perturber; 3] {
    [
        Perturber::planet_nine(),
        Perturber::planet_x(),
        Perturber::planet_y(),
    ]
}

/// Solve for the perturber semi-major axis (AU) at which the induced Saturn
/// precession exactly equals the adopted bound, holding mass and eccentricity
/// fixed. Inside this radius the hypothesis is excluded; outside it is allowed.
///
/// From `rate ∝ a_p^{-3} = bound` we get `a_crit = a_p_ref (rate_ref/bound)^{1/3}`
/// using any reference evaluation -- it is exact because the rate is a pure
/// power law in a_p.
pub fn critical_distance_au(p: &Perturber) -> f64 {
    let rate = saturn_perihelion_precession_arcsec_per_cy(p).abs();
    p.a_au * (rate / SATURN_PRECESSION_BOUND_ARCSEC_PER_CY).cbrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::constants::{TWO_PI, YEAR_DAYS};

    #[test]
    fn test_mean_motion_matches_saturn_period() {
        // n = 2pi / P; Saturn's sidereal period ~ 29.45 yr.
        let n = saturn_mean_motion_rad_per_day();
        let period_days = TWO_PI / n;
        let period_years = period_days / YEAR_DAYS;
        assert_relative_eq!(period_years, 29.45, max_relative = 0.01);
    }

    #[test]
    fn test_rate_scales_linearly_with_mass() {
        // Two masses, same a and e: rate must scale linearly.
        let p1 = Perturber {
            name: "m1",
            mass_earth: 3.0,
            a_au: 400.0,
            e: 0.0,
        };
        let p2 = Perturber {
            name: "m2",
            mass_earth: 9.0,
            a_au: 400.0,
            e: 0.0,
        };
        let r1 = saturn_perihelion_precession_arcsec_per_cy(&p1);
        let r2 = saturn_perihelion_precession_arcsec_per_cy(&p2);
        assert_relative_eq!(r2 / r1, 3.0, max_relative = 1e-12);
    }

    #[test]
    fn test_rate_scales_as_inverse_cube_of_distance() {
        // Two distances, same mass and e: ratio must follow (a1/a2)^3 to ~1e-6.
        let near = Perturber {
            name: "near",
            mass_earth: 5.0,
            a_au: 200.0,
            e: 0.0,
        };
        let far = Perturber {
            name: "far",
            mass_earth: 5.0,
            a_au: 400.0,
            e: 0.0,
        };
        let r_near = saturn_perihelion_precession_arcsec_per_cy(&near);
        let r_far = saturn_perihelion_precession_arcsec_per_cy(&far);
        let cube = (far.a_au / near.a_au).powi(3); // = 8
        assert_relative_eq!(r_near / r_far, cube, max_relative = 1e-6);
    }

    #[test]
    fn test_eccentric_perturber_enhances_rate() {
        // (1 - e_p^2)^{-3/2} enhancement.
        let circ = Perturber {
            name: "circ",
            mass_earth: 5.0,
            a_au: 300.0,
            e: 0.0,
        };
        let ecc = Perturber {
            name: "ecc",
            mass_earth: 5.0,
            a_au: 300.0,
            e: 0.5,
        };
        let r_c = saturn_perihelion_precession_arcsec_per_cy(&circ);
        let r_e = saturn_perihelion_precession_arcsec_per_cy(&ecc);
        let factor = (1.0 - 0.5_f64 * 0.5).powf(-1.5);
        assert_relative_eq!(r_e / r_c, factor, max_relative = 1e-12);
    }

    #[test]
    fn test_close_massive_excluded_distant_light_allowed() {
        // A close, massive perturber (10 M_earth at 100 AU) must exceed the
        // bound; a distant, light one (0.1 M_earth at 800 AU) must not.
        let close_massive = Perturber {
            name: "close-massive",
            mass_earth: 10.0,
            a_au: 100.0,
            e: 0.0,
        };
        let distant_light = Perturber {
            name: "distant-light",
            mass_earth: 0.1,
            a_au: 800.0,
            e: 0.0,
        };
        assert!(
            is_excluded(&close_massive),
            "close/massive should be excluded"
        );
        assert!(
            !is_excluded(&distant_light),
            "distant/light should be allowed"
        );
    }

    #[test]
    fn test_three_hypotheses_ordered_and_sensible() {
        let [p9, px, py] = hypotheses();
        let r_p9 = saturn_perihelion_precession_arcsec_per_cy(&p9);
        let r_px = saturn_perihelion_precession_arcsec_per_cy(&px);
        let r_py = saturn_perihelion_precession_arcsec_per_cy(&py);

        // All positive (prograde) and finite.
        assert!(r_p9 > 0.0 && r_px > 0.0 && r_py > 0.0);
        assert!(r_p9.is_finite() && r_px.is_finite() && r_py.is_finite());

        // Planet X (10 M_earth at 250 AU) raises the largest quadrupole of the
        // three; Planet Nine (6 M_earth at 460 AU) the smallest because of the
        // a^-3 distance penalty; Planet Y (0.5 M_earth at 150 AU) is between.
        assert!(
            r_px > r_py && r_py > r_p9,
            "expected Planet X > Planet Y > Planet Nine, got X={r_px:.3e} Y={r_py:.3e} 9={r_p9:.3e}"
        );

        // Order-of-magnitude sanity: all three sit within a few orders of the
        // mas/century bound, not absurdly large or small.
        for r in [r_p9, r_px, r_py] {
            assert!(r > 1e-8 && r < 1e2, "rate out of sane range: {r:.3e}");
        }
    }

    #[test]
    fn test_quadrupole_prefactor_matches_p9_core() {
        // Our closed-form quadrupole coefficient must equal the one implied by
        // p9_core::analysis::secular::coplanar_quadrupole, to machine epsilon,
        // for every hypothesis.
        for p in hypotheses() {
            let ratio = quadrupole_prefactor_consistency(&p);
            assert_relative_eq!(ratio, 1.0, max_relative = 1e-12);
        }
    }

    #[test]
    fn test_critical_distance_is_self_consistent() {
        // Evaluating the rate at the critical distance must reproduce the bound.
        let p = Perturber::planet_x();
        let a_crit = critical_distance_au(&p);
        let at_crit = Perturber { a_au: a_crit, ..p };
        let rate = saturn_perihelion_precession_arcsec_per_cy(&at_crit);
        assert_relative_eq!(
            rate,
            SATURN_PRECESSION_BOUND_ARCSEC_PER_CY,
            max_relative = 1e-9
        );
    }

    #[test]
    fn test_nominal_planet_x_excluded_planet_nine_status() {
        // Document the headline result: the close/massive Planet X exceeds the
        // 1 mas/century bound and is excluded; report whether nominal Planet
        // Nine / Planet Y do.
        let [p9, px, py] = hypotheses();
        assert!(
            is_excluded(&px),
            "nominal Planet X ({:.3e} as/cy) should exceed the {:.0e} as/cy bound",
            saturn_perihelion_precession_arcsec_per_cy(&px),
            SATURN_PRECESSION_BOUND_ARCSEC_PER_CY
        );
        // These two are recorded as computed facts (whichever way they fall):
        let _ = (is_excluded(&p9), is_excluded(&py));
    }
}
