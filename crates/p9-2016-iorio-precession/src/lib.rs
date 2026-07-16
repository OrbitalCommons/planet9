//! Reproduction of Iorio (2016), "Constraints on Planet Nine from the
//! anomalous perihelion precession of the planets" (arXiv:1512.05288).
//!
//! # Headline
//!
//! A distant, massive Planet Nine raises a tidal quadrupole at the location of
//! the known planets. Orbit-averaging that quadrupole over the perturber's
//! orbit (it is far enough to be treated as a fixed external ring) leaves a
//! static, axisymmetric quadrupole potential whose secular effect on each
//! planet is a slow, prograde precession of its longitude of perihelion. The
//! rate scales as the tidal-field gradient,
//!
//! ```text
//!   d(varpi)/dt  ∝  GM_p9 / a_p9^3 ,
//! ```
//!
//! i.e. linearly in the perturber mass and as the inverse cube of its distance.
//! For an inner planet of semi-major axis `a` and mean motion `n` the closed
//! form (coplanar, doubly orbit-averaged quadrupole) is
//!
//! ```text
//!   d(varpi)/dt = (3/4) n (m_p9/M_sun) (a/a_p9)^3 * G(e, e_p9),
//!   G(e, e_p9) = sqrt(1 - e^2) / (1 - e_p9^2)^{3/2}.
//! ```
//!
//! Planetary ephemerides (EPM2011/INPOP-class solutions; Pitjeva & Pitjev 2013,
//! Fienga et al.) bound any *anomalous* (unmodelled) secular apsidal precession
//! of each planet at the milliarcsecond-per-century level or tighter. Comparing
//! the P9-induced rate to those bounds excludes the closer/more-massive corner
//! of the (mass, distance) plane.
//!
//! # Relationship to the 2016 paper and the 2026 sibling crate
//!
//! Iorio (2016) is the *earlier* analysis: it predates the
//! `p9-2026-iorio-precession` crate (Iorio 2026, a Saturn-only Cassini-ranging
//! constraint). This crate reuses the 2026 crate's secular precession machinery
//! --- the `Perturber` type, its `gm`/`mass_solar` accessors, and the exact
//! same closed-form scaling --- but applies it to (a) the 2016-era *per-planet*
//! anomalous-precession bounds and (b) *multiple* planets (Saturn plus the
//! terrestrial planets and Jupiter), which is the distinguishing content of the
//! 2016 paper. We never reimplement the quadrupole; `p9-core`'s
//! `coplanar_quadrupole` is used as the consistency cross-check (re-exported via
//! the 2026 crate's `quadrupole_prefactor_consistency`).
//!
//! # Honest approximations
//!
//! * Leading **quadrupole** term only in a/a_p9. With the planets at <= 9.6 AU
//!   and P9 at >= 100 AU, a/a_p9 <= 0.1, so octupole corrections are O(a/a_p9)
//!   <~ 10% --- ample for an order-of-magnitude exclusion.
//! * Perturber taken **coplanar** with each planet (mutual inclination = 0),
//!   the conservative (largest in-plane apsidal) choice. A general inclination
//!   multiplies G by a (1/2, 1)-bounded function of cos^2 i.
//! * The adopted per-planet bounds (`PlanetBound`) are *reference
//!   constants* drawn from the published EPM2011/INPOP supplementary-precession
//!   estimates, not computed quantities; tightening them only sharpens the
//!   exclusions. See `planet_bounds()` for provenance.

use serde::{Deserialize, Serialize};

use p9_core::analysis::secular::perturber_perihelion_precession;
use p9_core::constants::GM_SUN;
use p9_core::units::{arcseconds, au, days, julian_centuries, radians, AngularVelocity, Length};

// Reuse the 2026 sibling's perturber model and consistency check verbatim.
pub use p9_2026_iorio_precession::{
    quadrupole_prefactor_consistency, saturn_mean_motion_rad_per_day, Perturber,
};

/// A known planet whose perihelion precession we constrain, with the adopted
/// observational bound on its *anomalous* secular apsidal rate.
///
/// Semi-major axes and eccentricities are JPL mean Keplerian elements (J2000).
/// The bound is the published 1-sigma magnitude of the supplementary
/// (unmodelled) perihelion precession.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PlanetBound {
    /// Planet name.
    pub name: &'static str,
    /// Semi-major axis (AU).
    pub a_au: f64,
    /// Orbital eccentricity.
    pub e: f64,
    /// Adopted upper bound on the anomalous secular perihelion precession,
    /// arcseconds per century.
    pub bound_arcsec_per_cy: f64,
}

impl PlanetBound {
    /// Keplerian mean motion n = sqrt(GM_sun / a^3) in radians/day.
    pub fn mean_motion_rad_per_day(&self) -> f64 {
        (GM_SUN / self.a_au.powi(3)).sqrt()
    }

    // ---- typed (uom) boundary: dimension-checked views of the f64 fields ----

    /// Semi-major axis as a [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a_au)
    }

    /// Adopted anomalous-precession bound as a typed [`AngularVelocity`]
    /// (arcseconds per Julian century).
    pub fn bound(&self) -> AngularVelocity {
        (arcseconds(self.bound_arcsec_per_cy) / julian_centuries(1.0)).into()
    }

    /// Keplerian mean motion as a typed [`AngularVelocity`].
    pub fn mean_motion(&self) -> AngularVelocity {
        (radians(self.mean_motion_rad_per_day()) / days(1.0)).into()
    }
}

/// The planets and 2016-era anomalous-precession bounds from the EPM2011
/// supplementary precessions Δϖ̇ of Pitjeva & Pitjev (2013, MNRAS 432,
/// 3431, Table 4) that Iorio (2016) draws on: each bound is that table's
/// 1σ formal uncertainty (mas/century = 1e-3 arcsec/century), kept as
/// labelled reference constants:
///
/// | planet  | a (AU) | e      | bound (mas/cy) | P&P 2013 Δϖ̇        |
/// |---------|--------|--------|----------------|---------------------|
/// | Mercury | 0.3871 | 0.2056 | 2.0            | −2.0 ± 3.0          |
/// | Venus   | 0.7233 | 0.0068 | 1.5            |  2.6 ± 1.6          |
/// | Earth   | 1.0000 | 0.0167 | 0.2            |  0.19 ± 0.19        |
/// | Mars    | 1.5237 | 0.0934 | 0.04           | −0.020 ± 0.037      |
/// | Jupiter | 5.2029 | 0.0484 | 28.3           |  58.7 ± 28.3        |
/// | Saturn  | 9.5371 | 0.0539 | 0.47           | −0.32 ± 0.47        |
///
/// (A previous table adopted Jupiter 280 and Saturn 1.0 mas/cy while citing
/// the same source — neither matches it; Saturn's critical distances were
/// ~28% smaller than the cited bound gives.) Mars carries the tightest
/// absolute bound (~0.04 mas/cy from MRO/Mars Express ranging); Saturn the
/// tightest *outer*-planet bound (Cassini ranging), which is why the 2016
/// analysis leans on Saturn for the most distant perturbers. Jupiter's
/// bound is comparatively loose. These are reference constants, not
/// computed values.
pub fn planet_bounds() -> [PlanetBound; 6] {
    [
        PlanetBound {
            name: "Mercury",
            a_au: 0.387_098,
            e: 0.205_630,
            bound_arcsec_per_cy: 2.0e-3,
        },
        PlanetBound {
            name: "Venus",
            a_au: 0.723_332,
            e: 0.006_772,
            bound_arcsec_per_cy: 1.5e-3,
        },
        PlanetBound {
            name: "Earth",
            a_au: 1.000_000,
            e: 0.016_710,
            bound_arcsec_per_cy: 0.2e-3,
        },
        PlanetBound {
            name: "Mars",
            a_au: 1.523_679,
            e: 0.093_400,
            bound_arcsec_per_cy: 0.04e-3,
        },
        PlanetBound {
            name: "Jupiter",
            a_au: 5.202_887,
            e: 0.048_386,
            bound_arcsec_per_cy: 28.3e-3,
        },
        PlanetBound {
            name: "Saturn",
            a_au: 9.537_070,
            e: 0.053_862,
            bound_arcsec_per_cy: 0.47e-3,
        },
    ]
}

/// Look up the adopted bound for a named planet (case-sensitive).
pub fn bound_for(name: &str) -> Option<PlanetBound> {
    planet_bounds().into_iter().find(|p| p.name == name)
}

/// Secular precession rate of a planet's longitude of perihelion induced by a
/// coplanar external perturber, in **arcseconds per century**.
///
/// Closed form (see module docs), identical to the 2026 sibling's Saturn rate
/// but parameterised on the target planet:
/// ```text
///   d(varpi)/dt = (3/4) n (m_p9/M_sun) (a/a_p9)^3 * G(e, e_p9),
///   G = sqrt(1 - e^2) / (1 - e_p9^2)^{3/2}.
/// ```
pub fn planet_perihelion_precession_arcsec_per_cy(planet: &PlanetBound, p9: &Perturber) -> f64 {
    perturber_perihelion_precession(planet.a_au, planet.e, p9.mass_solar(), p9.a_au, p9.e)
}

/// Predicted perihelion precession of `planet` induced by `p9`, as a typed
/// [`AngularVelocity`] (the dimension-checked view of
/// [`planet_perihelion_precession_arcsec_per_cy`]).
pub fn planet_perihelion_precession(planet: &PlanetBound, p9: &Perturber) -> AngularVelocity {
    (arcseconds(planet_perihelion_precession_arcsec_per_cy(planet, p9)) / julian_centuries(1.0))
        .into()
}

/// True if this perturber's predicted precession of `planet` exceeds the
/// planet's adopted ephemeris bound (and is therefore *excluded* by that
/// planet's ranging).
pub fn is_excluded_by(planet: &PlanetBound, p9: &Perturber) -> bool {
    planet_perihelion_precession_arcsec_per_cy(planet, p9).abs() > planet.bound_arcsec_per_cy
}

/// True if *any* planet's bound excludes this perturber.
pub fn is_excluded(p9: &Perturber) -> bool {
    planet_bounds().iter().any(|pl| is_excluded_by(pl, p9))
}

/// The planet whose bound gives the strongest (smallest critical-distance,
/// i.e. most constraining) exclusion for this perturber, together with its
/// predicted rate. Returns the planet for which rate/bound is largest.
pub fn most_constraining<'a>(p9: &Perturber, bounds: &'a [PlanetBound]) -> (&'a PlanetBound, f64) {
    bounds
        .iter()
        .map(|pl| {
            let r = planet_perihelion_precession_arcsec_per_cy(pl, p9);
            (pl, r / pl.bound_arcsec_per_cy)
        })
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .map(|(pl, ratio)| (pl, ratio))
        .unwrap()
}

/// Solve for the perturber semi-major axis (AU) at which the induced precession
/// of `planet` exactly equals that planet's adopted bound, holding mass and
/// eccentricity fixed. Inside this radius the perturber is excluded by that
/// planet; outside it is allowed. Exact, since the rate is a pure a_p9^{-3}
/// power law: a_crit = a_p9 * (rate / bound)^{1/3}.
pub fn critical_distance_au(planet: &PlanetBound, p9: &Perturber) -> f64 {
    let rate = planet_perihelion_precession_arcsec_per_cy(planet, p9).abs();
    p9.a_au * (rate / planet.bound_arcsec_per_cy).cbrt()
}

/// Critical perturber distance as a typed [`Length`] (the dimension-checked
/// view of [`critical_distance_au`]).
pub fn critical_distance(planet: &PlanetBound, p9: &Perturber) -> Length {
    au(critical_distance_au(planet, p9))
}

/// The three P9 hypotheses, re-exported from the 2026 sibling
/// (`Perturber::planet_nine / planet_x / planet_y`).
pub fn hypotheses() -> [Perturber; 3] {
    [
        Perturber::planet_nine(),
        Perturber::planet_x(),
        Perturber::planet_y(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn saturn() -> PlanetBound {
        bound_for("Saturn").unwrap()
    }

    #[test]
    fn typed_accessors_match_f64_sources() {
        use uom::si::angular_velocity::radian_per_second;
        use uom::si::length::astronomical_unit;

        let s = saturn();
        let p9 = Perturber::planet_nine();
        // arcsec/century -> rad/s
        let arcsec = std::f64::consts::PI / (180.0 * 3600.0);
        let cy_s = 36_525.0 * 86_400.0;
        assert_relative_eq!(s.semi_major_axis().get::<astronomical_unit>(), s.a_au);
        assert_relative_eq!(
            s.bound().get::<radian_per_second>(),
            s.bound_arcsec_per_cy * arcsec / cy_s,
            epsilon = 1e-30
        );
        assert_relative_eq!(
            s.mean_motion().get::<radian_per_second>(),
            s.mean_motion_rad_per_day() / 86_400.0,
            epsilon = 1e-30
        );
        assert_relative_eq!(
            planet_perihelion_precession(&s, &p9).get::<radian_per_second>(),
            planet_perihelion_precession_arcsec_per_cy(&s, &p9) * arcsec / cy_s,
            epsilon = 1e-30
        );
        assert_relative_eq!(
            critical_distance(&s, &p9).get::<astronomical_unit>(),
            critical_distance_au(&s, &p9)
        );
    }

    #[test]
    fn test_saturn_mean_motion_matches_sibling() {
        // Our per-planet mean motion for Saturn must equal the 2026 sibling's
        // dedicated Saturn mean motion (same GM_sun, same formula). The 2026
        // crate rounds a_Sat to 9.5371 AU; we carry 9.53707 AU, a ~4.7e-6
        // relative difference that propagates as ~7e-6 into n ∝ a^{-3/2}.
        let ours = saturn().mean_motion_rad_per_day();
        let theirs = saturn_mean_motion_rad_per_day();
        assert_relative_eq!(ours, theirs, max_relative = 1e-5);
    }

    #[test]
    fn test_rate_scales_linearly_with_mass() {
        let pl = saturn();
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
        let r1 = planet_perihelion_precession_arcsec_per_cy(&pl, &p1);
        let r2 = planet_perihelion_precession_arcsec_per_cy(&pl, &p2);
        assert_relative_eq!(r2 / r1, 3.0, max_relative = 1e-12);
    }

    #[test]
    fn test_rate_scales_as_inverse_cube_of_distance() {
        // Two distances, same mass and e: ratio follows (a1/a2)^3 to ~1e-6.
        let pl = saturn();
        let near = Perturber {
            name: "near",
            mass_earth: 5.0,
            a_au: 200.0,
            e: 0.0,
        };
        let far = Perturber {
            name: "far",
            mass_earth: 5.0,
            a_au: 600.0,
            e: 0.0,
        };
        let r_near = planet_perihelion_precession_arcsec_per_cy(&pl, &near);
        let r_far = planet_perihelion_precession_arcsec_per_cy(&pl, &far);
        let cube = (far.a_au / near.a_au).powi(3); // = 27
        assert_relative_eq!(r_near / r_far, cube, max_relative = 1e-6);
    }

    #[test]
    fn test_eccentric_perturber_enhances_rate() {
        let pl = saturn();
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
        let r_c = planet_perihelion_precession_arcsec_per_cy(&pl, &circ);
        let r_e = planet_perihelion_precession_arcsec_per_cy(&pl, &ecc);
        let factor = (1.0 - 0.5_f64 * 0.5).powf(-1.5);
        assert_relative_eq!(r_e / r_c, factor, max_relative = 1e-12);
    }

    #[test]
    fn test_close_massive_saturn_excluded_distant_light_allowed() {
        // For a nominal close/massive P9 (10 M_earth at 100 AU) Saturn's
        // predicted precession exceeds the adopted 1 mas/cy bound (excluded);
        // a distant low-mass one (0.1 M_earth at 800 AU) does not.
        let pl = saturn();
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
            is_excluded_by(&pl, &close_massive),
            "close/massive should be excluded by Saturn ({:.3e} as/cy > {:.0e})",
            planet_perihelion_precession_arcsec_per_cy(&pl, &close_massive),
            pl.bound_arcsec_per_cy
        );
        assert!(
            !is_excluded_by(&pl, &distant_light),
            "distant/light should be allowed by Saturn ({:.3e} as/cy)",
            planet_perihelion_precession_arcsec_per_cy(&pl, &distant_light)
        );
    }

    #[test]
    fn test_precession_ordering_across_planets() {
        // Fixed perturber; the induced rate must increase with planet
        // semi-major axis (rate ∝ a^3 * n = a^3 * a^{-3/2} = a^{3/2}, since
        // GM_p9, a_p9 are common). So Saturn > Jupiter > Mars > ... > Mercury.
        let p9 = Perturber::planet_x();
        let bounds = planet_bounds();
        let rates: Vec<(f64, f64)> = bounds
            .iter()
            .map(|pl| (pl.a_au, planet_perihelion_precession_arcsec_per_cy(pl, &p9)))
            .collect();
        // The raw induced rate is monotone increasing in a (ordering with a).
        for w in rates.windows(2) {
            assert!(
                w[1].1 > w[0].1,
                "rate not monotone in a: a={} r={:.3e} then a={} r={:.3e}",
                w[0].0,
                w[0].1,
                w[1].0,
                w[1].1
            );
        }
        // Verify the a^{3/2} scaling explicitly between Mercury and Saturn.
        let merc = &bounds[0];
        let sat = &bounds[5];
        let r_merc = planet_perihelion_precession_arcsec_per_cy(merc, &p9);
        let r_sat = planet_perihelion_precession_arcsec_per_cy(sat, &p9);
        let pred = (sat.a_au / merc.a_au).powf(1.5)
            * ((1.0 - sat.e * sat.e).sqrt() / (1.0 - merc.e * merc.e).sqrt());
        assert_relative_eq!(r_sat / r_merc, pred, max_relative = 1e-9);
    }

    #[test]
    fn test_quadrupole_prefactor_matches_p9_core() {
        // The 2026 sibling's cross-check against p9_core::coplanar_quadrupole
        // (used by reference) must return unity for every hypothesis. This
        // confirms our closed-form prefactor inherits p9-core's tidal scaling.
        for p in hypotheses() {
            let ratio = quadrupole_prefactor_consistency(&p);
            assert_relative_eq!(ratio, 1.0, max_relative = 1e-12);
        }
    }

    #[test]
    fn test_critical_distance_self_consistent() {
        // Evaluating the rate at the critical distance reproduces the bound.
        let pl = saturn();
        let p = Perturber::planet_x();
        let a_crit = critical_distance_au(&pl, &p);
        let at_crit = Perturber { a_au: a_crit, ..p };
        let rate = planet_perihelion_precession_arcsec_per_cy(&pl, &at_crit);
        assert_relative_eq!(rate, pl.bound_arcsec_per_cy, max_relative = 1e-9);
    }

    #[test]
    fn test_most_constraining_planet_is_saturn_for_p9() {
        // With the CITED Pitjeva & Pitjev (2013) bounds, Saturn (Cassini
        // ranging, 0.47 mas/cy) edges out Mars for the nominal Planet Nine —
        // the raw a^{3/2}(a/a9)^3 rate advantage of the outermost planet
        // beats Mars's tighter absolute bound. This matches why the
        // literature leans on Saturn for distant perturbers. (Under the
        // previous mistranscribed Saturn bound of 1.0 mas/cy, Mars won.)
        let p9 = Perturber::planet_nine();
        let bounds = planet_bounds();
        let (best, ratio) = most_constraining(&p9, &bounds);
        assert!(ratio.is_finite() && ratio > 0.0);
        assert_eq!(
            best.name, "Saturn",
            "expected Saturn most constraining, got {} (ratio {:.3e})",
            best.name, ratio
        );
    }

    #[test]
    fn test_nominal_hypotheses_exclusion_status() {
        // Headline: report which nominal hypotheses are excluded by the planet
        // suite. Planet X (10 M_earth at 250 AU) must be excluded; record the
        // others as computed facts.
        let [p9, px, py] = hypotheses();
        assert!(is_excluded(&px), "nominal Planet X should be excluded");
        // These are recorded as computed facts, whichever way they fall.
        let _ = (is_excluded(&p9), is_excluded(&py));
        // Every rate must be finite and prograde for Saturn.
        let sat = saturn();
        for p in [p9, px, py] {
            let r = planet_perihelion_precession_arcsec_per_cy(&sat, &p);
            assert!(r > 0.0 && r.is_finite());
        }
    }

    #[test]
    fn test_bounds_table_well_formed() {
        let bounds = planet_bounds();
        assert_eq!(bounds.len(), 6);
        for b in &bounds {
            assert!(b.a_au > 0.0 && b.e >= 0.0 && b.e < 1.0);
            assert!(b.bound_arcsec_per_cy > 0.0);
        }
        // Semi-major axes strictly increasing Mercury -> Saturn.
        for w in bounds.windows(2) {
            assert!(w[1].a_au > w[0].a_au);
        }
    }
}
