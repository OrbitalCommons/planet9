//! Film export for `p9-2026-iorio-precession`: the numbers its scene and ledger entry draw.
//!
//! The crate carries the orbit-averaged quadrupole precession of Saturn's
//! perihelion. Where the paper asks about a perturber's present position, the
//! same closed form is evaluated for a circular ring at the heliocentric
//! distance in question: the direction-averaged tide of a body parked there.
//! The candidate planets are the ones the paper tests (arXiv:2602.00802,
//! Sections 1 and 4).

use p9_2026_iorio_precession::{
    Perturber, SATURN_PRECESSION_BOUND_ARCSEC_PER_CY, critical_distance_au,
    saturn_perihelion_precession_arcsec_per_cy,
};
use p9_core::constants::{EARTH_MASS_SOLAR, GM_MERCURY, GM_SUN};
use serde_json::{Value, json};

const MAS_PER_ARCSEC: f64 = 1.0e3;

/// Formal uncertainty of Saturn's perihelion rate, rescaled tenfold as in the
/// paper's Section 4 (mas per century).
const PUBLISHED_BOUND_MAS_CY: f64 = 0.67;

/// Smallest heliocentric distance the paper leaves open to a 4.9 Earth-mass
/// Planet Nine, and smallest semi-major axis it leaves open to a Mercury-mass
/// Planet Y (AU).
const PUBLISHED_P9_LIGHT_MIN_DISTANCE_AU: f64 = 560.0;
const PUBLISHED_PLANET_Y_MIN_A_AU: f64 = 125.0;

/// A candidate planet as the paper specifies it, with the paper's verdict.
struct Candidate {
    name: &'static str,
    mass_earth: f64,
    a_au: f64,
    e: f64,
    /// Published admissible true-anomaly interval (degrees), if any survives.
    published_allowed_deg: Option<(f64, f64)>,
    /// Published minimum heliocentric distance (AU), where quoted.
    published_min_distance_au: Option<f64>,
    /// The paper's verdict in words.
    published_verdict: &'static str,
}

fn candidates(mercury_mass_earth: f64) -> [Candidate; 5] {
    [
        Candidate {
            name: "Planet Nine, 4.9 M⊕",
            mass_earth: 4.9,
            a_au: 520.0,
            e: 0.538,
            published_allowed_deg: Some((140.0, 220.0)),
            published_min_distance_au: Some(PUBLISHED_P9_LIGHT_MIN_DISTANCE_AU),
            published_verdict: "beyond 560 AU",
        },
        Candidate {
            name: "Planet Nine, 8.4 M⊕",
            mass_earth: 8.4,
            a_au: 520.0,
            e: 0.538,
            published_allowed_deg: Some((160.0, 200.0)),
            published_min_distance_au: Some(670.0),
            published_verdict: "beyond 670 AU",
        },
        Candidate {
            name: "Planet X, 4 M⊕",
            mass_earth: 4.0,
            a_au: 290.0,
            e: 0.28,
            published_allowed_deg: None,
            published_min_distance_au: None,
            published_verdict: "ruled out",
        },
        Candidate {
            name: "Planet Y, Earth mass",
            mass_earth: 1.0,
            a_au: 125.0,
            e: 0.2,
            published_allowed_deg: None,
            published_min_distance_au: None,
            published_verdict: "ruled out",
        },
        Candidate {
            name: "Planet Y, Mercury mass",
            mass_earth: mercury_mass_earth,
            a_au: 125.0,
            e: 0.2,
            published_allowed_deg: None,
            published_min_distance_au: None,
            published_verdict: "allowed from 125 AU",
        },
    ]
}

fn parked(mass_earth: f64, r_au: f64) -> Perturber {
    Perturber {
        name: "parked",
        mass_earth,
        a_au: r_au,
        e: 0.0,
    }
}

pub fn export() -> Value {
    let mercury_mass_earth = GM_MERCURY / GM_SUN / EARTH_MASS_SOLAR;
    let bound_mas = MAS_PER_ARCSEC * SATURN_PRECESSION_BOUND_ARCSEC_PER_CY;

    // The bound in the (distance, mass) plane: the distance inside which a
    // body of each mass drives Saturn's perihelion faster than the bound.
    let masses: Vec<f64> = (0..=60)
        .map(|k| 10f64.powf(-1.5 + 2.7 * k as f64 / 60.0))
        .collect();
    let boundary_au: Vec<f64> = masses
        .iter()
        .map(|&m| critical_distance_au(&parked(m, 100.0)))
        .collect();
    // The same boundary for the paper's bound: the rate falls as distance^-3.
    let published_scale = (bound_mas / PUBLISHED_BOUND_MAS_CY).cbrt();

    let cases: Vec<Value> = candidates(mercury_mass_earth)
        .iter()
        .map(|c| {
            let q = c.a_au * (1.0 - c.e);
            let big_q = c.a_au * (1.0 + c.e);
            let r_crit = critical_distance_au(&parked(c.mass_earth, c.a_au));
            let allowed_from = if r_crit <= q {
                Some(0.0)
            } else if r_crit >= big_q {
                None
            } else {
                let cos_f = (c.a_au * (1.0 - c.e * c.e) / r_crit - 1.0) / c.e;
                Some(cos_f.acos().to_degrees())
            };
            let f_deg: Vec<f64> = (0..=180).map(|k| 2.0 * k as f64).collect();
            let r_au: Vec<f64> = f_deg
                .iter()
                .map(|f| c.a_au * (1.0 - c.e * c.e) / (1.0 + c.e * f.to_radians().cos()))
                .collect();
            let rate: Vec<f64> = r_au
                .iter()
                .map(|&r| {
                    MAS_PER_ARCSEC
                        * saturn_perihelion_precession_arcsec_per_cy(&parked(c.mass_earth, r))
                })
                .collect();
            json!({
                "name": c.name,
                "mass_earth": c.mass_earth,
                "a_au": c.a_au,
                "e": c.e,
                "perihelion_au": q,
                "aphelion_au": big_q,
                "critical_distance_au": r_crit,
                "critical_distance_paper_bound_au": r_crit * published_scale,
                "allowed_from_deg": allowed_from,
                "allowed_to_deg": allowed_from.map(|f| 360.0 - f),
                "f_deg": f_deg,
                "r_au": r_au,
                "rate_mas_cy": rate,
                "published_allowed_deg": c.published_allowed_deg,
                "published_min_distance_au": c.published_min_distance_au,
                "published_verdict": c.published_verdict,
            })
        })
        .collect();

    json!({
        "bound_mas_cy": bound_mas,
        "boundary": {
            "mass_earth": masses,
            "distance_au": boundary_au,
            "published_bound_distance_au":
                boundary_au.iter().map(|r| r * published_scale).collect::<Vec<_>>(),
        },
        "cases": cases,
        "mercury_mass_earth": mercury_mass_earth,
        "p9_light_min_distance_au": critical_distance_au(&parked(4.9, 100.0)),
        "p9_heavy_min_distance_au": critical_distance_au(&parked(8.4, 100.0)),
        "planet_y_mercury_min_distance_au":
            critical_distance_au(&parked(mercury_mass_earth, 100.0)),
        "published": {
            "bound_mas_cy": PUBLISHED_BOUND_MAS_CY,
            "p9_light_min_distance_au": PUBLISHED_P9_LIGHT_MIN_DISTANCE_AU,
            "planet_y_mercury_min_a_au": PUBLISHED_PLANET_Y_MIN_A_AU,
        },
    })
}
