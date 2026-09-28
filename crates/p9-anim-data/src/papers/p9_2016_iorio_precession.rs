//! Film export for `p9-2016-iorio-precession`: the numbers its scene and ledger entry draw.
//!
//! Iorio (arXiv:1512.05288) holds the perturber fixed at its present position
//! and asks where along its orbit the induced secular precession of Saturn
//! stays inside the ephemeris bounds. The reproduction crate carries the
//! orbit-averaged quadrupole rate, so the true-anomaly dependence exported here
//! is that same closed form evaluated for a circular ring at the planet's
//! present heliocentric distance r(f): the direction-averaged tide of a body
//! parked at r(f).

use p9_2016_iorio_precession::{
    Perturber, bound_for, critical_distance_au, planet_bounds,
    planet_perihelion_precession_arcsec_per_cy,
};
use p9_core::data::ephemeris_constraint::brown_batygin_orbit;
use p9_core::types::{helio_distance_at_true_anomaly, true_to_mean_anomaly};
use serde_json::{Value, json};
use std::f64::consts::TAU;

/// Published admissible true-anomaly interval (degrees), arXiv:1512.05288 abstract.
const PUBLISHED_ALLOWED_DEG: (f64, f64) = (130.0, 240.0);

const MAS_PER_ARCSEC: f64 = 1.0e3;

pub fn export() -> Value {
    let orbit = brown_batygin_orbit();
    let averaged = Perturber {
        name: "Planet Nine",
        mass_earth: orbit.mass_earth,
        a_au: orbit.a,
        e: orbit.e,
    };
    let parked_at = |r_au: f64| Perturber {
        name: "Planet Nine",
        mass_earth: orbit.mass_earth,
        a_au: r_au,
        e: 0.0,
    };
    let perihelion = orbit.a * (1.0 - orbit.e);
    let aphelion = orbit.a * (1.0 + orbit.e);

    // Every planet with a bound: what the perturber does to it at perihelion,
    // orbit-averaged, and at aphelion, against what the ephemerides allow.
    let planets: Vec<Value> = planet_bounds()
        .iter()
        .map(|pl| {
            let rate =
                |p: &Perturber| MAS_PER_ARCSEC * planet_perihelion_precession_arcsec_per_cy(pl, p);
            json!({
                "name": pl.name,
                "a_au": pl.a_au,
                "bound_mas_cy": MAS_PER_ARCSEC * pl.bound_arcsec_per_cy,
                "rate_averaged_mas_cy": rate(&averaged),
                "rate_perihelion_mas_cy": rate(&parked_at(perihelion)),
                "rate_aphelion_mas_cy": rate(&parked_at(aphelion)),
            })
        })
        .collect();

    // Saturn carries the constraint: its rate against true anomaly.
    let saturn = bound_for("Saturn").expect("Saturn is in the bounds table");
    let f_deg: Vec<f64> = (0..=180).map(|k| 2.0 * k as f64).collect();
    let r_au: Vec<f64> = f_deg
        .iter()
        .map(|f| helio_distance_at_true_anomaly(&orbit, f.to_radians()))
        .collect();
    let rate_mas: Vec<f64> = r_au
        .iter()
        .map(|&r| {
            MAS_PER_ARCSEC * planet_perihelion_precession_arcsec_per_cy(&saturn, &parked_at(r))
        })
        .collect();

    // Distance at which the rate equals Saturn's bound, and the true anomalies
    // where the orbit crosses it.
    let r_crit = critical_distance_au(&saturn, &parked_at(orbit.a));
    let cos_f = (orbit.a * (1.0 - orbit.e * orbit.e) / r_crit - 1.0) / orbit.e;
    let f_lo = cos_f.clamp(-1.0, 1.0).acos().to_degrees();
    let f_hi = 360.0 - f_lo;
    // Share of the orbital period spent on the allowed arc (Kepler's second law).
    let m_of = |f_deg: f64| true_to_mean_anomaly(orbit.e, f_deg.to_radians());
    let allowed_time_fraction = (m_of(f_hi) - m_of(f_lo)) / TAU;

    json!({
        "orbit": {
            "mass_earth": orbit.mass_earth,
            "a_au": orbit.a,
            "e": orbit.e,
            "perihelion_au": perihelion,
            "aphelion_au": aphelion,
        },
        "planets": planets,
        "saturn": {
            "bound_mas_cy": MAS_PER_ARCSEC * saturn.bound_arcsec_per_cy,
            "f_deg": f_deg,
            "r_au": r_au,
            "rate_mas_cy": rate_mas,
        },
        "critical_distance_au": r_crit,
        "allowed_from_deg": f_lo,
        "allowed_to_deg": f_hi,
        "allowed_width_deg": f_hi - f_lo,
        "allowed_time_fraction": allowed_time_fraction,
        "published": {
            "allowed_from_deg": PUBLISHED_ALLOWED_DEG.0,
            "allowed_to_deg": PUBLISHED_ALLOWED_DEG.1,
        },
    })
}
