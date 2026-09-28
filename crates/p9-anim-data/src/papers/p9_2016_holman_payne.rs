//! Film export for `p9-2016-holman-payne`: the numbers its scene and ledger entry draw.

use p9_2016_holman_payne::published::{
    CASSINI_RANGE_PRECISION_M, PREFERRED_DEC_DEG, PREFERRED_RA_DEG, PREFERRED_SKY_HALF_EXTENT_DEG,
};
use p9_2016_holman_payne::signal::range_residual;
use p9_2016_holman_payne::sky::favored_sky_position;
use p9_core::coords::sky::ecliptic_vec_to_equatorial_deg;
use p9_core::data::ephemeris_constraint::brown_batygin_orbit;
use p9_core::types::{P9Params, position_at_true_anomaly};
use p9_core::units::meters;
use serde_json::{Value, json};

/// True-anomaly samples scanned for the favoured position.
const N_SCAN: usize = 1440;

/// Tidal parameter the Cassini residuals favour, in units of a 10 Earth-mass
/// planet at 700 AU: Holman & Payne (2016) Eq. 11.
const PUBLISHED_TIDAL_RANGE: (f64, f64) = (2.3, 11.4);
const TIDAL_UNIT_MASS_EARTH: f64 = 10.0;
const TIDAL_UNIT_DISTANCE_AU: f64 = 700.0;

pub fn export() -> Value {
    let orbit = brown_batygin_orbit();
    let unit_mass = P9Params {
        mass_earth: 1.0,
        ..orbit
    };

    // The reference orbit laid on the sky, with the range signal a planet at
    // each point would leave in the Cassini residuals and the heaviest planet
    // that would stay under the ranging precision there.
    let track: Vec<Value> = (0..180)
        .map(|k| {
            let nu = (2.0 * k as f64).to_radians();
            let geom = position_at_true_anomaly(&orbit, nu);
            let (ra, dec) = ecliptic_vec_to_equatorial_deg(&geom.position);
            let residual_m = (range_residual(&orbit, nu) / meters(1.0)).value;
            let per_earth_mass = (range_residual(&unit_mass, nu) / meters(1.0)).value;
            json!({
                "nu_deg": nu.to_degrees(),
                "ra_deg": ra,
                "dec_deg": dec,
                "r_au": geom.distance,
                "residual_m": residual_m,
                "excluded": residual_m > CASSINI_RANGE_PRECISION_M,
                "max_mass_earth": CASSINI_RANGE_PRECISION_M / per_earth_mass,
            })
        })
        .collect();

    let sky = favored_sky_position(&orbit, N_SCAN);
    let favored_per_earth_mass = (range_residual(&unit_mass, sky.true_anomaly) / meters(1.0)).value;

    json!({
        "orbit": {"mass_earth": orbit.mass_earth, "a_au": orbit.a, "e": orbit.e},
        "precision_m": CASSINI_RANGE_PRECISION_M,
        "track": track,
        "favored": {
            "ra_deg": sky.ra_deg,
            "dec_deg": sky.dec_deg,
            "nu_deg": sky.true_anomaly.to_degrees(),
            "r_au": sky.distance_au,
            "max_mass_earth": CASSINI_RANGE_PRECISION_M / favored_per_earth_mass,
        },
        "favored_ra_deg": sky.ra_deg,
        "favored_dec_deg": sky.dec_deg,
        "offset_deg": sky.separation_deg(PREFERRED_RA_DEG, PREFERRED_DEC_DEG),
        "published": {
            "ra_deg": PREFERRED_RA_DEG,
            "dec_deg": PREFERRED_DEC_DEG,
            "half_extent_deg": PREFERRED_SKY_HALF_EXTENT_DEG,
            "tidal_range": PUBLISHED_TIDAL_RANGE,
            "tidal_unit_mass_earth": TIDAL_UNIT_MASS_EARTH,
            "tidal_unit_distance_au": TIDAL_UNIT_DISTANCE_AU,
        },
    })
}
