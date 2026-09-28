//! Film export for `p9-2025-russell-albedo`: the numbers its scene and ledger entry draw.

use p9_2025_russell_albedo::composition::{Composition, LARGE_COLD, SMALL_WARM};
use p9_2025_russell_albedo::photometry::{
    MOST_LIKELY_DISTANCE_AU, absolute_h, apparent_v, apparent_v_most_likely,
};
use p9_2025_russell_albedo::{
    PUBLISHED_ALBEDO_MAX, PUBLISHED_ALBEDO_MIN, PUBLISHED_APPARENT_BRIGHT,
    PUBLISHED_APPARENT_FAINT, PUBLISHED_ENVELOPE_FRAC_MAX, PUBLISHED_ENVELOPE_FRAC_MIN,
    PUBLISHED_H_BRIGHT, PUBLISHED_H_FAINT, PUBLISHED_MASS_EARTH, PUBLISHED_MASS_MINUS,
    PUBLISHED_MASS_PLUS, PUBLISHED_RADIUS_MAX_EARTH, PUBLISHED_RADIUS_MIN_EARTH,
};
use p9_core::analysis::photometry::{
    ALBEDO_NEPTUNE, absolute_magnitude, apparent_magnitude, mass_radius_neptunian, opposition_delta,
};
use p9_core::analysis::surveys::limiting_magnitude;
use p9_core::constants::{AU_KM, EARTH_RADIUS_KM, PC_AU};
use serde_json::{Value, json};

/// Angular diameter (arcsec) of a body of `radius_earth` seen from `distance_au`.
fn angular_diameter_arcsec(radius_earth: f64, distance_au: f64) -> f64 {
    2.0 * radius_earth * EARTH_RADIUS_KM / (distance_au * AU_KM) * PC_AU
}

fn endpoint(name: &str, comp: &Composition, distances: &[f64]) -> Value {
    json!({
        "name": name,
        "radius_earth": comp.radius_earth,
        "albedo": comp.albedo,
        "envelope_fraction": comp.envelope_fraction,
        "h": absolute_h(comp),
        "v": apparent_v_most_likely(comp),
        "v_curve": distances.iter().map(|&d| apparent_v(comp, d)).collect::<Vec<_>>(),
        "disk_arcsec": angular_diameter_arcsec(comp.radius_earth, MOST_LIKELY_DISTANCE_AU),
    })
}

pub fn export() -> Value {
    let distances: Vec<f64> = (0..=70).map(|k| 300.0 + 10.0 * k as f64).collect();

    // The assumption the survey papers used before: a scaled-down Neptune.
    let neptune_radius = mass_radius_neptunian(PUBLISHED_MASS_EARTH);
    let neptune_h = absolute_magnitude(neptune_radius * EARTH_RADIUS_KM, ALBEDO_NEPTUNE);
    let neptune_v = |d: f64| apparent_magnitude(neptune_h, d, opposition_delta(d));

    let bright = endpoint("larger, colder, more reflective", &LARGE_COLD, &distances);
    let faint = endpoint("smaller, warmer, darker", &SMALL_WARM, &distances);

    json!({
        "distance_au": distances,
        "most_likely_distance_au": MOST_LIKELY_DISTANCE_AU,
        "bright": bright,
        "faint": faint,
        "h_bright": absolute_h(&LARGE_COLD),
        "h_faint": absolute_h(&SMALL_WARM),
        "v_bright": apparent_v_most_likely(&LARGE_COLD),
        "v_faint": apparent_v_most_likely(&SMALL_WARM),
        "scaled_neptune": {
            "mass_earth": PUBLISHED_MASS_EARTH,
            "radius_earth": neptune_radius,
            "albedo": ALBEDO_NEPTUNE,
            "h": neptune_h,
            "v": neptune_v(MOST_LIKELY_DISTANCE_AU),
            "v_curve": distances.iter().map(|&d| neptune_v(d)).collect::<Vec<_>>(),
        },
        "depths": {
            "ps1": limiting_magnitude("PS1 3pi"),
            "des": limiting_magnitude("DES"),
        },
        "published": {
            "mass_earth": PUBLISHED_MASS_EARTH,
            "mass_plus": PUBLISHED_MASS_PLUS,
            "mass_minus": PUBLISHED_MASS_MINUS,
            "radius_earth": [PUBLISHED_RADIUS_MIN_EARTH, PUBLISHED_RADIUS_MAX_EARTH],
            "envelope_fraction": [PUBLISHED_ENVELOPE_FRAC_MIN, PUBLISHED_ENVELOPE_FRAC_MAX],
            "albedo": [PUBLISHED_ALBEDO_MIN, PUBLISHED_ALBEDO_MAX],
            "h": [PUBLISHED_H_BRIGHT, PUBLISHED_H_FAINT],
            "v": [PUBLISHED_APPARENT_BRIGHT, PUBLISHED_APPARENT_FAINT],
        },
    })
}
