//! Film export for `p9-2025-simons-forecast`: the numbers its scene and ledger entry draw.

use p9_2025_simons_forecast::bands::{
    ACT_SENSITIVITY, MmSensitivity, SO_SENSITIVITY, SO_SENSITIVITY_SHALLOW,
};
use p9_2025_simons_forecast::forecast::{
    ReferenceBox, detection_significance, max_detectable_distance,
};
use p9_2025_simons_forecast::published;
use p9_2025_simons_forecast::thermal::P9Thermal;
use p9_core::data::reference_population::{generate_reference_population, heliocentric_distance};
use p9_core::types::P9Params;
use p9_core::units::au;
use rand::SeedableRng;
use serde_json::{Value, json};

/// Reference-population members drawn on the mass-distance plane.
const N_POP: usize = 500;

fn reach_au(mass_earth: f64, sens: &MmSensitivity) -> f64 {
    (max_detectable_distance(mass_earth, sens) / au(1.0)).value
}

pub fn export() -> Value {
    let surveys = [
        ("act", &ACT_SENSITIVITY),
        ("so_shallow", &SO_SENSITIVITY_SHALLOW),
        ("so_deep", &SO_SENSITIVITY),
    ];
    let bx = ReferenceBox::nominal();

    // The Brown & Batygin (2021) reference population on the mass-distance
    // plane, for context.
    let mut rng = rand::rngs::StdRng::seed_from_u64(2025);
    let population: Vec<(f64, f64)> = generate_reference_population(N_POP, &mut rng)
        .iter()
        .map(|obj| {
            let r = heliocentric_distance(&P9Params {
                mass_earth: obj.mass,
                a: obj.a,
                e: obj.e,
                i: obj.i,
                omega: obj.omega,
                omega_big: obj.omega_big,
                mean_anomaly: obj.mean_anomaly,
            });
            (obj.mass, r)
        })
        .collect();

    let masses: Vec<f64> = (0..=52).map(|k| 2.0 + 0.25 * k as f64).collect();
    let distances: Vec<f64> = (0..=110).map(|k| 250.0 + 10.0 * k as f64).collect();

    let mut out = serde_json::Map::new();
    for (key, sens) in surveys {
        out.insert(
            key.to_string(),
            json!({
                "label": sens.survey,
                "nu_ghz": sens.nu_hz / 1e9,
                "flux_limit_mjy": sens.flux_limit_mjy,
                "reach_5": reach_au(5.0, sens),
                "reach_10": reach_au(10.0, sens),
                "reach_au": masses.iter().map(|&m| reach_au(m, sens)).collect::<Vec<_>>(),
                "snr_5": distances
                    .iter()
                    .map(|&d| detection_significance(5.0, d, sens))
                    .collect::<Vec<_>>(),
                "box_fraction": bx.detectable_fraction(sens),
            }),
        );
    }

    let body = P9Thermal::new(5.0, 500.0);
    out.insert(
        "so_deep_reach_5".into(),
        json!(reach_au(5.0, &SO_SENSITIVITY)),
    );
    out.insert("act_reach_5".into(), json!(reach_au(5.0, &ACT_SENSITIVITY)));
    out.insert("mass_earth".into(), json!(masses));
    out.insert("distance_au".into(), json!(distances));
    out.insert(
        "body".into(),
        json!({
            "temp_k": body.effective_temp(),
            "radius_earth": body.radius_m() / p9_core::analysis::thermal::R_EARTH_M,
            "flux_500_mjy": body.flux_mjy(SO_SENSITIVITY.nu_hz),
            "flux_900_mjy": P9Thermal::new(5.0, 900.0).flux_mjy(SO_SENSITIVITY.nu_hz),
        }),
    );
    out.insert(
        "reference_box".into(),
        json!({
            "mass_min": bx.mass_min_earth,
            "mass_max": bx.mass_max_earth,
            "dist_min": bx.dist_min_au,
            "dist_max": bx.dist_max_au,
        }),
    );
    out.insert(
        "population".into(),
        json!(
            population
                .iter()
                .map(|(m, r)| json!({"mass_earth": m, "dist_au": r}))
                .collect::<Vec<_>>()
        ),
    );
    out.insert(
        "published".into(),
        json!({
            "so_deep_5": published::SO_REACH_5ME_DEEP_AU,
            "so_shallow_5": published::SO_REACH_5ME_SHALLOW_AU,
            "act_5_min": published::ACT_REACH_5ME_MIN_AU,
            "act_5_max": published::ACT_REACH_5ME_MAX_AU,
            "act_10_min": published::ACT_REACH_10ME_MIN_AU,
            "act_10_max": published::ACT_REACH_10ME_MAX_AU,
        }),
    );
    Value::Object(out)
}
