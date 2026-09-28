//! Film export for `p9-2016-linder-evolution`: the numbers its scene and ledger entry draw.

use p9_2016_linder_evolution::reference::{
    Q_10ME_700AU, Q_VS_MASS_APHELION, T_EFF_10ME_700AU_K, V_10ME_700AU, V_VS_MASS_APHELION,
    W4_10ME_700AU,
};
use p9_2016_linder_evolution::reflected_v_magnitude;
use p9_2016_linder_evolution::thermal::{Band, P9};
use p9_core::analysis::thermal::R_EARTH_M;
use serde_json::{Value, json};

/// The paper's orbit: perihelion, semi-major axis and aphelion distances (AU).
const PERIHELION_AU: f64 = 280.0;
const NOMINAL_AU: f64 = 700.0;
const APHELION_AU: f64 = 1120.0;

/// Published radius of the nominal 10 Earth-mass planet (Earth radii).
const PUBLISHED_RADIUS_EARTH: f64 = 3.7;

pub fn export() -> Value {
    let masses: Vec<f64> = V_VS_MASS_APHELION.iter().map(|&(m, _)| m).collect();
    let distances: Vec<f64> = (0..=50).map(|k| 200.0 + 20.0 * k as f64).collect();

    // Reflected-light V and intrinsic Q against distance, one curve per mass.
    let curves: Vec<Value> = masses
        .iter()
        .map(|&m| {
            let v: Vec<f64> = distances
                .iter()
                .map(|&d| reflected_v_magnitude(m, d))
                .collect();
            let q: Vec<f64> = distances
                .iter()
                .map(|&d| {
                    P9 {
                        mass_earth: m,
                        distance_au: d,
                    }
                    .thermal_magnitude(Band::Q)
                })
                .collect();
            json!({"mass_earth": m, "v_mag": v, "q_mag": q})
        })
        .collect();

    // The aphelion table: reproduced against published, per mass.
    let aphelion: Vec<Value> = V_VS_MASS_APHELION
        .iter()
        .zip(Q_VS_MASS_APHELION.iter())
        .map(|(&(m, v_pub), &(_, q_pub))| {
            let body = P9 {
                mass_earth: m,
                distance_au: APHELION_AU,
            };
            json!({
                "mass_earth": m,
                "t_eff_k": body.effective_temperature(),
                "radius_earth": body.radius_m() / R_EARTH_M,
                "v_mag": reflected_v_magnitude(m, APHELION_AU),
                "v_published": v_pub,
                "q_mag": body.thermal_magnitude(Band::Q),
                "q_published": q_pub,
            })
        })
        .collect();

    let nominal = P9 {
        mass_earth: 10.0,
        distance_au: NOMINAL_AU,
    };

    json!({
        "perihelion_au": PERIHELION_AU,
        "nominal_au": NOMINAL_AU,
        "aphelion_au": APHELION_AU,
        "distance_au": distances,
        "curves": curves,
        "aphelion": aphelion,
        "t_eff_10me_k": nominal.effective_temperature(),
        "t_eff_published_k": T_EFF_10ME_700AU_K,
        "radius_10me_earth": nominal.radius_m() / R_EARTH_M,
        "radius_published_earth": PUBLISHED_RADIUS_EARTH,
        "v_10me_700au": reflected_v_magnitude(10.0, NOMINAL_AU),
        "v_10me_700au_published": V_10ME_700AU,
        "q_10me_700au": nominal.thermal_magnitude(Band::Q),
        "q_10me_700au_published": Q_10ME_700AU,
        "w4_10me_700au": nominal.thermal_magnitude(Band::W4),
        "w4_10me_700au_published": W4_10ME_700AU,
    })
}
