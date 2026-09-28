//! Film export for `p9-2016-obliquity-gomes`: the numbers its scene and ledger entry draw.

use p9_2016_obliquity_gomes::precession_model::{
    GomesParams, forced_inclination_typed, l_p9, l_planets, max_obliquity_deg,
    obliquity_over_age_deg, obliquity_series, precession_period_typed,
};
use p9_2016_obliquity_gomes::reference::SOLAR_SYSTEM_AGE_GYR;
use p9_core::constants::GYR_DAYS;
use p9_core::units::{days, radians};
use serde_json::{Value, json};
use std::f64::consts::PI;

/// Tilt of the giant planets' invariant plane to the solar equator that Gomes
/// et al. (2016) set out to match (degrees).
const TARGET_TILT_DEG: f64 = 5.9;

/// Semi-major axes at which the paper quotes the eccentricity needed by a
/// 10 Earth-mass, 30 degree perturber, with those eccentricities.
const PUBLISHED_REQUIRED_E: [(f64, f64); 2] = [(600.0, 0.71), (700.0, 0.80)];

/// Eccentricity of the Batygin & Brown (2016) perturber at 700 AU.
const NOMINAL_E: f64 = 0.6;

/// Planet Nine masses whose tilt curves are compared (Earth masses).
const MASSES_EARTH: [f64; 3] = [10.0, 15.0, 20.0];

fn case(p: &GomesParams) -> Value {
    let i_p = (forced_inclination_typed(p) / radians(1.0)).value;
    let period_days = (precession_period_typed(p) / days(1.0)).value;
    let age_days = SOLAR_SYSTEM_AGE_GYR * GYR_DAYS;
    let series = obliquity_series(p, age_days, 90);

    // Track of the planetary-plane pole about the total angular momentum: a
    // circle of radius i_p through the primordial pole, where the solar spin
    // axis stays. Coordinates are degrees, centred on the total angular
    // momentum, with the primordial pole on the +x axis.
    let r = i_p.to_degrees();
    let pole: Vec<(f64, f64)> = series
        .iter()
        .map(|s| {
            let phase = 2.0 * PI * s.t / period_days;
            (r * phase.cos(), -r * phase.sin())
        })
        .collect();

    json!({
        "mass_earth": p.m9_earth,
        "a_au": p.a9,
        "e": p.e9,
        "i9_deg": p.i9.to_degrees(),
        "forced_inclination_deg": r,
        "precession_period_gyr": period_days / GYR_DAYS,
        "phase_swept_deg": 360.0 * age_days / period_days,
        "max_tilt_deg": max_obliquity_deg(p),
        "tilt_today_deg": obliquity_over_age_deg(p),
        "angular_momentum_ratio": l_p9(p) / l_planets(),
        "t_gyr": series.iter().map(|s| s.t / GYR_DAYS).collect::<Vec<_>>(),
        "tilt_deg": series.iter().map(|s| s.obliquity.to_degrees()).collect::<Vec<_>>(),
        "pole_track_deg": pole,
    })
}

/// Tilt after 4.5 Gyr against eccentricity, and the eccentricity at which it
/// first reaches the target.
fn eccentricity_scan(base: &GomesParams) -> Value {
    let e: Vec<f64> = (0..=180).map(|k| 0.005 * k as f64).collect();
    let tilt: Vec<f64> = e
        .iter()
        .map(|&e9| obliquity_over_age_deg(&GomesParams { e9, ..*base }))
        .collect();
    let required = e
        .windows(2)
        .zip(tilt.windows(2))
        .find(|(_, t)| t[0] < TARGET_TILT_DEG && t[1] >= TARGET_TILT_DEG)
        .map(|(e, t)| e[0] + (e[1] - e[0]) * (TARGET_TILT_DEG - t[0]) / (t[1] - t[0]));
    json!({
        "mass_earth": base.m9_earth,
        "a_au": base.a9,
        "i9_deg": base.i9.to_degrees(),
        "e": e,
        "tilt_deg": tilt,
        "required_e": required,
    })
}

pub fn export() -> Value {
    let nominal = GomesParams {
        e9: NOMINAL_E,
        ..GomesParams::nominal()
    };

    let scans: Vec<Value> = PUBLISHED_REQUIRED_E
        .iter()
        .map(|&(a9, published_e)| {
            let mut scan = eccentricity_scan(&GomesParams { a9, ..nominal });
            scan["published_required_e"] = json!(published_e);
            scan
        })
        .collect();
    let mass_scans: Vec<Value> = MASSES_EARTH
        .iter()
        .map(|&m9_earth| {
            eccentricity_scan(&GomesParams {
                m9_earth,
                ..nominal
            })
        })
        .collect();

    let required_e_700 = scans
        .iter()
        .find(|s| s["a_au"].as_f64() == Some(nominal.a9))
        .and_then(|s| s["required_e"].as_f64());
    let solved = required_e_700.map(|e9| case(&GomesParams { e9, ..nominal }));

    json!({
        "target_tilt_deg": TARGET_TILT_DEG,
        "age_gyr": SOLAR_SYSTEM_AGE_GYR,
        "nominal": case(&nominal),
        "nominal_tilt_deg": obliquity_over_age_deg(&nominal),
        "solved": solved,
        "required_e_700": required_e_700,
        "eccentricity_scans": scans,
        "mass_scans": mass_scans,
    })
}
