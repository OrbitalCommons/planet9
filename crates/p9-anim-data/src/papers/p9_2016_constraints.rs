//! Film export for `p9-2016-constraints`: the numbers its scene and ledger entry draw.

use p9_2016_constraints::clustering_metric::{observed_etno_r_bar, observed_six_kbo_r_bar};
use p9_2016_constraints::detection_limits::{
    brightness_curve, detectable_fraction, estimate_radius, predicted_v_magnitude, sky_position,
    survey_depths,
};
use p9_2016_constraints::parameter_grid::generate_grid;
use p9_core::analysis::hansen::mean_to_true_anomaly;
use p9_core::coords::candidate_pair::orbital_period_years;
use p9_core::coords::sky::ecliptic_to_equatorial;
use p9_core::types::{P9Params, true_to_mean_anomaly};
use p9_core::units::km;
use serde_json::{Value, json};
use std::f64::consts::PI;

use super::p9_2016_evidence::p9_orbit_json;

/// Geometric albedo assumed for the brightness curve.
const ALBEDO: f64 = 0.5;
/// Samples around the orbit.
const N_ORBIT: usize = 180;
/// Equal-time snapshots around one period.
const N_SNAPSHOT: usize = 60;
/// Mass slice of the paper's grid that the scene draws.
const GRID_MASS: f64 = 10.0;
/// The wide survey whose depth splits the orbit into searched and unsearched.
const WIDE_SURVEY: &str = "Pan-STARRS 3π";

pub fn export() -> Value {
    let p9 = P9Params::nominal_2016();
    let wide_depth = survey_depths()
        .iter()
        .find(|s| s.0 == WIDE_SURVEY)
        .map(|s| s.1)
        .expect("wide survey is in the crate's depth table");

    // Brightness and sky position around the nominal orbit, with the share of
    // the orbital period spent in each step (the planet lingers at aphelion).
    let curve = brightness_curve(&p9, ALBEDO, N_ORBIT);
    let step = 2.0 * PI / N_ORBIT as f64;
    let mut time_total = 0.0;
    let mut time_faint = 0.0;
    let orbit: Vec<Value> = curve
        .iter()
        .map(|&(nu, r, v_mag)| {
            let (lon, lat) = sky_position(&p9, nu);
            let (ra, dec) = ecliptic_to_equatorial(lon, lat);
            let dm = (true_to_mean_anomaly(p9.e, nu + 0.5 * step)
                - true_to_mean_anomaly(p9.e, nu - 0.5 * step))
            .rem_euclid(2.0 * PI);
            time_total += dm;
            if v_mag > wide_depth {
                time_faint += dm;
            }
            json!({
                "nu_deg": nu.to_degrees(),
                "r_au": r,
                "v_mag": v_mag,
                "ra_deg": ra.to_degrees().rem_euclid(360.0),
                "dec_deg": dec.to_degrees(),
                "time_share": dm / (2.0 * PI),
            })
        })
        .collect();

    let surveys: Vec<Value> = survey_depths()
        .iter()
        .map(|&(name, depth)| {
            json!({
                "name": name,
                "depth": depth,
                "orbit_fraction": detectable_fraction(&p9, depth, ALBEDO),
            })
        })
        .collect();

    // Snapshots equally spaced in time (mean anomaly): where the planet is
    // after each 1/N_SNAPSHOT of its period. They bunch up at aphelion.
    let radius_km = (estimate_radius(p9.mass_earth) / km(1.0)).value;
    let snapshots: Vec<Value> = (0..N_SNAPSHOT)
        .map(|k| {
            let mean_anomaly = 2.0 * PI * (k as f64 + 0.5) / N_SNAPSHOT as f64;
            let nu = mean_to_true_anomaly(mean_anomaly, p9.e);
            let (lon, lat) = sky_position(&p9, nu);
            let (ra, dec) = ecliptic_to_equatorial(lon, lat);
            json!({
                "ra_deg": ra.to_degrees().rem_euclid(360.0),
                "dec_deg": dec.to_degrees(),
                "v_mag": predicted_v_magnitude(&p9, nu, ALBEDO, radius_km),
            })
        })
        .collect();

    let v_peri = curve.iter().map(|c| c.2).fold(f64::INFINITY, f64::min);
    let v_apo = curve.iter().map(|c| c.2).fold(f64::NEG_INFINITY, f64::max);

    // One mass slice of the paper's planar survey grid.
    let grid = generate_grid();
    let slice: Vec<Value> = grid
        .iter()
        .filter(|g| (g.mass_earth - GRID_MASS).abs() < 1e-9)
        .map(|g| json!({"a": g.a, "e": g.e, "q": g.perihelion}))
        .collect();

    json!({
        "mass": p9.mass_earth,
        "a": p9.a,
        "e": p9.e,
        "i": p9.i.to_degrees(),
        "p9": p9_orbit_json("Planet Nine", &p9),
        "albedo": ALBEDO,
        "orbit": orbit,
        "v_perihelion": v_peri,
        "v_aphelion": v_apo,
        "wide_survey": WIDE_SURVEY,
        "wide_depth": wide_depth,
        "orbit_fraction_wide": detectable_fraction(&p9, wide_depth, ALBEDO),
        "time_share_too_faint": time_faint / time_total,
        "snapshots": snapshots,
        "period_yr": orbital_period_years(p9.a),
        "surveys": surveys,
        "grid": slice,
        "grid_total": grid.len(),
        "grid_mass": GRID_MASS,
        "r_bar_six": observed_six_kbo_r_bar(),
        "r_bar_ten": observed_etno_r_bar(),
    })
}
