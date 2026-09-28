//! Film export for `p9-2023-mond-efe`: the numbers its scene and ledger entry draw.

use nalgebra::Vector3;
use p9_2023_mond_efe::efe::{
    MOND_A0_M_S2, efe_acceleration, galactic_center_ecliptic, galactic_center_ecliptic_lat,
    galactic_center_ecliptic_lon,
};
use p9_2023_mond_efe::secular::{
    PlaneGeometry, apsidal_line_separation, apsidal_torque, averaged_disturbing, forced_varpi,
    precession_rate,
};
use p9_core::analysis::circular::{circular_mean, mean_resultant_length};
use p9_core::constants::{AU_M, DAY_S, GM_SUN, TWO_PI, YEAR_DAYS};
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::units::degrees;
use serde_json::{Value, json};

/// Test orbit of the crate's headline tests.
const A_AU: f64 = 350.0;
const E0: f64 = 0.7;
const N_QUAD: usize = 128;
/// Half-width (AU) and points per side of the tidal-field grid.
const FIELD_HALF_AU: f64 = 600.0;
const FIELD_N: usize = 9;

/// The MOND transition radius r_M = sqrt(GM_sun / a0), in AU.
fn mond_radius_au(a0_au_day2: f64) -> f64 {
    (GM_SUN / a0_au_day2).sqrt()
}

/// Longitude offset from the axis (rad) at which the averaged disturbing
/// function of the e = `E0` orbit equals that of a circular orbit: the edge of
/// the lobe inside which an apsidal line cannot leave the axis.
fn lobe_half_width(a_efe: f64, n: &Vector3<f64>, geom: &PlaneGeometry, axis: f64) -> f64 {
    let floor = averaged_disturbing(A_AU, 0.0, axis, a_efe, n, geom, N_QUAD);
    let excess =
        |off: f64| averaged_disturbing(A_AU, E0, axis + off, a_efe, n, geom, N_QUAD) - floor;
    let (mut lo, mut hi) = (0.0, std::f64::consts::FRAC_PI_2);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if excess(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

pub fn export() -> Value {
    let n = galactic_center_ecliptic();
    let geom = PlaneGeometry::from_normal(Vector3::new(0.0, 0.0, 1.0));
    let gc_lon = geom.in_plane_longitude(&n);

    // The crate carries the tidal amplitude as a free parameter. The natural
    // MOND scale is a0 / r_M; the dimensionless factor multiplying it depends
    // on the interpolating function and is not computed by the crate.
    let a0 = MOND_A0_M_S2 * DAY_S * DAY_S / AU_M;
    let r_m = mond_radius_au(a0);
    let a_efe = a0 / r_m;

    let forced = forced_varpi(A_AU, E0, a_efe, &n, &geom, N_QUAD);
    let drift = precession_rate(A_AU, E0, forced, a_efe, &n, &geom, N_QUAD);
    let lobe = lobe_half_width(a_efe, &n, &geom, forced);

    // Averaged disturbing function and torque against apsidal longitude, in
    // units of A a^2.
    let scale = a_efe * A_AU * A_AU;
    let lon_deg: Vec<f64> = (0..=180).map(|k| 2.0 * k as f64).collect();
    let (rbar, torque): (Vec<f64>, Vec<f64>) = lon_deg
        .iter()
        .map(|&d| {
            let v = d.to_radians();
            (
                averaged_disturbing(A_AU, E0, v, a_efe, &n, &geom, N_QUAD) / scale,
                apsidal_torque(A_AU, E0, v, a_efe, &n, &geom, N_QUAD) / scale,
            )
        })
        .unzip();
    let rbar_circular = averaged_disturbing(A_AU, 0.0, forced, a_efe, &n, &geom, N_QUAD) / scale;

    // The tidal acceleration on a grid in the ecliptic plane, in units of A r.
    let mut field = Vec::new();
    for ix in 0..FIELD_N {
        for iy in 0..FIELD_N {
            let at = |k: usize| FIELD_HALF_AU * (2.0 * k as f64 / (FIELD_N - 1) as f64 - 1.0);
            let r = Vector3::new(at(ix), at(iy), 0.0);
            if r.norm() < 1.0 {
                continue;
            }
            let acc = efe_acceleration(a_efe, &n, &r) / (a_efe * FIELD_HALF_AU);
            field.push(json!({"x_au": r.x, "y_au": r.y, "ax": acc.x, "ay": acc.y}));
        }
    }

    // The observed sample against the predicted axis.
    let varpis: Vec<f64> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| o.longitude_of_perihelion())
        .collect();
    let seps: Vec<f64> = varpis
        .iter()
        .map(|&v| apsidal_line_separation(v, gc_lon).to_degrees())
        .collect();
    let objects: Vec<Value> = BROWN_2017_SAMPLE
        .iter()
        .zip(varpis.iter().zip(&seps))
        .map(|(o, (&v, &sep))| {
            json!({
                "name": o.name,
                "a_au": o.a,
                "e": o.e,
                "varpi_deg": v.to_degrees().rem_euclid(360.0),
                "line_sep_deg": sep,
                "in_lobe": sep < lobe.to_degrees(),
            })
        })
        .collect();
    let mean_varpi = circular_mean(&varpis).unwrap_or(0.0);
    let n_in_lobe = seps.iter().filter(|&&s| s < lobe.to_degrees()).count();

    json!({
        "a_au": A_AU,
        "e": E0,
        "a0_m_s2": MOND_A0_M_S2,
        "mond_radius_au": r_m,
        "a_efe_per_day2": a_efe,
        "gc_lon_deg": (galactic_center_ecliptic_lon() / degrees(1.0)).value,
        "gc_lat_deg": (galactic_center_ecliptic_lat() / degrees(1.0)).value,
        "forced_varpi_deg": forced.to_degrees(),
        "forced_line_sep_deg": apsidal_line_separation(forced, gc_lon).to_degrees(),
        "drift_period_myr": (TWO_PI / drift).abs() / YEAR_DAYS / 1e6,
        "lobe_half_width_deg": lobe.to_degrees(),
        "lobe_chance_fraction": lobe.to_degrees() / 90.0,
        "lon_deg": lon_deg,
        "rbar": rbar,
        "rbar_circular": rbar_circular,
        "torque": torque,
        "field": field,
        "field_half_au": FIELD_HALF_AU,
        "objects": objects,
        "n_objects": seps.len(),
        "n_in_lobe": n_in_lobe,
        "fraction_in_lobe": n_in_lobe as f64 / seps.len() as f64,
        "mean_varpi_deg": mean_varpi.to_degrees().rem_euclid(360.0),
        "r_bar_varpi": mean_resultant_length(&varpis),
        "mean_direction_sep_deg": apsidal_line_separation(mean_varpi, gc_lon).to_degrees(),
        "mean_line_sep_deg": seps.iter().sum::<f64>() / seps.len() as f64,
        "random_line_sep_deg": 45.0,
    })
}
