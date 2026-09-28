//! Film export for `scenes/preface/preface_e_detection.py`: the numbers its scenes draw.
//!
//! - `reflected`: the V magnitude of a Planet Nine vs distance, the Neptune and
//!   Pluto anchors, and where each survey's depth line cuts the curve (P09).
//! - `thermal`: sunlight-set vs internal-heat temperature, a family of Planck
//!   spectra from the Sun's 5778 K down to 40 K, and reflected vs thermal flux
//!   fall-off with distance (P10).
//! - `movers`: apparent drift rate vs distance near opposition (P11Movers).
//! - `bias`: a synthetic isotropic population pushed through a survey wedge and
//!   a depth cut, showing that detections cluster where the survey looked (P11Bias).

use rand::{Rng, SeedableRng};
use serde_json::{Value, json};

use p9_2023_lsst_strategy::strategy::published::SINGLE_VISIT_DEPTH_R;
use p9_2025_iras_akari::survey_model::AkariFisSurvey;
use p9_core::analysis::circular::mean_resultant_length;
use p9_core::analysis::hansen::mean_to_true_anomaly;
use p9_core::analysis::photometry::{
    ALBEDO_NEPTUNE, NEPTUNE_MASS_EARTH, SOLAR_V_MINUS_R, absolute_magnitude, apparent_magnitude,
    mass_radius_neptunian, opposition_delta, planet_apparent_magnitude,
};
use p9_core::analysis::stacking::orbit_metric::apparent_sky_rate_at_opposition;
use p9_core::analysis::surveys::limiting_magnitude;
use p9_core::analysis::thermal::{
    C_LIGHT, H_PLANCK, K_BOLTZ, T_SUN, WIEN_X_NU, effective_temp, max_detectable_distance,
    planck_bnu, reflected_flux_jy, solar_equilibrium_temp, thermal_flux_jy,
};
use p9_core::constants::EARTH_RADIUS_KM;
use p9_core::types::{P9Params, solve_kepler};

/// Internal-heat temperature floor of a cold ~5-10 Earth-mass giant (K), the
/// value the workspace's thermal reproductions and the viability map adopt.
const INTERNAL_TEMP_K: f64 = 40.0;
/// Pluto's mean radius (km, New Horizons) and V geometric albedo.
const PLUTO_RADIUS_KM: f64 = 1188.3;
const PLUTO_ALBEDO: f64 = 0.575;
/// Pluto's heliocentric distance in the mid-2020s (AU).
const PLUTO_NOW_AU: f64 = 35.0;
/// Neptune's semimajor axis (AU) and observed effective temperature (K).
const NEPTUNE_AU: f64 = 30.07;
const NEPTUNE_T_OBS_K: f64 = 59.0;
/// The Moon's mean apparent diameter (arcsec), the everyday angular anchor.
const MOON_DIAMETER_ARCSEC: f64 = 1865.0;

fn log_grid(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    let (a, b) = (lo.log10(), hi.log10());
    (0..n)
        .map(|k| 10f64.powf(a + (b - a) * k as f64 / (n - 1) as f64))
        .collect()
}

fn radius_m(mass_earth: f64) -> f64 {
    mass_radius_neptunian(mass_earth) * EARTH_RADIUS_KM * 1e3
}

/// Reflected-light magnitudes, anchors, and survey reaches for P09.
fn reflected() -> Value {
    let mass = P9Params::mcmc_2021().mass_earth;
    let v_at = |d: f64| planet_apparent_magnitude(mass, ALBEDO_NEPTUNE, d);
    let dist = log_grid(20.0, 1500.0, 160);
    let v_mag: Vec<f64> = dist.iter().map(|&d| v_at(d)).collect();

    // Survey depths as V-equivalent limits (r-band depths shifted by the
    // solar V - r colour, so they sit on the same axis as the V curve).
    let ztf_v = limiting_magnitude("ZTF").expect("ZTF depth");
    let ps1_r = limiting_magnitude("PS1 3pi").expect("PS1 depth");
    let des_r = limiting_magnitude("DES").expect("DES depth");
    let depths = [
        ("ZTF", "V", ztf_v, ztf_v),
        ("Pan-STARRS", "r", ps1_r, ps1_r + SOLAR_V_MINUS_R),
        ("DES", "r", des_r, des_r + SOLAR_V_MINUS_R),
        (
            "Rubin (1 visit)",
            "r",
            SINGLE_VISIT_DEPTH_R,
            SINGLE_VISIT_DEPTH_R + SOLAR_V_MINUS_R,
        ),
    ];
    let surveys: Vec<Value> = depths
        .iter()
        .map(|&(name, band, depth, v_lim)| {
            json!({
                "name": name,
                "band": band,
                "depth": depth,
                "v_limit": v_lim,
                "reach_au": max_detectable_distance(20.0, 5000.0, v_lim, v_at),
            })
        })
        .collect();

    let pluto_h = absolute_magnitude(PLUTO_RADIUS_KM, PLUTO_ALBEDO);
    json!({
        "mass_earth": mass,
        "albedo": ALBEDO_NEPTUNE,
        "radius_earth": mass_radius_neptunian(mass),
        "distance_au": dist,
        "v_mag": v_mag,
        "sunlight_vs_earth": dist.iter().map(|d| 1.0 / (d * d)).collect::<Vec<_>>(),
        "neptune": {
            "r_au": NEPTUNE_AU,
            "v": planet_apparent_magnitude(NEPTUNE_MASS_EARTH, ALBEDO_NEPTUNE, NEPTUNE_AU),
        },
        "pluto": {
            "r_au": PLUTO_NOW_AU,
            "v": apparent_magnitude(pluto_h, PLUTO_NOW_AU, opposition_delta(PLUTO_NOW_AU)),
        },
        "v_300": v_at(300.0),
        "v_600": v_at(600.0),
        "surveys": surveys,
    })
}

/// Temperatures, Planck spectra and flux fall-off for P10.
fn thermal() -> Value {
    // Temperature set by sunlight alone vs the internal-heat floor.
    let dist = log_grid(10.0, 1500.0, 120);
    let t_eq: Vec<f64> = dist
        .iter()
        .map(|&d| solar_equilibrium_temp(d, ALBEDO_NEPTUNE))
        .collect();
    // Distance where sunlight alone can no longer hold 40 K.
    let crossover = max_detectable_distance(1.0, 5000.0, -INTERNAL_TEMP_K, |d| {
        -solar_equilibrium_temp(d, ALBEDO_NEPTUNE)
    });

    // Planck family: B_nu normalised to its own peak, from the Sun to 40 K.
    let wl_um = log_grid(0.1, 3000.0, 180);
    let temps = log_grid(T_SUN, INTERNAL_TEMP_K, 60);
    let spectra: Vec<Vec<f64>> = temps
        .iter()
        .map(|&t| {
            let raw: Vec<f64> = wl_um
                .iter()
                .map(|um| planck_bnu(t, C_LIGHT / (um * 1e-6)))
                .collect();
            let peak = planck_bnu(t, WIEN_X_NU * K_BOLTZ * t / H_PLANCK);
            raw.iter().map(|b| b / peak).collect()
        })
        .collect();
    let peak_um: Vec<f64> = temps
        .iter()
        .map(|&t| 1e6 * C_LIGHT * H_PLANCK / (WIEN_X_NU * K_BOLTZ * t))
        .collect();

    // Reflected (V band) vs thermal (AKARI 90 um) flux of the same planet,
    // each normalised to its value at 100 AU.
    let mass = P9Params::mcmc_2021().mass_earth;
    let r_m = radius_m(mass);
    let nu_v = C_LIGHT / 0.55e-6;
    let nu_fir = C_LIGHT / 90e-6;
    let fdist = log_grid(100.0, 1500.0, 80);
    let refl = |d: f64| reflected_flux_jy(ALBEDO_NEPTUNE, r_m, d, nu_v);
    let therm = |d: f64| {
        thermal_flux_jy(
            effective_temp(d, ALBEDO_NEPTUNE, INTERNAL_TEMP_K),
            r_m,
            d,
            nu_fir,
        )
    };
    let (r0, t0) = (refl(100.0), therm(100.0));

    json!({
        "distance_au": dist,
        "t_eq": t_eq,
        "t_internal": INTERNAL_TEMP_K,
        "crossover_au": crossover,
        "t_eq_600": solar_equilibrium_temp(600.0, ALBEDO_NEPTUNE),
        "neptune": {
            "r_au": NEPTUNE_AU,
            "t_obs": NEPTUNE_T_OBS_K,
            "t_eq": solar_equilibrium_temp(NEPTUNE_AU, ALBEDO_NEPTUNE),
        },
        "wavelength_um": wl_um,
        "temps": temps,
        "spectra": spectra,
        "peak_um": peak_um,
        "bands": [
            {"name": "eye", "um": 0.55},
            {"name": "WISE", "um": 3.4},
            {"name": "IRAS", "um": 60.0},
            {"name": "AKARI", "um": 90.0},
            {"name": "ACT (mm)", "um": 1e6 * C_LIGHT / 229e9},
        ],
        "akari_limit_jy": AkariFisSurvey::default().sensitivity_jy,
        "flux_distance_au": fdist,
        "reflected_rel": fdist.iter().map(|&d| refl(d) / r0).collect::<Vec<_>>(),
        "thermal_rel": fdist.iter().map(|&d| therm(d) / t0).collect::<Vec<_>>(),
        "reflected_jy_600": refl(600.0),
        "thermal_jy_600": therm(600.0),
    })
}

/// Apparent drift rate near opposition vs heliocentric distance (P11Movers).
fn movers() -> Value {
    let dist = log_grid(20.0, 1500.0, 120);
    let rate: Vec<f64> = dist
        .iter()
        .map(|&d| apparent_sky_rate_at_opposition(d))
        .collect();
    let anchor = |name: &str, d: f64| {
        let rate = apparent_sky_rate_at_opposition(d);
        json!({
            "name": name,
            "r_au": d,
            "arcsec_per_day": rate,
            "days_per_moon_width": MOON_DIAMETER_ARCSEC / rate,
        })
    };
    json!({
        "distance_au": dist,
        "arcsec_per_day": rate,
        "anchors": [
            anchor("Pluto", PLUTO_NOW_AU),
            anchor("Sedna", 83.0),
            anchor("Planet Nine", 600.0),
        ],
        "moon_arcsec": MOON_DIAMETER_ARCSEC,
    })
}

/// Synthetic selection-bias experiment (P11Bias).
///
/// An isotropic (uniform longitude of perihelion) planar population of
/// distant objects is observed by a survey that only looks within a wedge of
/// ecliptic longitude and only down to a limiting magnitude. Magnitudes use
/// the core H + 5 log(r Delta) law at opposition; positions come from the
/// core Kepler solver at a uniform-random mean anomaly (the object's phase
/// today). Detected objects are overwhelmingly caught near perihelion, so
/// their perihelion longitudes pile up inside the wedge.
fn bias() -> Value {
    const N: usize = 30_000;
    const H: f64 = 6.0;
    const DEPTH: f64 = 24.0;
    let wedges = [(40.0_f64, 25.0_f64), (220.0, 25.0)]; // (centre, half-width) deg
    let mut rng = rand::rngs::StdRng::seed_from_u64(9);

    let in_wedge = |lon_deg: f64, (c, hw): (f64, f64)| -> bool {
        let d = (lon_deg - c + 540.0).rem_euclid(360.0) - 180.0;
        d.abs() <= hw
    };

    let mut objects = Vec::with_capacity(N);
    let mut detected: Vec<Vec<usize>> = vec![Vec::new(); wedges.len()];
    for k in 0..N {
        let a = 10f64.powf(rng.gen_range(250f64.log10()..700f64.log10()));
        let q = rng.gen_range(33.0..50.0);
        let e = 1.0 - q / a;
        let varpi = rng.gen_range(0.0..360.0_f64);
        let m = rng.gen_range(0.0..std::f64::consts::TAU);
        let r = a * (1.0 - e * solve_kepler(e, m).cos());
        let nu = (mean_to_true_anomaly(m, e) + std::f64::consts::PI)
            .rem_euclid(std::f64::consts::TAU)
            - std::f64::consts::PI;
        let lon = (varpi + nu.to_degrees()).rem_euclid(360.0);
        let v = apparent_magnitude(H, r, opposition_delta(r));
        for (w, &wedge) in wedges.iter().enumerate() {
            if v <= DEPTH && in_wedge(lon, wedge) {
                detected[w].push(k);
            }
        }
        objects.push(json!({
            "a": a, "e": e, "varpi_deg": varpi, "nu_deg": nu.to_degrees(),
            "r_au": r, "lon_deg": lon, "v": v,
        }));
    }

    // Histogram of perihelion longitude: everything vs what each wedge caught.
    let nbins = 18;
    let hist = |idx: &mut dyn Iterator<Item = f64>| -> Vec<usize> {
        let mut h = vec![0usize; nbins];
        for v in idx {
            h[((v / 360.0 * nbins as f64) as usize).min(nbins - 1)] += 1;
        }
        h
    };
    let all_hist = hist(&mut objects.iter().map(|o| o["varpi_deg"].as_f64().unwrap()));
    let wedge_json: Vec<Value> = wedges
        .iter()
        .zip(&detected)
        .map(|(&(c, hw), idx)| {
            let h = hist(
                &mut idx
                    .iter()
                    .map(|&k| objects[k]["varpi_deg"].as_f64().unwrap()),
            );
            let near_peri = idx
                .iter()
                .filter(|&&k| objects[k]["nu_deg"].as_f64().unwrap().abs() < 45.0)
                .count();
            let varpis: Vec<f64> = idx
                .iter()
                .map(|&k| objects[k]["varpi_deg"].as_f64().unwrap().to_radians())
                .collect();
            json!({
                "r_bar": mean_resultant_length(&varpis),
                "centre_deg": c,
                "half_width_deg": hw,
                "detected": idx.iter().map(|&k| objects[k].clone()).collect::<Vec<_>>(),
                "hist": h,
                "near_perihelion_frac": near_peri as f64 / idx.len().max(1) as f64,
            })
        })
        .collect();

    json!({
        "n": N,
        "h": H,
        "depth": DEPTH,
        // heliocentric distance inside which an H = 6 body beats the depth
        "r_limit_au": max_detectable_distance(1.5, 5000.0, DEPTH, |r| {
            apparent_magnitude(H, r, opposition_delta(r))
        }),
        "background": objects.iter().take(160).cloned().collect::<Vec<_>>(),
        "bins": nbins,
        "all_hist": all_hist,
        "all_r_bar": mean_resultant_length(
            &objects
                .iter()
                .map(|o| o["varpi_deg"].as_f64().unwrap().to_radians())
                .collect::<Vec<_>>(),
        ),
        "wedges": wedge_json,
    })
}

pub fn export() -> Value {
    json!({
        "reflected": reflected(),
        "thermal": thermal(),
        "movers": movers(),
        "bias": bias(),
    })
}
