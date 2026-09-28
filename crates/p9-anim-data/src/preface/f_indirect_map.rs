//! Film export for `scenes/preface/preface_f_indirect_map.py`: the numbers its scenes draw.
//!
//! - `ranging`: Planet Nine's tidal tug on Saturn and the Cassini range residual
//!   it would leave, around the Batygin & Brown (2016) orbit (P12Ranging).
//! - `obliquity`: the Bailey, Batygin & Brown (2016) secular spin-orbit
//!   integration over 4.5 Gyr, as the tips of the three axes seen down the
//!   total angular momentum, plus the solar obliquity vs time (P12Obliquity).
//! - `map`: Brown & Batygin (2021) posterior draws in mass x distance, each
//!   pushed through the ZTF, DES and Pan-STARRS detection models, with each
//!   survey's reach curve and the Rubin single-visit forecast (P13Map).

use rand::{Rng, SeedableRng};
use serde_json::{Value, json};

use p9_2016_cassini_ranging::perturbation::{
    INPOP_RESIDUAL_FLOOR_KM, differential_acceleration, favored_true_anomaly, prefit_amplitude,
    range_perturbation_amplitude,
};
use p9_2016_obliquity::parameter_survey::find_required_inclination;
use p9_2016_obliquity::secular_hamiltonian::{SecularParams, SpinOrbitState, integrate_obliquity};
use p9_2018_wise_search::detectability::max_detectable_distance as wise_reach;
use p9_2018_wise_search::survey_model::WiseSurvey;
use p9_2021_ztf::detection_efficiency::detection_probability_for_orbit as ztf_probability;
use p9_2021_ztf::survey_model::ZtfSurvey;
use p9_2022_des::color_models::fiducial as des_fiducial_colors;
use p9_2022_des::survey_model::DesSurvey;
use p9_2023_lsst_strategy::strategy::published::SINGLE_VISIT_DEPTH_R;
use p9_2024_panstarrs::detection_pipeline::detection_probability_for_orbit as ps1_probability;
use p9_2024_panstarrs::survey_model::Ps1Survey;
use p9_2025_iras_akari::survey_model::AkariFisSurvey;
use p9_core::analysis::photometry::{
    ALBEDO_NEPTUNE, BB21_ALBEDO_MAX, BB21_ALBEDO_MIN, SOLAR_V_MINUS_R, bb21_apparent_magnitude,
    mass_radius_neptunian,
};
use p9_core::analysis::surveys::limiting_magnitude;
use p9_core::analysis::thermal::{
    AU_M, C_LIGHT, effective_temp, flux_to_magnitude, max_detectable_distance, thermal_flux_jy,
};
use p9_core::constants::{
    DEG2RAD, EARTH_MASS_SOLAR, EARTH_RADIUS_KM, GYR_DAYS, RAD2DEG, YEAR_DAYS,
};
use p9_core::data::ephemeris_constraint::FAVORED_INTERVAL_DEG;
use p9_core::data::reference_population::{generate_reference_population, heliocentric_distance};
use p9_core::initial_conditions::giant_planets::{
    giant_planet_angular_momentum, p9_angular_momentum,
};
use p9_core::types::{OrbitalElements, P9Params, helio_distance_at_true_anomaly};

/// Internal-heat temperature floor for the thermal reaches (K), as in the
/// viability map and the thermal reproductions.
const INTERNAL_TEMP_K: f64 = 40.0;
/// Saturn's semimajor axis (AU), for the tidal-lever anchor.
const SATURN_AU: f64 = 9.58;

/// Cassini ranging: tidal acceleration, pre/post-fit range residual and
/// heliocentric distance around the BB16 orbit.
fn ranging() -> Value {
    let p = P9Params::nominal_2016();
    let au_per_day2_to_m_s2 = AU_M / (86_400.0 * 86_400.0);
    let nu_deg: Vec<f64> = (0..180).map(|k| k as f64 * 2.0).collect();
    let rows: Vec<(f64, f64, f64, f64)> = nu_deg
        .iter()
        .map(|&nd| {
            let nu = nd * DEG2RAD;
            (
                helio_distance_at_true_anomaly(&p, nu),
                prefit_amplitude(&p, nu),
                range_perturbation_amplitude(&p, nu),
                differential_acceleration(&p, nu).norm() * au_per_day2_to_m_s2,
            )
        })
        .collect();
    let favored = favored_true_anomaly(&p, 360);
    let d_fav = helio_distance_at_true_anomaly(&p, favored);
    json!({
        "mass_earth": p.mass_earth,
        "a": p.a,
        "e": p.e,
        "nu_deg": nu_deg,
        "distance_au": rows.iter().map(|r| r.0).collect::<Vec<_>>(),
        "prefit_km": rows.iter().map(|r| r.1).collect::<Vec<_>>(),
        "postfit_km": rows.iter().map(|r| r.2).collect::<Vec<_>>(),
        "tidal_m_s2": rows.iter().map(|r| r.3).collect::<Vec<_>>(),
        // The paper's Fig 6 zones: a pre-fit signal below the INPOP floor is
        // undetectable; a detectable one whose post-fit residual still
        // exceeds the floor would have spoiled the fit (excluded).
        "status": rows
            .iter()
            .map(|r| {
                if r.1 < INPOP_RESIDUAL_FLOOR_KM {
                    "undetectable"
                } else if r.2 > INPOP_RESIDUAL_FLOOR_KM {
                    "excluded"
                } else {
                    "allowed"
                }
            })
            .collect::<Vec<_>>(),
        "floor_km": INPOP_RESIDUAL_FLOOR_KM,
        "favored_nu_deg": favored * RAD2DEG,
        "favored_distance_au": d_fav,
        "favored_interval_deg": [FAVORED_INTERVAL_DEG.0, FAVORED_INTERVAL_DEG.1],
        // Direct pull on the Sun at the favoured position, for the "tidal
        // difference is tiny" comparison.
        "sun_pull_m_s2": p.gm() / (d_fav * d_fav) * au_per_day2_to_m_s2,
        "saturn_au": SATURN_AU,
    })
}

/// Solar obliquity from an inclined Planet Nine (Bailey et al. 2016).
fn obliquity() -> Value {
    let (m_earth, a9, e9) = (10.0, 500.0, 0.5);
    let m9 = m_earth * EARTH_MASS_SOLAR;
    let (i9_req, _) =
        find_required_inclination(m_earth, a9, e9, 6.0, 0.5).expect("a 6 deg solution exists");
    let l_gp = giant_planet_angular_momentum();
    let l_9 = p9_angular_momentum(m9, a9, e9);
    let init = SpinOrbitState::from_inclinations(i9_req, std::f64::consts::PI, l_gp, l_9);
    let params = SecularParams {
        m9_solar: m9,
        a9,
        e9,
        t_total: 4.5 * GYR_DAYS,
        dt: 5e4 * YEAR_DAYS,
    };
    let snaps = integrate_obliquity(init, &params, 2.5e7 * YEAR_DAYS);

    // Tips of the unit vectors seen down the total angular momentum (z).
    // The planets' plane and Planet Nine's orbit balance about z (the Sun's
    // spin carries negligible angular momentum), so the giant-planet normal
    // sits opposite Planet Nine's node at tilt asin((L9/Lgp) sin i9). The
    // Sun's spin tip is fixed by its node longitude and its angle to the
    // planets' normal (the obliquity); of the two solutions we follow the
    // one continuous with the start (spin aligned with the planets).
    let mut tips = Vec::with_capacity(snaps.len());
    let mut prev_theta = (l_9 / l_gp * i9_req.sin()).asin();
    for s in &snaps {
        let i_gp = (l_9 / l_gp * s.i_9.sin()).asin();
        let phi_gp = s.omega_big_9 + std::f64::consts::PI;
        let (a, b) = (i_gp.cos(), i_gp.sin() * (s.omega_big_sun - phi_gp).cos());
        let r = a.hypot(b);
        let delta = b.atan2(a);
        let spread = (s.obliquity.cos() / r).clamp(-1.0, 1.0).acos();
        let theta = [delta + spread, delta - spread]
            .into_iter()
            .map(|t| if t < 0.0 { -t } else { t })
            .min_by(|x, y| {
                (x - prev_theta)
                    .abs()
                    .partial_cmp(&(y - prev_theta).abs())
                    .unwrap()
            })
            .unwrap();
        prev_theta = theta;
        let tip = |tilt: f64, az: f64| [tilt.sin() * az.cos(), tilt.sin() * az.sin()];
        tips.push(json!({
            "t_gyr": s.t / GYR_DAYS,
            "obliquity_deg": s.obliquity * RAD2DEG,
            "p9": tip(s.i_9, s.omega_big_9),
            "planets": tip(i_gp, phi_gp),
            "sun": tip(theta, s.omega_big_sun),
        }));
    }
    json!({
        "mass_earth": m_earth,
        "a": a9,
        "e": e9,
        "i9_required_deg": i9_req * RAD2DEG,
        "observed_deg": 6.0,
        "snapshots": tips,
    })
}

/// Posterior draws vs survey reach in mass x current distance.
fn map() -> Value {
    const N: usize = 900;
    let mut rng = rand::rngs::StdRng::seed_from_u64(2024);
    let population = generate_reference_population(N, &mut rng);
    let ztf = ZtfSurvey::default();
    let des = DesSurvey::default();
    let des_colors = des_fiducial_colors();
    let ps1 = Ps1Survey::default();
    let rubin_v = SINGLE_VISIT_DEPTH_R + SOLAR_V_MINUS_R;

    let draws: Vec<Value> = population
        .iter()
        .map(|o| {
            let elements = OrbitalElements {
                a: o.a,
                e: o.e,
                i: o.i,
                omega: o.omega,
                omega_big: o.omega_big,
                mean_anomaly: o.mean_anomaly,
            };
            let r = heliocentric_distance(&P9Params {
                mass_earth: o.mass,
                a: o.a,
                e: o.e,
                i: o.i,
                omega: o.omega,
                omega_big: o.omega_big,
                mean_anomaly: o.mean_anomaly,
            });
            let mags = des_colors.band_magnitudes(o.mass, r);
            let p_ztf = ztf_probability(&ztf, &elements, o.v_magnitude);
            let p_des = des.detection_probability_for_orbit(&elements, &mags);
            let p_ps1 = ps1_probability(&ps1, &elements, o.v_magnitude);
            json!({
                "mass": o.mass,
                "r_au": r,
                "v": o.v_magnitude,
                // cumulative detection probability after ZTF, +DES, +PS1
                "p_cum": [
                    p_ztf,
                    1.0 - (1.0 - p_ztf) * (1.0 - p_des),
                    1.0 - (1.0 - p_ztf) * (1.0 - p_des) * (1.0 - p_ps1),
                ],
                // fixed per-draw uniform: a draw is drawn as excluded once its
                // cumulative detection probability exceeds it
                "u": rng.gen_range(0.0..1.0_f64),
            })
        })
        .collect();
    let p_cum = |d: &Value, k: usize| d["p_cum"][k].as_f64().unwrap();
    let excluded_frac: Vec<f64> = (0..3)
        .map(|k| draws.iter().map(|d| p_cum(d, k)).sum::<f64>() / N as f64)
        .collect();
    // Of what survives all three surveys, the share bright enough for a
    // single Rubin visit (probability-weighted).
    let survive = |d: &Value| 1.0 - p_cum(d, 2);
    let rubin_bright_frac = draws
        .iter()
        .filter(|d| d["v"].as_f64().unwrap() <= rubin_v)
        .map(survive)
        .sum::<f64>()
        / draws.iter().map(survive).sum::<f64>();

    let mass_grid: Vec<f64> = (0..=76).map(|k| 1.0 + k as f64 * 0.25).collect();
    let albedo_mid = 0.5 * (BB21_ALBEDO_MIN + BB21_ALBEDO_MAX);
    let optical = |v_lim: f64| -> Vec<Option<f64>> {
        mass_grid
            .iter()
            .map(|&m| {
                max_detectable_distance(60.0, 20_000.0, v_lim, |d| {
                    bb21_apparent_magnitude(m, albedo_mid, d)
                })
            })
            .collect()
    };
    let ztf_v = limiting_magnitude("ZTF").expect("ZTF depth");
    let ps1_v = limiting_magnitude("PS1 3pi").expect("PS1 depth") + SOLAR_V_MINUS_R;
    let des_v = limiting_magnitude("DES").expect("DES depth") + SOLAR_V_MINUS_R;

    // Best far-infrared all-sky reach: WISE W1, IRAS 60 um, AKARI 90 um.
    let wise = WiseSurvey::default();
    let akari_jy = AkariFisSurvey::default().sensitivity_jy;
    let thermal_reach = |m: f64, nu_hz: f64, limit_jy: f64| -> Option<f64> {
        let r_m = mass_radius_neptunian(m) * EARTH_RADIUS_KM * 1e3;
        max_detectable_distance(60.0, 20_000.0, flux_to_magnitude(limit_jy, 1.0), |d| {
            let t = effective_temp(d, ALBEDO_NEPTUNE, INTERNAL_TEMP_K);
            flux_to_magnitude(thermal_flux_jy(t, r_m, d, nu_hz), 1.0)
        })
    };
    let infrared: Vec<Option<f64>> = mass_grid
        .iter()
        .map(|&m| {
            [
                wise_reach(&wise, m),
                thermal_reach(m, C_LIGHT / 60e-6, 0.2),
                thermal_reach(m, C_LIGHT / 90e-6, akari_jy),
            ]
            .into_iter()
            .flatten()
            .reduce(f64::max)
        })
        .collect();

    json!({
        "n": N,
        "draws": draws,
        "mass_grid": mass_grid,
        "albedo_mid": albedo_mid,
        "reach": {
            "ZTF": optical(ztf_v),
            "DES": optical(des_v),
            "Pan-STARRS": optical(ps1_v),
            "Rubin": optical(rubin_v),
            "infrared": infrared,
        },
        "rubin_v_limit": rubin_v,
        "excluded_frac": excluded_frac,
        "rubin_bright_frac": rubin_bright_frac,
    })
}

pub fn export() -> Value {
    json!({
        "ranging": ranging(),
        "obliquity": obliquity(),
        "map": map(),
    })
}
