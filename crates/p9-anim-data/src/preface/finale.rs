//! Film export for the `scenes/finale/` scenes: the numbers its scenes draw.
//!
//! One seeded Brown & Batygin (2021) reference population (the same draw the
//! `exclusion` section and the Pan-STARRS1 paper export use) is placed on
//! tonight's sky and pushed through the reproduced ZTF, DES and PS1 searches
//! and the Rubin/LSST baseline, so the finale can carve the predicted sky
//! away survey by survey and say where, and how faint, the survivors are.

use rand::SeedableRng;
use serde_json::{Value, json};

use p9_2021_ztf::detection_efficiency::detection_probability_for_orbit as ztf_probability;
use p9_2021_ztf::survey_model::ZtfSurvey;
use p9_2022_des::color_models::fiducial;
use p9_2022_des::survey_model::DesSurvey;
use p9_2023_lsst_strategy::{LsstStrategy, P9Sample, discovery_probability};
use p9_2024_panstarrs::detection_pipeline::detection_probability_for_orbit as ps1_probability;
use p9_2024_panstarrs::survey_model::Ps1Survey;
use p9_core::analysis::photometry::{ALBEDO_NEPTUNE, planet_apparent_magnitude};
use p9_core::constants::{A_NEPTUNE_AU, DEG2RAD, GM_SUN, RAD2DEG};
use p9_core::coords::sky::{ecliptic_vec_to_equatorial_deg, equatorial_to_galactic};
use p9_core::data::reference_population::{generate_reference_population, heliocentric_distance};
use p9_core::types::{OrbitalElements, P9Params, solve_kepler};

const N_ORBITS: usize = 3000;
const SEED: u64 = 2024;
/// Galactic latitude inside which the ground surveys' fields are crowded.
const PLANE_B_DEG: f64 = 10.0;

struct Member {
    ra: f64,
    dec: f64,
    b: f64,
    dist: f64,
    a: f64,
    v: f64,
    p: [f64; 3],
    p_rubin: f64,
}

impl Member {
    fn survival(&self) -> f64 {
        self.p.iter().map(|p| 1.0 - p).product()
    }
}

fn population() -> Vec<Member> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(SEED);
    let ztf = ZtfSurvey::default();
    let des = DesSurvey::default();
    let colours = fiducial();
    let ps1 = Ps1Survey::default();
    let rubin = LsstStrategy::baseline();
    generate_reference_population(N_ORBITS, &mut rng)
        .into_iter()
        .map(|o| {
            let elements = OrbitalElements {
                a: o.a,
                e: o.e,
                i: o.i,
                omega: o.omega,
                omega_big: o.omega_big,
                mean_anomaly: o.mean_anomaly,
            };
            let dist = heliocentric_distance(&P9Params {
                mass_earth: o.mass,
                a: o.a,
                e: o.e,
                i: o.i,
                omega: o.omega,
                omega_big: o.omega_big,
                mean_anomaly: o.mean_anomaly,
            });
            let (ra, dec) = ecliptic_vec_to_equatorial_deg(&elements.to_state_vector(GM_SUN).pos);
            let (_, b) = equatorial_to_galactic(ra * DEG2RAD, dec * DEG2RAD);
            let mags = colours.band_magnitudes(o.mass, dist);
            let sample = P9Sample {
                elements,
                v_magnitude: o.v_magnitude,
            };
            Member {
                ra,
                dec,
                b: b * RAD2DEG,
                dist,
                a: o.a,
                v: o.v_magnitude,
                p: [
                    ztf_probability(&ztf, &elements, o.v_magnitude),
                    des.detection_probability_for_orbit(&elements, &mags),
                    ps1_probability(&ps1, &elements, o.v_magnitude),
                ],
                p_rubin: discovery_probability(&rubin, &sample),
            }
        })
        .collect()
}

fn round(x: f64, digits: i32) -> f64 {
    let s = 10f64.powi(digits);
    (x * s).round() / s
}

/// Weighted share of `pop` for which `pick` holds.
fn share(pop: &[Member], weight: impl Fn(&Member) -> f64, pick: impl Fn(&Member) -> bool) -> f64 {
    let total: f64 = pop.iter().map(&weight).sum();
    pop.iter().filter(|m| pick(m)).map(&weight).sum::<f64>() / total
}

/// Weighted median of `value` over `pop`.
fn median(pop: &[Member], weight: impl Fn(&Member) -> f64, value: impl Fn(&Member) -> f64) -> f64 {
    let mut rows: Vec<(f64, f64)> = pop.iter().map(|m| (value(m), weight(m))).collect();
    rows.sort_by(|x, y| x.0.total_cmp(&y.0));
    let half = rows.iter().map(|r| r.1).sum::<f64>() / 2.0;
    let mut acc = 0.0;
    for (v, w) in rows {
        acc += w;
        if acc >= half {
            return v;
        }
    }
    f64::NAN
}

/// Histogram of `value` with each member weighted by `weight / N`, so bars
/// are fractions of the whole prediction.
fn histogram(
    pop: &[Member],
    edges: &[f64],
    weight: impl Fn(&Member) -> f64,
    value: impl Fn(&Member) -> f64,
) -> Vec<f64> {
    let mut counts = vec![0.0; edges.len() - 1];
    for m in pop {
        let x = value(m);
        if let Some(k) = edges.windows(2).position(|w| x >= w[0] && x < w[1]) {
            counts[k] += weight(m) / pop.len() as f64;
        }
    }
    counts.into_iter().map(|c| round(c, 5)).collect()
}

/// Observed opposition magnitude of Pluto, the familiar faint anchor, read
/// from the film's shared solar-system reference table.
fn pluto_v() -> f64 {
    crate::solar_system()["magnitudes"]
        .as_array()
        .and_then(|rows| rows.iter().find(|r| r["name"] == "Pluto"))
        .and_then(|r| r["app_mag"].as_f64())
        .expect("Pluto in the solar-system magnitude table")
}

/// DES footprint as 2 x 2 degree cells (RA, Dec centres) for the sky map.
fn des_cells() -> Vec<[f64; 2]> {
    let des = DesSurvey::default();
    let mut cells = Vec::new();
    for j in 0..45 {
        let dec = -89.0 + 2.0 * j as f64;
        for k in 0..180 {
            let ra = 1.0 + 2.0 * k as f64;
            if des.is_in_footprint(ra, dec) {
                cells.push([ra, dec]);
            }
        }
    }
    cells
}

/// The Brown & Batygin (2021) best-fit orbit as a clock: positions at equal
/// time steps in the orbital plane (perihelion along +x), the share of the
/// period spent beyond the semi-major axis, and how much fainter aphelion is
/// than perihelion in reflected light.
fn orbit_clock() -> Value {
    let p = P9Params::mcmc_2021();
    let steps = 24;
    let points: Vec<[f64; 2]> = (0..steps)
        .map(|k| {
            let m = std::f64::consts::TAU * k as f64 / steps as f64;
            let ea = solve_kepler(p.e, m);
            [
                round(p.a * (ea.cos() - p.e), 2),
                round(p.a * (1.0 - p.e * p.e).sqrt() * ea.sin(), 2),
            ]
        })
        .collect();
    let fine = 3600;
    let beyond = (0..fine)
        .filter(|&k| {
            let m = std::f64::consts::TAU * (k as f64 + 0.5) / fine as f64;
            p.a * (1.0 - p.e * solve_kepler(p.e, m).cos()) > p.a
        })
        .count();
    let (q, big_q) = (p.a * (1.0 - p.e), p.a * (1.0 + p.e));
    let v = |r: f64| planet_apparent_magnitude(p.mass_earth, ALBEDO_NEPTUNE, r);
    json!({
        "a": p.a,
        "e": p.e,
        "q": q,
        "big_q": big_q,
        "period_yr": p.a.powf(1.5),
        "points": points,
        "time_beyond_a": beyond as f64 / fine as f64,
        "v_peri": v(q),
        "v_aph": v(big_q),
    })
}

pub fn export() -> Value {
    let pop = population();
    let n = pop.len() as f64;

    // Per-orbit OR in the papers' ZTF -> DES -> PS1 order.
    let mut cumulative = Vec::new();
    for k in 1..=3 {
        let excluded: f64 = pop
            .iter()
            .map(|m| 1.0 - m.p[..k].iter().map(|p| 1.0 - p).product::<f64>())
            .sum();
        cumulative.push(excluded / n);
    }

    let prior = |_: &Member| 1.0;
    let surv = |m: &Member| m.survival();
    let in_plane = |m: &Member| m.b.abs() < PLANE_B_DEG;
    let rubin = LsstStrategy::baseline();
    let north_of_rubin = |m: &Member| m.dec > rubin.dec_max_deg;
    let beyond_a = |m: &Member| m.dist > m.a;

    let survivors: f64 = pop.iter().map(surv).sum();
    let unseen = |m: &Member| m.survival() * (1.0 - m.p_rubin);
    let rubin_of_survivors = pop.iter().map(|m| m.survival() * m.p_rubin).sum::<f64>() / survivors;
    let v_med = median(&pop, surv, |m| m.v);
    let d_med = median(&pop, surv, |m| m.dist);
    let pluto = pluto_v();

    let v_edges: Vec<f64> = (0..=22).map(|k| 15.0 + 0.5 * k as f64).collect();
    let d_edges: Vec<f64> = (0..=24).map(|k| 50.0 * k as f64).collect();

    let ztf = ZtfSurvey::default();
    let des = DesSurvey::default();
    let ps1 = Ps1Survey::default();

    json!({
        "n": pop.len(),
        "seed": SEED,
        "orbits": {
            "ra": pop.iter().map(|m| round(m.ra, 2)).collect::<Vec<_>>(),
            "dec": pop.iter().map(|m| round(m.dec, 2)).collect::<Vec<_>>(),
            "b": pop.iter().map(|m| round(m.b, 2)).collect::<Vec<_>>(),
            "dist": pop.iter().map(|m| round(m.dist, 1)).collect::<Vec<_>>(),
            "v": pop.iter().map(|m| round(m.v, 3)).collect::<Vec<_>>(),
            "p_ztf": pop.iter().map(|m| round(m.p[0], 4)).collect::<Vec<_>>(),
            "p_des": pop.iter().map(|m| round(m.p[1], 4)).collect::<Vec<_>>(),
            "p_ps1": pop.iter().map(|m| round(m.p[2], 4)).collect::<Vec<_>>(),
            "p_rubin": pop.iter().map(|m| round(m.p_rubin, 4)).collect::<Vec<_>>(),
        },
        "surveys": [
            {"key": "p_ztf", "name": "ZTF", "year": 2021, "depth_v": ztf.depth_limit,
             "dec_min": ztf.dec_limit_deg, "cumulative": cumulative[0]},
            {"key": "p_des", "name": "DES", "year": 2022, "depth_r": des.depth_r,
             "area_deg2": des.footprint_area, "cumulative": cumulative[1]},
            {"key": "p_ps1", "name": "Pan-STARRS1", "year": 2024, "depth_v": ps1.depth_limit,
             "dec_min": ps1.dec_limit_deg, "cumulative": cumulative[2]},
        ],
        "des_cells": des_cells(),
        "hiding": {
            "remaining": 1.0 - cumulative[2],
            "plane_b_deg": PLANE_B_DEG,
            "plane_share_prior": share(&pop, prior, in_plane),
            "plane_share_survivors": share(&pop, surv, in_plane),
            "north_share_prior": share(&pop, prior, north_of_rubin),
            "north_share_survivors": share(&pop, surv, north_of_rubin),
            "beyond_a_share_prior": share(&pop, prior, beyond_a),
            "beyond_a_share_survivors": share(&pop, surv, beyond_a),
            "v_median_prior": median(&pop, prior, |m| m.v),
            "v_median_survivors": v_med,
            "dist_median_prior": median(&pop, prior, |m| m.dist),
            "dist_median_survivors": d_med,
            "pluto_v": pluto,
            "fainter_than_pluto": 10f64.powf(0.4 * (v_med - pluto)),
            "neptune_au": A_NEPTUNE_AU,
            "times_neptune": d_med / A_NEPTUNE_AU,
            "v_edges": v_edges,
            "v_prior": histogram(&pop, &v_edges, prior, |m| m.v),
            "v_survivors": histogram(&pop, &v_edges, surv, |m| m.v),
            "dist_edges": d_edges,
            "dist_prior": histogram(&pop, &d_edges, prior, |m| m.dist),
            "dist_survivors": histogram(&pop, &d_edges, surv, |m| m.dist),
        },
        "orbit_clock": orbit_clock(),
        "rubin": {
            "left_plane_share": share(&pop, unseen, in_plane),
            "left_north_share": share(&pop, unseen, north_of_rubin),
            "left_either_share": share(&pop, unseen, |m| in_plane(m) || north_of_rubin(m)),
            "dec_min": rubin.dec_min_deg,
            "dec_max": rubin.dec_max_deg,
            "gal_b_min": rubin.galactic_lat_min_deg,
            "depth_r": rubin.single_visit_depth,
            "of_survivors": rubin_of_survivors,
            "of_prediction": rubin_of_survivors * survivors / n,
            "left_after": (survivors - pop.iter().map(|m| m.survival() * m.p_rubin).sum::<f64>()) / n,
        },
    })
}
