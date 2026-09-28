//! Film export for `p9-2024-siraj-orbit`: the numbers its scene and ledger entry draw.

use p9_2024_siraj_orbit::confinement::kappa_from_angles;
use p9_2024_siraj_orbit::forcing::mass_for_strength;
use p9_2024_siraj_orbit::posterior::{KAPPA_REFERENCE, OrbitPoint, calibrated_posterior};
use p9_2024_siraj_orbit::{
    BB2021_A_AU, BB2021_A_SIGMA_AU, BB2021_I_DEG, BB2021_MASS_EARTH, BB2021_MASS_SIGMA_EARTH,
    infer_from_etnos,
};
use p9_core::analysis::circular::mean_resultant_length;
use p9_core::data::ephemeris_constraint::{
    SIRAJ_2024_A_AU, SIRAJ_2024_A_SIGMA_AU, SIRAJ_2024_I_DEG, SIRAJ_2024_MASS_EARTH,
    SIRAJ_2024_MASS_SIGMA_EARTH,
};
use p9_core::data::etno::longitudes_of_perihelion;
use p9_core::types::P9Params;
use serde_json::{Value, json};

/// Grid of the posterior map over (semi-major axis, mass).
const A_RANGE: (f64, f64) = (150.0, 650.0);
const M_RANGE: (f64, f64) = (0.0, 14.0);
const N_A: usize = 50;
const N_M: usize = 35;
/// Published best-fit eccentricity (the crate models only the mass-distance
/// plane, so the abstract's value is carried as a reference).
const PAPER_E: f64 = 0.29;
/// Published share of the Brown & Batygin (2021) reference population that
/// lands within 1 sigma of this paper's maximum.
const PAPER_BB21_OVERLAP: f64 = 0.0006;

pub fn export() -> Value {
    let varpis = longitudes_of_perihelion();
    let post = calibrated_posterior(&varpis);
    let map = infer_from_etnos(&varpis);

    let a_centres: Vec<f64> = (0..N_A)
        .map(|k| A_RANGE.0 + (A_RANGE.1 - A_RANGE.0) * (k as f64 + 0.5) / N_A as f64)
        .collect();
    let m_centres: Vec<f64> = (0..N_M)
        .map(|k| M_RANGE.0 + (M_RANGE.1 - M_RANGE.0) * (k as f64 + 0.5) / N_M as f64)
        .collect();
    let best = post.neg2_log_post(map);
    // Relative posterior density, row-major (mass rows x semi-major-axis columns).
    let density: Vec<f64> = m_centres
        .iter()
        .flat_map(|&mass_earth| {
            let post = &post;
            a_centres.iter().map(move |&a_p| {
                (-0.5 * (post.neg2_log_post(OrbitPoint { mass_earth, a_p }) - best)).exp()
            })
        })
        .collect();

    // The confinement ridge m = S a^3: what the observed sample implies, and
    // what a sample as weakly clustered as the paper's would imply.
    let ridge_a: Vec<f64> = (0..=100).map(|k| 150.0 + 5.0 * k as f64).collect();
    let ridge_m: Vec<f64> = ridge_a
        .iter()
        .map(|&a| mass_for_strength(post.s_obs, a))
        .collect();
    let kappa = kappa_from_angles(&varpis);
    let s_paper = post.s_obs * KAPPA_REFERENCE / kappa;
    let ridge_m_paper: Vec<f64> = ridge_a
        .iter()
        .map(|&a| mass_for_strength(s_paper, a))
        .collect();

    let bb = P9Params::mcmc_2021();

    json!({
        "n_sample": varpis.len(),
        "r_bar": mean_resultant_length(&varpis),
        "kappa": kappa,
        "kappa_reference": KAPPA_REFERENCE,
        "a_prior": post.a_ref,
        "a_prior_sigma": post.sigma_a,
        "map_mass": map.mass_earth,
        "map_a": map.a_p,
        "grid": {
            "a": a_centres,
            "mass": m_centres,
            "density": density,
        },
        "ridge": {"a": ridge_a, "mass": ridge_m, "mass_paper_sample": ridge_m_paper},
        "mass": SIRAJ_2024_MASS_EARTH,
        "mass_sigma": SIRAJ_2024_MASS_SIGMA_EARTH,
        "a": SIRAJ_2024_A_AU,
        "a_sigma": SIRAJ_2024_A_SIGMA_AU,
        "e": PAPER_E,
        "i": SIRAJ_2024_I_DEG,
        "paper_bb21_overlap": PAPER_BB21_OVERLAP,
        "bb21": {
            "mass": BB2021_MASS_EARTH,
            "mass_sigma": BB2021_MASS_SIGMA_EARTH,
            "a": BB2021_A_AU,
            "a_sigma": BB2021_A_SIGMA_AU,
            "e": bb.e,
            "i": BB2021_I_DEG,
        },
    })
}
