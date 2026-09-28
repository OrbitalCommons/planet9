//! Film export for `p9-2021-des-catalog`: the numbers its scene and ledger entry draw.

use p9_2021_des_catalog::catalog::{Angle, DES_EXTREME_TNOS};
use p9_2021_des_catalog::clustering::{
    ClusteringResult, NullModel, clustering_battery, ecliptic_to_equatorial, perihelion_direction,
};
use p9_2021_des_catalog::completeness::{
    apparent_mag, des_depth_r, magnitude_completeness, sky_coverage_fraction,
};
use p9_2021_des_catalog::reference::{DEPTH_R, FOOTPRINT_DEG2, N_EXTREME_TNOS, N_NEW_TNOS, N_TNOS};
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::analysis::surveys::{
    DES_FOOTPRINT_BANDS, des_footprint_contains, des_footprint_solid_angle_deg2,
};
use p9_core::constants::{DEG2RAD, RAD2DEG};
use serde_json::{Value, json};

const SEED: u64 = 2021;
const MC_ITERS: usize = 30_000;
const COVERAGE_SAMPLES: usize = 400_000;
/// Galactic-latitude mask used for the coverage fraction (degrees).
const B_CUT_DEG: f64 = 10.0;

fn battery_json(results: &[ClusteringResult]) -> Vec<Value> {
    results
        .iter()
        .map(|r| {
            json!({
                "angle": r.angle.label(),
                "r_bar": r.r_bar,
                "rayleigh_p": r.rayleigh_p,
                "mc_p": r.mc_p,
            })
        })
        .collect()
}

fn varpi_p(results: &[ClusteringResult]) -> f64 {
    results
        .iter()
        .find(|r| r.angle == Angle::Varpi)
        .map(|r| r.mc_p)
        .unwrap()
}

pub fn export() -> Value {
    let objects: Vec<Value> = DES_EXTREME_TNOS
        .iter()
        .map(|o| {
            let (lambda, beta) =
                perihelion_direction(o.i_deg * DEG2RAD, Angle::ArgPeri.of(o), Angle::Node.of(o));
            let (ra, dec) = ecliptic_to_equatorial(lambda, beta);
            let (ra_deg, dec_deg) = (ra * RAD2DEG, dec * RAD2DEG);
            let q = o.a * (1.0 - o.e);
            json!({
                "name": o.name,
                "a": o.a,
                "e": o.e,
                "q": q,
                "i_deg": o.i_deg,
                "h_mag": o.h_mag,
                "varpi_deg": Angle::Varpi.of(o) * RAD2DEG,
                "node_deg": Angle::Node.of(o) * RAD2DEG,
                "ra_deg": ra_deg,
                "dec_deg": dec_deg,
                "in_footprint": des_footprint_contains(ra_deg, dec_deg),
                "r_at_perihelion": apparent_mag(o.h_mag, q),
            })
        })
        .collect();

    let footprint: Vec<Value> = DES_FOOTPRINT_BANDS
        .iter()
        .map(|b| {
            json!({
                "ra_start_deg": b.ra_start_deg,
                "ra_end_deg": b.ra_end_deg,
                "dec_min_deg": b.dec_min_deg,
                "dec_max_deg": b.dec_max_deg,
            })
        })
        .collect();

    // Detection efficiency against apparent magnitude.
    let r_mag: Vec<f64> = (0..=80).map(|k| 19.0 + 0.1 * k as f64).collect();
    let efficiency: Vec<f64> = r_mag.iter().map(|&m| magnitude_completeness(m)).collect();

    let flat = clustering_battery(&DES_EXTREME_TNOS, NullModel::Uniform, SEED, MC_ITERS);
    let sel = clustering_battery(&DES_EXTREME_TNOS, NullModel::DesSelection, SEED, MC_ITERS);
    let p_flat = varpi_p(&flat);
    let p_sel = varpi_p(&sel);

    // The strongest of the three angle tests, read with the number of tests in
    // mind: the chance that the smallest of `n` independent p-values is at
    // least this small.
    let min_p = |r: &[ClusteringResult]| r.iter().map(|c| c.mc_p).fold(1.0_f64, f64::min);
    let after_trials = |p: f64| 1.0 - (1.0 - p).powi(sel.len() as i32);

    json!({
        "n_tnos": N_TNOS,
        "n_new_tnos": N_NEW_TNOS,
        "n_extreme": N_EXTREME_TNOS,
        "n_extreme_encoded": DES_EXTREME_TNOS.len(),
        "paper_footprint_deg2": FOOTPRINT_DEG2,
        "footprint_deg2": des_footprint_solid_angle_deg2(),
        "paper_depth_r": DEPTH_R,
        "depth_r": des_depth_r(),
        "sky_fraction": sky_coverage_fraction(B_CUT_DEG, SEED, COVERAGE_SAMPLES),
        "b_cut_deg": B_CUT_DEG,
        "objects": objects,
        "footprint": footprint,
        "efficiency": {"r_mag": r_mag, "fraction": efficiency},
        "battery_flat": battery_json(&flat),
        "battery_selection": battery_json(&sel),
        "p_varpi_flat": p_flat,
        "p_varpi_selection": p_sel,
        "min_p_flat": min_p(&flat),
        "min_p_selection": min_p(&sel),
        "p_after_trials": after_trials(min_p(&sel)),
        "sigma_flat": p_value_to_sigma(after_trials(min_p(&flat))),
        "sigma": p_value_to_sigma(after_trials(min_p(&sel))),
    })
}
