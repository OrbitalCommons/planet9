//! Film export for `p9-2025-stellar-flybys`: the numbers its scene and ledger entry draw.

use p9_2025_stellar_flybys::encounter::{
    CLUSTER_DENSITY_PER_PC3, CLUSTER_RESIDENCE_TIME_YR, CLUSTER_VELOCITY_DISPERSION_KMS,
    PERTURBER_MASS_SOLAR, Q_STAR_MAX_AU, expected_close_encounters, focused_cross_section_au2,
    p_close_encounter,
};
use p9_2025_stellar_flybys::geometry::{
    COPLANAR_TOL_RAD, SYMMETRIC_TOL_RAD, coplanar_fraction, geometry_fraction, symmetric_fraction,
};
use p9_2025_stellar_flybys::{
    F_SUCCESS, INCLINATION_CONSTRAINT_DEG, PUBLISHED_MAX_PROBABILITY, total_probability,
};
use p9_core::constants::{KMS_TO_AUDAY, PC_AU, YEAR_DAYS};
use rand::{Rng, SeedableRng};
use serde_json::{Value, json};
use std::f64::consts::PI;

/// Isotropic encounter orientations drawn for the acceptance map.
const N_ORIENTATIONS: usize = 1500;

/// Probability of at least one encounter inside `q_au` for a cluster of the
/// given density, with the crate's focused cross-section, velocity dispersion
/// and residence time.
fn p_encounter_within(q_au: f64, density_per_pc3: f64) -> f64 {
    let sigma = focused_cross_section_au2(
        q_au,
        CLUSTER_VELOCITY_DISPERSION_KMS,
        1.0 + PERTURBER_MASS_SOLAR,
    );
    let rate =
        density_per_pc3 / PC_AU.powi(3) * sigma * CLUSTER_VELOCITY_DISPERSION_KMS * KMS_TO_AUDAY;
    1.0 - (-rate * CLUSTER_RESIDENCE_TIME_YR * YEAR_DAYS).exp()
}

/// Which acceptance band an orientation (inclination, argument of periastron)
/// falls in, by the crate's tolerances.
fn band(incl: f64, arg_peri: f64) -> &'static str {
    let coplanar = !(COPLANAR_TOL_RAD..=PI - COPLANAR_TOL_RAD).contains(&incl);
    let near_node = |w: f64| {
        let d = w.rem_euclid(PI);
        d.min(PI - d) < SYMMETRIC_TOL_RAD
    };
    let symmetric = incl.cos().abs() < SYMMETRIC_TOL_RAD.sin() && near_node(arg_peri);
    if coplanar {
        "coplanar"
    } else if symmetric {
        "symmetric"
    } else {
        "rejected"
    }
}

pub fn export() -> Value {
    let mut rng = rand::rngs::StdRng::seed_from_u64(2025_0516);
    let orientations: Vec<Value> = (0..N_ORIENTATIONS)
        .map(|_| {
            let cos_i: f64 = rng.gen_range(-1.0..1.0);
            let arg_peri: f64 = rng.gen_range(0.0..2.0 * PI);
            json!({
                "cos_i": cos_i,
                "arg_peri_deg": arg_peri.to_degrees(),
                "band": band(cos_i.acos(), arg_peri),
            })
        })
        .collect();
    let accepted = orientations
        .iter()
        .filter(|o| o["band"] != "rejected")
        .count();

    let q_grid: Vec<f64> = (0..=58).map(|k| 100.0 + 50.0 * k as f64).collect();
    let densities = [30.0, CLUSTER_DENSITY_PER_PC3, 300.0];
    let p_vs_q: Vec<Value> = densities
        .iter()
        .map(|&n| {
            let p: Vec<f64> = q_grid.iter().map(|&q| p_encounter_within(q, n)).collect();
            json!({"density_per_pc3": n, "p": p})
        })
        .collect();

    let p_enc = p_close_encounter();
    let f_geom = geometry_fraction();
    let total = total_probability();
    let geometric = PI * Q_STAR_MAX_AU * Q_STAR_MAX_AU;

    json!({
        "q_star_max_au": Q_STAR_MAX_AU,
        "inclination_limit_deg": INCLINATION_CONSTRAINT_DEG,
        "cluster": {
            "density_per_pc3": CLUSTER_DENSITY_PER_PC3,
            "velocity_dispersion_kms": CLUSTER_VELOCITY_DISPERSION_KMS,
            "residence_myr": CLUSTER_RESIDENCE_TIME_YR / 1e6,
            "perturber_mass_solar": PERTURBER_MASS_SOLAR,
        },
        "focusing_factor": focused_cross_section_au2(
            Q_STAR_MAX_AU,
            CLUSTER_VELOCITY_DISPERSION_KMS,
            1.0 + PERTURBER_MASS_SOLAR,
        ) / geometric,
        "expected_encounters": expected_close_encounters(),
        "p_encounter": p_enc,
        "coplanar_tol_deg": COPLANAR_TOL_RAD.to_degrees(),
        "symmetric_tol_deg": SYMMETRIC_TOL_RAD.to_degrees(),
        "f_coplanar": coplanar_fraction(COPLANAR_TOL_RAD),
        "f_symmetric": symmetric_fraction(SYMMETRIC_TOL_RAD),
        "f_geometry": f_geom,
        "f_success": F_SUCCESS,
        "p_after_geometry": p_enc * f_geom,
        "p_total": total,
        "published_max": PUBLISHED_MAX_PROBABILITY,
        "orientations": orientations,
        "n_orientations": N_ORIENTATIONS,
        "n_accepted": accepted,
        "sampled_fraction": accepted as f64 / N_ORIENTATIONS as f64,
        "q_grid_au": q_grid,
        "p_vs_q": p_vs_q,
    })
}
