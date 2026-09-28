//! Film export for `p9-2024-panstarrs`: the numbers its scene and ledger entry draw.

use p9_2021_ztf::detection_efficiency::detection_probability_for_orbit as ztf_probability;
use p9_2021_ztf::survey_model::ZtfSurvey;
use p9_2022_des::color_models::fiducial;
use p9_2022_des::survey_model::DesSurvey;
use p9_2024_panstarrs::combined_exclusion::{
    CombinedExclusion, UpdatedParameters, compute_combined_from_population,
};
use p9_2024_panstarrs::detection_pipeline::{
    Ps1Result, detection_probability_for_orbit, equatorial_position_deg,
};
use p9_2024_panstarrs::survey_model::Ps1Survey;
use p9_core::analysis::surveys::DES_FOOTPRINT_BANDS;
use p9_core::data::reference_population::{generate_reference_population, heliocentric_distance};
use p9_core::types::{OrbitalElements, P9Params};
use rand::SeedableRng;
use serde_json::{Value, json};

/// Size and seed of the scored population (the same draw the exclusion
/// fractions are computed from).
const N_POP: usize = 3000;
const SEED: u64 = 2024;
/// Members drawn on the sky panel.
const N_SKY: usize = 800;

/// Median of `values` under `weights`.
fn weighted_median(values: &[f64], weights: &[f64]) -> f64 {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&i, &j| values[i].total_cmp(&values[j]));
    let half = 0.5 * weights.iter().sum::<f64>();
    let mut acc = 0.0;
    for k in order {
        acc += weights[k];
        if acc >= half {
            return values[k];
        }
    }
    f64::NAN
}

pub fn export() -> Value {
    let ex = compute_combined_from_population(N_POP, SEED);
    let paper = CombinedExclusion::paper_values();
    let paper_counts = Ps1Result::paper_values();
    let updated = UpdatedParameters::paper_values();
    let orbit = updated.to_p9_params();

    let ztf = ZtfSurvey::default();
    let des = DesSurvey::default();
    let des_color = fiducial();
    let ps1 = Ps1Survey::default();

    // The same seeded Brown & Batygin (2021) population the fractions above
    // are averaged over, scored survey by survey.
    let mut rng = rand::rngs::StdRng::seed_from_u64(SEED);
    let population = generate_reference_population(N_POP, &mut rng);

    let mut sky = Vec::with_capacity(N_SKY);
    let (mut a, mut mass, mut v_mag, mut survive, mut all) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut ps1_total = 0.0;
    for (k, obj) in population.iter().enumerate() {
        let elements = OrbitalElements {
            a: obj.a,
            e: obj.e,
            i: obj.i,
            omega: obj.omega,
            omega_big: obj.omega_big,
            mean_anomaly: obj.mean_anomaly,
        };
        let r_au = heliocentric_distance(&P9Params {
            mass_earth: obj.mass,
            a: obj.a,
            e: obj.e,
            i: obj.i,
            omega: obj.omega,
            omega_big: obj.omega_big,
            mean_anomaly: obj.mean_anomaly,
        });
        let p_ztf = ztf_probability(&ztf, &elements, obj.v_magnitude);
        let p_des = des
            .detection_probability_for_orbit(&elements, &des_color.band_magnitudes(obj.mass, r_au));
        let p_ps1 = detection_probability_for_orbit(&ps1, &elements, obj.v_magnitude);
        let p_before = 1.0 - (1.0 - p_ztf) * (1.0 - p_des);
        let p_after = 1.0 - (1.0 - p_before) * (1.0 - p_ps1);
        ps1_total += p_ps1;

        a.push(obj.a);
        mass.push(obj.mass);
        v_mag.push(obj.v_magnitude);
        survive.push(1.0 - p_after);
        all.push(1.0);

        if k < N_SKY {
            let (ra, dec) = equatorial_position_deg(&elements);
            sky.push(json!({
                "ra_deg": ra,
                "dec_deg": dec,
                "dist_au": r_au,
                "v_mag": obj.v_magnitude,
                "p_before": p_before,
                "p_ps1": p_ps1,
                "p_after": p_after,
            }));
        }
    }

    let mags: Vec<f64> = (0..=60).map(|k| 17.0 + 0.1 * k as f64).collect();
    let ps1_eff: Vec<f64> = mags
        .iter()
        .map(|&m| ps1.detection_probability(m, 10.0))
        .collect();
    let ztf_eff: Vec<f64> = mags.iter().map(|&m| ztf.detection_efficiency(m)).collect();

    let des_boxes: Vec<Value> = DES_FOOTPRINT_BANDS
        .iter()
        .map(|b| {
            json!({
                "ra_lo": b.ra_start_deg,
                "ra_hi": b.ra_end_deg,
                "dec_lo": b.dec_min_deg,
                "dec_hi": b.dec_max_deg,
            })
        })
        .collect();

    json!({
        "ztf": ex.ztf_frac,
        "des_unique": ex.des_unique,
        "ps1_unique": ex.ps1_unique,
        "before": ex.ztf_frac + ex.des_unique,
        "cumulative": ex.combined,
        "remaining": 1.0 - ex.combined,
        "ps1_total": ps1_total / N_POP as f64,
        "published": {
            "ztf": paper.ztf_frac,
            "des_unique": paper.des_unique,
            "ps1_unique": paper.ps1_unique,
            "cumulative": paper.combined,
            "ps1_total": paper_counts.n_detected as f64 / 100_000.0,
            "a_au": updated.a_median,
            "mass_earth": updated.mass_earth_median,
            "v_mag": updated.v_mag_median,
        },
        "ps1": {
            "depth_v": ps1.depth_limit,
            "dec_limit_deg": ps1.dec_limit_deg,
            "n_epochs": ps1.n_epochs,
            "linking_threshold": ps1.linking_threshold,
            "mag_max": ps1.quality_cuts.mag_max,
        },
        "ztf_depth_v": ztf.depth_limit,
        "des_boxes": des_boxes,
        "population": sky,
        "efficiency": {"v_mag": mags, "ps1": ps1_eff, "ztf": ztf_eff},
        "prior": {
            "a_au": weighted_median(&a, &all),
            "mass_earth": weighted_median(&mass, &all),
            "v_mag": weighted_median(&v_mag, &all),
        },
        "survivors": {
            "a_au": weighted_median(&a, &survive),
            "mass_earth": weighted_median(&mass, &survive),
            "v_mag": weighted_median(&v_mag, &survive),
        },
        "survivor_a": weighted_median(&a, &survive),
        "orbit_mass": orbit.mass_earth,
        "orbit_a": orbit.a,
        "orbit_e": orbit.e,
        "orbit_i_deg": orbit.i.to_degrees(),
    })
}
