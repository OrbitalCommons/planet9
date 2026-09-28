//! Film export for `p9-2022-des`: the numbers its scene and ledger entry draw.

use p9_2021_ztf::detection_efficiency::detection_probability_for_orbit;
use p9_2021_ztf::survey_model::ZtfSurvey;
use p9_2022_des::color_models::{
    ColorModel, fiducial, methane_40k, neptune_like, super_ganymede, super_kbo,
};
use p9_2022_des::recovery_analysis::{RecoveryResult, compute_recovery_for_model};
use p9_2022_des::survey_model::{DesBand, DesSurvey};
use p9_2024_panstarrs::combined_exclusion::compute_combined_from_population;
use p9_core::analysis::surveys::DES_FOOTPRINT_BANDS;
use p9_core::constants::GM_SUN;
use p9_core::coords::sky::ecliptic_vec_to_equatorial_deg;
use p9_core::data::reference_population::generate_reference_population;
use p9_core::types::OrbitalElements;
use rand::SeedableRng;
use serde_json::{Value, json};

/// Synthetic Planet Nines drawn for the sky panel.
const N_SKY: usize = 900;

/// Synthetic Planet Nines injected for the recovery statistics.
const N_RECOVERY: usize = 4000;

/// Published size of the injected catalogue and the fraction of it that
/// crosses the footprint (11,709 of 100,000).
const PUBLISHED_CROSSING_FRACTION: f64 = 0.117_09;

/// Published additional exclusion after removing what ZTF already covered
/// ("rules out an additional 5% of the parameter space").
const PUBLISHED_DES_UNIQUE: f64 = 0.05;

fn recovery(model: &ColorModel) -> Value {
    let r = compute_recovery_for_model(model, N_RECOVERY, 2022);
    json!({
        "name": model.name,
        "albedo": if model.per_object_albedo { Value::Null } else { json!(model.albedo) },
        "g_minus_r": model.g_minus_r,
        "recovery": r.recovery_frac,
        "published_recovery": model.paper_recovery,
        "n_crossing": r.n_crossing,
        "n_recovered": r.n_recovered,
    })
}

pub fn export() -> Value {
    let ex = compute_combined_from_population(3000, 2024);
    let des = DesSurvey::default();
    let ztf = ZtfSurvey::default();
    let colours = fiducial();

    let mut rng = rand::rngs::StdRng::seed_from_u64(2022);
    let population: Vec<Value> = generate_reference_population(N_SKY, &mut rng)
        .iter()
        .map(|obj| {
            let elements = OrbitalElements {
                a: obj.a,
                e: obj.e,
                i: obj.i,
                omega: obj.omega,
                omega_big: obj.omega_big,
                mean_anomaly: obj.mean_anomaly,
            };
            let pos = elements.to_state_vector(GM_SUN).pos;
            let (ra, dec) = ecliptic_vec_to_equatorial_deg(&pos);
            let mags = colours.band_magnitudes_with_albedo(obj.mass, pos.norm(), obj.albedo);
            json!({
                "ra_deg": ra,
                "dec_deg": dec,
                "dist_au": pos.norm(),
                "v_mag": obj.v_magnitude,
                "r_mag": mags.magnitude(DesBand::R),
                "in_footprint": des.crosses_footprint(&elements),
                "p_des": des.detection_probability_for_orbit(&elements, &mags),
                "p_ztf": detection_probability_for_orbit(&ztf, &elements, obj.v_magnitude),
            })
        })
        .collect();

    let footprint: Vec<Value> = DES_FOOTPRINT_BANDS
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

    let fid = compute_recovery_for_model(&colours, N_RECOVERY, 2022);
    let paper = RecoveryResult::paper_result();
    let models: Vec<Value> = [
        fiducial(),
        methane_40k(),
        super_ganymede(),
        neptune_like(),
        super_kbo(),
    ]
    .iter()
    .map(recovery)
    .collect();

    let mags: Vec<f64> = (0..=50).map(|k| 21.0 + 0.1 * k as f64).collect();
    let completeness: Vec<f64> = mags
        .iter()
        .map(|&m| des.completeness(m, DesBand::R))
        .collect();

    json!({
        "des_unique": ex.des_unique,
        "published_des_unique": PUBLISHED_DES_UNIQUE,
        "cumulative": ex.ztf_frac + ex.des_unique,
        "ztf": ex.ztf_frac,
        "population": population,
        "footprint": footprint,
        "footprint_area_deg2": des.solid_angle_deg2(),
        "footprint_sky_fraction": des.solid_angle_deg2() / 41_252.96,
        "depth_r": des.depth_r,
        "min_nights": des.min_nights,
        "n_injected": N_RECOVERY,
        "n_crossing": fid.n_crossing,
        "n_recovered": fid.n_recovered,
        "crossing_fraction": fid.n_crossing as f64 / N_RECOVERY as f64,
        "recovery": fid.recovery_frac,
        "published_crossing_fraction": PUBLISHED_CROSSING_FRACTION,
        "published_recovery": paper.recovery_frac,
        "published_n_crossing": paper.n_crossing,
        "published_n_recovered": paper.n_recovered,
        "models": models,
        "completeness": {"r_mag": mags, "fraction": completeness},
    })
}
