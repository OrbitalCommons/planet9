//! Film export for `p9-2021-ztf`: the numbers its scene and ledger entry draw.

use p9_2021_ztf::detection_efficiency::detection_probability_for_orbit;
use p9_2021_ztf::survey_model::ZtfSurvey;
use p9_2024_panstarrs::combined_exclusion::compute_combined_from_population;
use p9_core::constants::GM_SUN;
use p9_core::coords::sky::ecliptic_vec_to_equatorial_deg;
use p9_core::data::reference_population::generate_reference_population;
use p9_core::types::{OrbitalElements, elements_to_cartesian};
use rand::SeedableRng;
use serde_json::{Value, json};

/// Synthetic Planet Nines drawn for the sky panel.
const N_SKY: usize = 700;

pub fn export() -> Value {
    let ex = compute_combined_from_population(3000, 2024);
    let ztf = ZtfSurvey::default();

    // The Brown & Batygin (2021) reference population, each member scored by
    // the ZTF survey model: where it is, how bright, how likely ZTF caught it.
    let mut rng = rand::rngs::StdRng::seed_from_u64(2021);
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
            let pos = elements_to_cartesian(&elements, GM_SUN).pos;
            let (ra, dec) = ecliptic_vec_to_equatorial_deg(&pos);
            json!({
                "ra_deg": ra,
                "dec_deg": dec,
                "dist_au": pos.norm(),
                "v_mag": obj.v_magnitude,
                "p_detect": detection_probability_for_orbit(&ztf, &elements, obj.v_magnitude),
            })
        })
        .collect();

    // Detection efficiency against magnitude inside the footprint.
    let mags: Vec<f64> = (0..=60).map(|k| 17.0 + 0.1 * k as f64).collect();
    let efficiency: Vec<f64> = mags.iter().map(|&m| ztf.detection_efficiency(m)).collect();

    json!({
        "ztf": ex.ztf_frac,
        "cumulative": ex.ztf_frac,
        "dec_limit_deg": ztf.dec_limit_deg,
        "population": population,
        "efficiency": {"v_mag": mags, "fraction": efficiency},
    })
}
