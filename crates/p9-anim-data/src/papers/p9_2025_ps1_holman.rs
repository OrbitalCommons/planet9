//! Film export for `p9-2025-ps1-holman`: the numbers its scene and ledger entry draw.

use nalgebra::Vector3;
use p9_2021_ztf::detection_efficiency::detection_probability_for_orbit as ztf_probability;
use p9_2021_ztf::survey_model::ZtfSurvey;
use p9_2022_des::color_models::fiducial;
use p9_2022_des::survey_model::DesSurvey;
use p9_2024_panstarrs::detection_pipeline::{
    detection_probability_for_orbit, equatorial_position_deg,
};
use p9_2024_panstarrs::survey_model::Ps1Survey;
use p9_2025_ps1_holman::exclusion::{compute_exclusion, reference_targets};
use p9_2025_ps1_holman::photometry::max_detectable_distance_typed;
use p9_2025_ps1_holman::survey_model::Ps1StackSurvey;
use p9_core::coords::sky::equatorial_to_galactic_matrix;
use p9_core::data::reference_population::{generate_reference_population, heliocentric_distance};
use p9_core::types::{OrbitalElements, P9Params};
use p9_core::units::au;
use rand::SeedableRng;
use serde_json::{Value, json};

/// Size and seed of the scored population: the draw the ZTF, DES and
/// Pan-STARRS1 exclusion fractions are averaged over.
const N_POP: usize = 3000;
const SEED: u64 = 2024;
/// Members drawn on the sky panel.
const N_SKY: usize = 800;
/// Half-width of the galactic-plane band used to summarise where the
/// survivors sit (degrees).
const PLANE_B_DEG: f64 = 10.0;

/// Holman et al. (2025), Sec. VIII: reference-population members the survey
/// simulator recovers, and the number Brown, Holman & Batygin (2024) ruled
/// out, both out of 100,000.
const PUBLISHED_RECOVERED: f64 = 75_769.0;
const PUBLISHED_EARLIER_PS1: f64 = 68_745.0;
const PUBLISHED_POPULATION: f64 = 100_000.0;
/// Holman et al. (2025), abstract and Sec. II.
const PUBLISHED_OBJECTS: u32 = 692;
const PUBLISHED_TNOS: u32 = 642;
const PUBLISHED_DWARF_PLANETS: u32 = 23;
const PUBLISHED_NEW: u32 = 109;
const PUBLISHED_EXPOSURES: u32 = 708_554;
const PUBLISHED_DEPTH_W: f64 = 22.5;
const PUBLISHED_DISTANCE_AU: (f64, f64) = (80.0, 1600.0);

fn galactic_latitude_deg(ra_deg: f64, dec_deg: f64) -> f64 {
    let (ra, dec) = (ra_deg.to_radians(), dec_deg.to_radians());
    let eq = Vector3::new(dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin());
    let gal = equatorial_to_galactic_matrix() * eq;
    (gal.z / gal.norm()).asin().to_degrees()
}

pub fn export() -> Value {
    let survey = Ps1StackSurvey::default();
    let targets = reference_targets(N_POP, SEED);
    let alone = compute_exclusion(&survey, &targets);

    let ztf = ZtfSurvey::default();
    let des = DesSurvey::default();
    let des_color = fiducial();
    let ps1 = Ps1Survey::default();

    let mut rng = rand::rngs::StdRng::seed_from_u64(SEED);
    let population = generate_reference_population(N_POP, &mut rng);

    let mut sky = Vec::with_capacity(N_SKY);
    let (mut before, mut after, mut new) = (0.0, 0.0, 0.0);
    let (mut plane_all, mut plane_left, mut left) = (0.0, 0.0, 0.0);
    for (k, (obj, target)) in population.iter().zip(&targets).enumerate() {
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
        let p_before = 1.0 - (1.0 - p_ztf) * (1.0 - p_des) * (1.0 - p_ps1);
        let p_here = survey.detection_probability(target.r_magnitude, target.dec_deg);
        let p_after = 1.0 - (1.0 - p_before) * (1.0 - p_here);

        let (ra, dec) = equatorial_position_deg(&elements);
        let in_plane = galactic_latitude_deg(ra, dec).abs() < PLANE_B_DEG;
        before += p_before;
        after += p_after;
        new += p_here * (1.0 - p_before);
        left += 1.0 - p_after;
        if in_plane {
            plane_all += 1.0;
            plane_left += 1.0 - p_after;
        }

        if k < N_SKY {
            sky.push(json!({
                "ra_deg": ra,
                "dec_deg": dec,
                "r_mag": target.r_magnitude,
                "p_before": p_before,
                "p_here": p_here,
                "p_after": p_after,
            }));
        }
    }
    let n = N_POP as f64;

    let mags: Vec<f64> = (0..=60).map(|k| 18.0 + 0.1 * k as f64).collect();
    let efficiency: Vec<f64> = mags
        .iter()
        .map(|&m| survey.recovery_efficiency(m))
        .collect();
    let masses: Vec<f64> = (0..=32).map(|k| 2.0 + 0.5 * k as f64).collect();
    let reach: Vec<f64> = masses
        .iter()
        .map(|&m| (max_detectable_distance_typed(&survey, m) / au(1.0)).value)
        .collect();

    json!({
        "alone": alone.fraction_excluded,
        "before": before / n,
        "cumulative": after / n,
        "new": new / n,
        "remaining": left / n,
        "plane_b_deg": PLANE_B_DEG,
        "plane_share_of_population": plane_all / n,
        "plane_share_of_survivors": plane_left / left,
        "survey": {
            "single_epoch_depth": survey.single_epoch_depth,
            "depth_gain": survey.stack_depth_gain,
            "effective_depth": survey.effective_depth(),
            "dec_limit_deg": survey.dec_limit_deg,
        },
        "population": sky,
        "efficiency": {"r_mag": mags, "fraction": efficiency},
        "reach": {"mass_earth": masses, "distance_au": reach},
        "published": {
            "alone": PUBLISHED_RECOVERED / PUBLISHED_POPULATION,
            "earlier_ps1": PUBLISHED_EARLIER_PS1 / PUBLISHED_POPULATION,
            "objects": PUBLISHED_OBJECTS,
            "tnos": PUBLISHED_TNOS,
            "dwarf_planets": PUBLISHED_DWARF_PLANETS,
            "new": PUBLISHED_NEW,
            "exposures": PUBLISHED_EXPOSURES,
            "depth_w": PUBLISHED_DEPTH_W,
            "distance_au": [PUBLISHED_DISTANCE_AU.0, PUBLISHED_DISTANCE_AU.1],
        },
    })
}
