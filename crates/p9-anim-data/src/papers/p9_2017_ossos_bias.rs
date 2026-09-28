//! Film export for `p9-2017-ossos-bias`: the numbers its scene and ledger entry draw.

use p9_2017_ossos_bias::bias::run_default;
use p9_2017_ossos_bias::detection::{DetectionOutcome, SurveyModel, detection_outcome};
use p9_2017_ossos_bias::population::{
    PopulationParams, SyntheticEtno, generate_population, longitudes_of_perihelion,
};
use p9_2019_clustering::ossos_comparison::ossos_sample;
use p9_core::analysis::circular::mean_resultant_length;
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::constants::GM_SUN;
use p9_core::coords::sky::{ecliptic_to_equatorial_deg, ecliptic_vec_to_equatorial_deg};
use p9_core::types::{OrbitalElements, elements_to_cartesian};
use rand::SeedableRng;
use rand::seq::SliceRandom;
use serde_json::{Value, json};

/// Size of the isotropic parent population.
const N_PARENT: usize = 60_000;
const SEED: u64 = 42;
/// Parent members drawn on the sky panel.
const N_SKY_PARENT: usize = 700;
/// Detections drawn on the sky panel.
const N_SKY_DETECTED: usize = 300;
/// Longitude-of-perihelion histogram bins.
const N_BINS: usize = 36;
/// Synthetic surveys drawn for the chance-alignment test.
const N_DRAWS: usize = 100_000;

/// Where an orbit's perihelion sits on the sky (RA, Dec in degrees).
fn perihelion_radec(elements: &OrbitalElements) -> (f64, f64) {
    let at_perihelion = OrbitalElements {
        mean_anomaly: 0.0,
        ..*elements
    };
    ecliptic_vec_to_equatorial_deg(&elements_to_cartesian(&at_perihelion, GM_SUN).pos)
}

fn sky_points(pop: &[&SyntheticEtno], n: usize) -> Vec<Value> {
    pop.iter()
        .take(n)
        .map(|o| {
            let (ra, dec) = perihelion_radec(&o.elements);
            json!({"ra_deg": ra, "dec_deg": dec})
        })
        .collect()
}

fn varpi_histogram(varpis: &[f64]) -> Vec<f64> {
    let mut counts = vec![0.0; N_BINS];
    for v in varpis {
        let k =
            ((v.to_degrees().rem_euclid(360.0) / 360.0 * N_BINS as f64) as usize).min(N_BINS - 1);
        counts[k] += 1.0 / varpis.len() as f64;
    }
    counts
}

pub fn export() -> Value {
    let model = SurveyModel::default();
    let experiment = run_default(N_PARENT, SEED);

    // The same parent the experiment draws, kept so the scene can show the
    // objects and why each one was or was not found.
    let mut rng = rand::rngs::StdRng::seed_from_u64(SEED);
    let parent = generate_population(N_PARENT, &PopulationParams::default(), &mut rng);
    let outcomes: Vec<DetectionOutcome> = parent
        .iter()
        .map(|o| detection_outcome(o, &model))
        .collect();
    let share = |which: DetectionOutcome| {
        outcomes.iter().filter(|&&o| o == which).count() as f64 / N_PARENT as f64
    };
    let detected: Vec<&SyntheticEtno> = parent
        .iter()
        .zip(&outcomes)
        .filter(|(_, o)| **o == DetectionOutcome::Detected)
        .map(|(p, _)| p)
        .collect();
    let bright: Vec<&SyntheticEtno> = parent
        .iter()
        .zip(&outcomes)
        .filter(|(_, o)| **o != DetectionOutcome::TooFaint)
        .map(|(p, _)| p)
        .collect();

    let parent_varpi = longitudes_of_perihelion(&parent);
    let detected_varpi: Vec<f64> = detected
        .iter()
        .map(|o| o.elements.longitude_of_perihelion())
        .collect();

    // Outline of each searched block: ecliptic longitude window by the
    // ecliptic latitude band, in equatorial coordinates.
    let band = model.ecliptic_lat_band.to_degrees();
    let blocks: Vec<Value> = model
        .windows
        .iter()
        .map(|w| {
            let (c, hw) = (w.center.to_degrees(), w.half_width.to_degrees());
            let mut outline = Vec::new();
            let steps = 12;
            for k in 0..=steps {
                let f = k as f64 / steps as f64;
                outline.push((c - hw + 2.0 * hw * f, -band));
            }
            for k in 0..=steps {
                let f = k as f64 / steps as f64;
                outline.push((c + hw, -band + 2.0 * band * f));
            }
            for k in 0..=steps {
                let f = k as f64 / steps as f64;
                outline.push((c + hw - 2.0 * hw * f, band));
            }
            for k in 0..=steps {
                let f = k as f64 / steps as f64;
                outline.push((c - hw, band - 2.0 * band * f));
            }
            let radec: Vec<[f64; 2]> = outline
                .iter()
                .map(|&(lon, lat)| {
                    let (ra, dec) = ecliptic_to_equatorial_deg(lon.rem_euclid(360.0), lat);
                    [ra.rem_euclid(360.0), dec]
                })
                .collect();
            json!({"lon_deg": c, "half_width_deg": hw, "outline": radec})
        })
        .collect();

    // The real survey: the OSSOS objects beyond 230 AU. How often does a
    // survey of the same size, run on the uniform parent, come back at least
    // as aligned?
    let ossos = ossos_sample();
    let ossos_varpi: Vec<f64> = ossos
        .iter()
        .map(|k| k.elements.longitude_of_perihelion())
        .collect();
    let r_bar_ossos = mean_resultant_length(&ossos_varpi);
    let mut hits = 0usize;
    for _ in 0..N_DRAWS {
        let drawn: Vec<f64> = detected_varpi
            .choose_multiple(&mut rng, ossos.len())
            .copied()
            .collect();
        if mean_resultant_length(&drawn) >= r_bar_ossos {
            hits += 1;
        }
    }
    let p_chance = hits as f64 / N_DRAWS as f64;

    let edges: Vec<f64> = (0..=N_BINS)
        .map(|k| 360.0 * k as f64 / N_BINS as f64)
        .collect();

    json!({
        "n_parent": N_PARENT,
        "n_detected": experiment.detected.n,
        "r_bar_parent": experiment.intrinsic.r_bar,
        "r_bar_detected": experiment.detected.r_bar,
        "rayleigh_p_parent": experiment.intrinsic.rayleigh_p,
        "rayleigh_p_detected": experiment.detected.rayleigh_p,
        "mean_varpi_detected_deg": experiment.detected.mean_varpi.map(|m| m.to_degrees().rem_euclid(360.0)),
        "limiting_mag": model.limiting_mag,
        "dec_min_deg": model.dec_min.to_degrees(),
        "lat_band_deg": band,
        "blocks": blocks,
        "lost": {
            "too_faint": share(DetectionOutcome::TooFaint),
            "below_declination": share(DetectionOutcome::BelowDeclination),
            "off_ecliptic": share(DetectionOutcome::OffEclipticBand),
            "outside_blocks": share(DetectionOutcome::OutsideLongitudeWindows),
            "detected": share(DetectionOutcome::Detected),
        },
        "sky_bright": sky_points(&bright, N_SKY_PARENT),
        "sky_detected": sky_points(&detected, N_SKY_DETECTED),
        "varpi": {
            "edges_deg": edges,
            "parent": varpi_histogram(&parent_varpi),
            "detected": varpi_histogram(&detected_varpi),
        },
        "ossos": ossos.iter().map(|k| json!({
            "name": k.name,
            "a": k.elements.a,
            "varpi_deg": k.elements.longitude_of_perihelion().to_degrees().rem_euclid(360.0),
        })).collect::<Vec<_>>(),
        "n_ossos": ossos.len(),
        "r_bar_ossos": r_bar_ossos,
        "p_chance": p_chance,
        "sigma": p_value_to_sigma(p_chance).max(0.0),
    })
}
