//! Film export for `p9-2018-wise-search`: the numbers its scene and ledger entry draw.

use p9_2018_wise_search::detectability::{detection_probability, max_detectable_distance};
use p9_2018_wise_search::sky::galactic_latitude_deg;
use p9_2018_wise_search::survey_model::WiseSurvey;
use p9_2018_wise_search::thermal_model::P9Thermal;
use p9_core::analysis::thermal::max_detectable_distance as crossing_distance;
use p9_core::constants::GM_SUN;
use p9_core::coords::sky::{ecliptic_vec_to_equatorial_deg, equatorial_to_galactic};
use p9_core::data::reference_population::generate_reference_population;
use p9_core::types::OrbitalElements;
use rand::SeedableRng;
use serde_json::{Value, json};

/// Synthetic Planet Nines placed on the sky and in the mass-distance plane.
const N_POPULATION: usize = 600;

/// Published 90%-completeness depth at high galactic latitude (W1, Vega).
const PUBLISHED_DEPTH_W1: f64 = 16.7;

/// Published search footprint: 31,480 square degrees, 76% of the sky.
const PUBLISHED_AREA_DEG2: f64 = 31_480.0;
const PUBLISHED_SKY_FRACTION: f64 = 0.76;

/// The most W1-luminous Fortney et al. (2016) model atmosphere as quoted by
/// Meisner et al.: W1 = 16.1 at 622 AU, fading as the inverse square.
const LUMINOUS_W1: f64 = 16.1;
const LUMINOUS_AT_AU: f64 = 622.0;

/// Published reach for that model (AU): "800-900 AU depending on latitude".
const PUBLISHED_REACH_AU: (f64, f64) = (800.0, 900.0);

/// Nearest distance the linking was designed to recover (AU).
const PUBLISHED_MIN_DISTANCE_AU: f64 = 250.0;

/// Cell size of the footprint grid (degrees).
const CELL_DEG: f64 = 5.0;

fn luminous_w1(distance_au: f64) -> f64 {
    LUMINOUS_W1 + 5.0 * (distance_au / LUMINOUS_AT_AU).log10()
}

pub fn export() -> Value {
    let survey = WiseSurvey::default();
    let published = WiseSurvey {
        w1_depth: PUBLISHED_DEPTH_W1,
        ..WiseSurvey::default()
    };

    // The masked galactic band on an RA/Dec grid (1 = masked).
    let n_ra = (360.0 / CELL_DEG) as usize;
    let n_dec = (180.0 / CELL_DEG) as usize;
    let ra_centres: Vec<f64> = (0..n_ra).map(|k| (k as f64 + 0.5) * CELL_DEG).collect();
    let dec_centres: Vec<f64> = (0..n_dec)
        .map(|k| -90.0 + (k as f64 + 0.5) * CELL_DEG)
        .collect();
    let mut masked = Vec::with_capacity(n_ra * n_dec);
    for &dec in &dec_centres {
        for &ra in &ra_centres {
            let (_, b) = equatorial_to_galactic(ra.to_radians(), dec.to_radians());
            masked.push(if survey.in_footprint(b.to_degrees()) {
                0.0
            } else {
                1.0
            });
        }
    }

    let reach_luminous = crossing_distance(50.0, 5000.0, PUBLISHED_DEPTH_W1, luminous_w1);

    let mut rng = rand::rngs::StdRng::seed_from_u64(2017);
    let mut sum_reflected = 0.0;
    let mut n_luminous = 0usize;
    let population: Vec<Value> = generate_reference_population(N_POPULATION, &mut rng)
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
            let dist = pos.norm();
            let b = galactic_latitude_deg(&elements);
            let body = P9Thermal::new(obj.mass, dist);
            let p_reflected = detection_probability(&published, &body, b);
            let in_luminous_reach =
                published.in_footprint(b) && reach_luminous.is_some_and(|r| dist <= r);
            sum_reflected += p_reflected;
            n_luminous += usize::from(in_luminous_reach);
            json!({
                "ra_deg": ra,
                "dec_deg": dec,
                "dist_au": dist,
                "mass_earth": obj.mass,
                "galactic_lat_deg": b,
                "w1_reflected": body.w1_magnitude(),
                "p_reflected": p_reflected,
                "in_luminous_reach": in_luminous_reach,
            })
        })
        .collect();

    let masses: Vec<f64> = (0..=34).map(|k| 3.0 + 0.5 * k as f64).collect();
    let reach_reflected: Vec<Option<f64>> = masses
        .iter()
        .map(|&m| max_detectable_distance(&published, m))
        .collect();

    // W1 brightness of a 10 Earth-mass planet against distance: sunlight only
    // (the crate's thermal model) and the self-luminous model atmosphere.
    let distances: Vec<f64> = (0..=50).map(|k| 100.0 + 20.0 * k as f64).collect();
    let w1_reflected_10me: Vec<f64> = distances
        .iter()
        .map(|&d| P9Thermal::new(10.0, d).w1_magnitude())
        .collect();
    let w1_luminous: Vec<f64> = distances.iter().map(|&d| luminous_w1(d)).collect();

    json!({
        "curve": {
            "distance_au": distances,
            "w1_reflected_10me": w1_reflected_10me,
            "w1_luminous": w1_luminous,
        },
        "depth_w1": survey.w1_depth,
        "published_depth_w1": PUBLISHED_DEPTH_W1,
        "galactic_mask_deg": survey.galactic_mask_deg,
        "sky_fraction": survey.sky_coverage_fraction(),
        "published_sky_fraction": PUBLISHED_SKY_FRACTION,
        "published_area_deg2": PUBLISHED_AREA_DEG2,
        "mask": {"ra_centres": ra_centres, "dec_centres": dec_centres, "masked": masked},
        "population": population,
        "fraction_reflected": sum_reflected / N_POPULATION as f64,
        "fraction_luminous": n_luminous as f64 / N_POPULATION as f64,
        "reach_luminous_au": reach_luminous,
        "published_reach_au": [PUBLISHED_REACH_AU.0, PUBLISHED_REACH_AU.1],
        "published_min_distance_au": PUBLISHED_MIN_DISTANCE_AU,
        "reach_reflected_10me_au": max_detectable_distance(&published, 10.0),
        "reach_vs_mass": {"mass_earth": masses, "reflected_au": reach_reflected},
        "luminous_anchor": {"w1": LUMINOUS_W1, "distance_au": LUMINOUS_AT_AU},
    })
}
