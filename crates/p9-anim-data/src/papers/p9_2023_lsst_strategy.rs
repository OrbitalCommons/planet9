//! Film export for `p9-2023-lsst-strategy`: the numbers its scene and ledger entry draw.

use p9_2023_lsst_strategy::strategy::published::{
    SINGLE_VISIT_DEPTH_R, STACK_DEPTH_R, VISITS_FOR_LINKING, VISITS_PER_FIELD,
};
use p9_2023_lsst_strategy::{
    LsstStrategy, P9Sample, dec_and_galactic_lat_deg, discoverable_fraction, discovery_probability,
};
use p9_core::analysis::photometry::SOLAR_V_MINUS_R;
use p9_core::constants::GM_SUN;
use p9_core::coords::sky::ecliptic_vec_to_equatorial_deg;
use p9_core::data::reference_population::generate_reference_population;
use p9_core::types::OrbitalElements;
use rand::SeedableRng;
use serde_json::{Value, json};

/// Synthetic Planet Nines scored for the strategy scans.
const N_POPULATION: usize = 3000;

/// Of those, the number drawn on the sky panel.
const N_SKY: usize = 800;

/// Published spread of the solar-system metrics across most of the simulated
/// cadences (Schwamb et al. 2023 abstract): within ±5%.
const PUBLISHED_TYPICAL_SPREAD: f64 = 0.05;

pub fn export() -> Value {
    let mut rng = rand::rngs::StdRng::seed_from_u64(2023);
    let population: Vec<P9Sample> = generate_reference_population(N_POPULATION, &mut rng)
        .into_iter()
        .map(|p| P9Sample {
            elements: OrbitalElements {
                a: p.a,
                e: p.e,
                i: p.i,
                omega: p.omega,
                omega_big: p.omega_big,
                mean_anomaly: p.mean_anomaly,
            },
            v_magnitude: p.v_magnitude,
        })
        .collect();

    let baseline = LsstStrategy::baseline();
    let base = discoverable_fraction(&baseline, &population);
    let fraction = |s: &LsstStrategy| discoverable_fraction(s, &population).fraction;

    let sky: Vec<Value> = population
        .iter()
        .take(N_SKY)
        .map(|s| {
            let pos = s.elements.to_state_vector(GM_SUN).pos;
            let (ra, dec) = ecliptic_vec_to_equatorial_deg(&pos);
            let (_, b) = dec_and_galactic_lat_deg(&s.elements);
            json!({
                "ra_deg": ra,
                "dec_deg": dec,
                "galactic_lat_deg": b,
                "r_mag": s.v_magnitude - SOLAR_V_MINUS_R,
                "in_footprint": baseline.footprint().accepts(dec, b),
                "p_discover": discovery_probability(&baseline, s),
            })
        })
        .collect();

    let dec_limits: Vec<f64> = (0..=24).map(|k| -20.0 + 2.5 * k as f64).collect();
    let by_dec: Vec<f64> = dec_limits
        .iter()
        .map(|&d| fraction(&baseline.clone().with_dec_max(d)))
        .collect();

    let depths: Vec<f64> = (0..=30).map(|k| 21.0 + 0.125 * k as f64).collect();
    let by_depth: Vec<f64> = depths
        .iter()
        .map(|&m| fraction(&baseline.clone().with_single_visit_depth(m)))
        .collect();

    let visits: Vec<u32> = (1..=VISITS_PER_FIELD).collect();
    let by_visits: Vec<f64> = visits
        .iter()
        .map(|&n| fraction(&baseline.clone().with_visits_for_linking(n)))
        .collect();

    json!({
        "fraction": base.fraction,
        "fraction_in_footprint": base.fraction_in_footprint,
        "fraction_stack": fraction(&baseline.clone().with_stack_depth()),
        "linked_given_footprint": base.fraction / base.fraction_in_footprint,
        "footprint": {
            "dec_lo": baseline.dec_min_deg,
            "dec_hi": baseline.dec_max_deg,
            "galactic_lat_min_deg": baseline.galactic_lat_min_deg,
            "coverage": baseline.coverage_fraction,
        },
        "single_visit_depth": SINGLE_VISIT_DEPTH_R,
        "stack_depth": STACK_DEPTH_R,
        "visits_for_linking": VISITS_FOR_LINKING,
        "visits_per_field": VISITS_PER_FIELD,
        "population": sky,
        "by_dec_limit": {"dec_max_deg": dec_limits, "fraction": by_dec},
        "by_depth": {"depth_r": depths, "fraction": by_depth},
        "by_visits": {"visits": visits, "fraction": by_visits},
        "published_typical_spread": PUBLISHED_TYPICAL_SPREAD,
    })
}
