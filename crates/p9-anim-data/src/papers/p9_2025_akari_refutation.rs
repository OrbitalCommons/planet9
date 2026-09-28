//! Film export for `p9-2025-akari-refutation`: the numbers its scene and ledger entry draw.

use p9_2025_akari_refutation::expected_motion::{
    parallax_range_arcmin, parallax_six_month_arcmin, proper_motion_range_arcmin,
    proper_motion_six_month_arcmin, proper_motion_six_month_circular_arcmin,
};
use p9_2025_akari_refutation::published;
use p9_2025_iras_akari::orbital_constraints::CandidatePair;
use p9_2025_iras_akari::survey_model::angular_separation_arcmin;
use serde_json::{Value, json};

/// Search region of Chen et al. (2025), Sec. 3.3 (degrees).
const REGION_RA_DEG: (f64, f64) = (30.0, 50.0);
const REGION_DEC_DEG: (f64, f64) = (-20.0, 20.0);
/// Distance range over which Chen et al. quote the six-month motion (AU).
const SEARCH_DISTANCE_AU: (f64, f64) = (300.0, 800.0);
/// Eccentricity of the fastest orbit in the quoted proper-motion range.
const FAST_ECCENTRICITY: f64 = 0.6;

/// One of the two candidates of Chen et al. (2025), Table 2.
struct Candidate {
    name: &'static str,
    ra_deg: f64,
    dec_deg: f64,
    flux_jy: (f64, f64),
    epoch: &'static str,
}

const CANDIDATES: [Candidate; 2] = [
    Candidate {
        name: "FISSSDL J0250422-150114",
        ra_deg: 42.676,
        dec_deg: -15.021,
        flux_jy: (0.61, 1.62),
        epoch: "2007-07-28",
    },
    Candidate {
        name: "FISSSDL J0301112-164240",
        ra_deg: 45.297,
        dec_deg: -16.711,
        flux_jy: (1.27, 0.51),
        epoch: "2006-07-30",
    },
];

/// Source counts after each selection step (Chen et al. 2025, Secs. 3.3-3.7).
const FUNNEL: [(&str, u64); 6] = [
    ("AKARI single-scan detections", 5_274_338),
    ("inside the search region", 50_033),
    ("not in 9 other catalogues", 29_901),
    ("bright, low background", 1_726),
    ("gone after six months, not cosmic rays", 13),
    ("candidates", 2),
];

pub fn export() -> Value {
    // Six-month motion of a bound body against distance.
    let distances: Vec<f64> = (0..=55).map(|k| 250.0 + 10.0 * k as f64).collect();
    let parallax: Vec<f64> = distances
        .iter()
        .map(|&d| parallax_six_month_arcmin(d))
        .collect();
    let pm_circular: Vec<f64> = distances
        .iter()
        .map(|&d| proper_motion_six_month_circular_arcmin(d))
        .collect();
    let pm_perihelion: Vec<f64> = distances
        .iter()
        .map(|&d| {
            proper_motion_six_month_arcmin(d, d / (1.0 - FAST_ECCENTRICITY), FAST_ECCENTRICITY)
        })
        .collect();
    let (par_lo, par_hi) = parallax_range_arcmin(SEARCH_DISTANCE_AU.0, SEARCH_DISTANCE_AU.1);
    let (pm_lo, pm_hi) = proper_motion_range_arcmin(
        SEARCH_DISTANCE_AU.0,
        SEARCH_DISTANCE_AU.1,
        FAST_ECCENTRICITY,
    );

    // The two candidates against each other.
    let (a, b) = (&CANDIDATES[0], &CANDIDATES[1]);
    let pair_separation = angular_separation_arcmin(a.ra_deg, a.dec_deg, b.ra_deg, b.dec_deg);

    // Where the earlier IRAS-AKARI candidate (Phan et al. 2025) sits.
    let earlier = CandidatePair::paper_candidate();

    json!({
        "region": {
            "ra_lo": REGION_RA_DEG.0,
            "ra_hi": REGION_RA_DEG.1,
            "dec_lo": REGION_DEC_DEG.0,
            "dec_hi": REGION_DEC_DEG.1,
        },
        "candidates": CANDIDATES.iter().map(|c| json!({
            "name": c.name,
            "ra_deg": c.ra_deg,
            "dec_deg": c.dec_deg,
            "flux_jy": [c.flux_jy.0, c.flux_jy.1],
            "epoch": c.epoch,
        })).collect::<Vec<_>>(),
        "pair_separation_arcmin": pair_separation,
        "funnel": FUNNEL.iter().map(|(label, n)| json!({"label": label, "n": n}))
            .collect::<Vec<_>>(),
        "motion": {
            "distance_au": distances,
            "parallax_arcmin": parallax,
            "pm_circular_arcmin": pm_circular,
            "pm_perihelion_arcmin": pm_perihelion,
            "eccentricity": FAST_ECCENTRICITY,
        },
        "search_distance_au": [SEARCH_DISTANCE_AU.0, SEARCH_DISTANCE_AU.1],
        "parallax_range_arcmin": [par_lo, par_hi],
        "proper_motion_range_arcmin": [pm_lo, pm_hi],
        "parallax_max": par_hi,
        "published": {
            "parallax_arcmin": [
                published::PARALLAX_6MO_MIN_ARCMIN,
                published::PARALLAX_6MO_MAX_ARCMIN,
            ],
            "proper_motion_arcmin": [
                published::PROPER_MOTION_6MO_MIN_ARCMIN,
                published::PROPER_MOTION_6MO_MAX_ARCMIN,
            ],
        },
        "earlier_pair": {
            "iras_ra_deg": earlier.iras_ra_deg,
            "iras_dec_deg": earlier.iras_dec_deg,
        },
    })
}
