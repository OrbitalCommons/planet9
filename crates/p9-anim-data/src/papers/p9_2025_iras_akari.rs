//! Film export for `p9-2025-iras-akari`: the numbers its scene and ledger entry draw.

use p9_2025_iras_akari::candidate_search::GalacticDensityModel;
use p9_2025_iras_akari::candidate_search::{
    PAPER_POST_CUT_AKARI, PAPER_POST_CUT_IRAS, PaperCuts, expected_chance_pairs,
    false_alarm_probability,
};
use p9_2025_iras_akari::orbital_constraints::{
    CandidatePair, derive_constraints, implied_distance_circular,
};
use p9_2025_iras_akari::survey_model::{AkariFisSurvey, IrasSurvey, epoch_baseline};
use p9_2025_iras_akari::thermal_model::{P9ThermalParams, flux_ratio_60_90};
use p9_core::coords::candidate_pair::{
    E_MAX, annual_proper_motion_circular, distance_range_for_separation, implied_distance,
    separation_window,
};
use serde_json::{Value, json};

/// Candidate pairs surviving the paper's selection cuts, and the number that
/// survive image inspection (Phan et al. 2025, abstract).
const PUBLISHED_PAIRS: u32 = 13;
const PUBLISHED_GOOD: u32 = 1;
/// The paper's quoted search range.
const PUBLISHED_DISTANCE_AU: (f64, f64) = (500.0, 700.0);
const PUBLISHED_MASS_EARTH: (f64, f64) = (7.0, 17.0);
const PUBLISHED_BASELINE_YEARS: f64 = 23.0;

pub fn export() -> Value {
    let iras = IrasSurvey::default();
    let akari = AkariFisSurvey::default();
    let cuts = PaperCuts::default();
    let baseline = epoch_baseline(&iras, &akari);
    let (sep_lo, sep_hi) = cuts.sep_window_arcmin;

    // The candidate, from the two catalogue positions.
    let pair = CandidatePair::paper_candidate();
    let motion = derive_constraints(&pair);

    // How far a bound body can move between the two epochs, against distance.
    let distances: Vec<f64> = (0..=24).map(|k| 350.0 + 25.0 * k as f64).collect();
    let mut window_min = Vec::new();
    let mut window_max = Vec::new();
    let mut circular = Vec::new();
    for &d in &distances {
        let (lo, hi) = separation_window(d, E_MAX, baseline);
        window_min.push(lo);
        window_max.push(hi);
        circular.push(annual_proper_motion_circular(d) * baseline);
    }
    let (d_near, d_far) = distance_range_for_separation(sep_lo, sep_hi, baseline);

    // Far-infrared brightness against distance for the paper's mass range.
    let temps = [30.0, 40.0, 50.0];
    let flux_d: Vec<f64> = (0..=40).map(|k| 300.0 + 15.0 * k as f64).collect();
    let mut flux = Vec::new();
    for &mass in &[PUBLISHED_MASS_EARTH.0, 10.0, PUBLISHED_MASS_EARTH.1] {
        for &t_eff in &temps {
            let at = |d: f64, wavelength_um: f64| {
                P9ThermalParams {
                    mass_earth: mass,
                    distance_au: d,
                    t_eff,
                }
                .flux_density_jy(wavelength_um * 1.0e-6)
            };
            flux.push(json!({
                "mass_earth": mass,
                "t_eff": t_eff,
                "f60_jy": flux_d.iter().map(|&d| at(d, iras.wavelength_um)).collect::<Vec<_>>(),
                "f90_jy": flux_d.iter().map(|&d| at(d, akari.wavelength_um)).collect::<Vec<_>>(),
            }));
        }
    }

    let chance = expected_chance_pairs(
        PAPER_POST_CUT_IRAS,
        PAPER_POST_CUT_AKARI,
        sep_lo,
        sep_hi,
        &GalacticDensityModel::default(),
    );

    json!({
        "baseline_years": baseline,
        "iras": {
            "wavelength_um": iras.wavelength_um,
            "limit_jy": iras.sensitivity_jy,
            "position_error_arcsec": iras.position_error_arcsec,
            "year": 1983,
        },
        "akari": {
            "wavelength_um": akari.wavelength_um,
            "limit_jy": akari.sensitivity_jy,
            "position_error_arcsec": akari.position_error_arcsec,
            "year": 2006,
        },
        "window_arcmin": [sep_lo, sep_hi],
        "window_distance_au": [d_near, d_far],
        "separation": {
            "distance_au": distances,
            "min_arcmin": window_min,
            "max_arcmin": window_max,
            "circular_arcmin": circular,
        },
        "candidate": {
            "iras_ra_deg": pair.iras_ra_deg,
            "iras_dec_deg": pair.iras_dec_deg,
            "akari_ra_deg": pair.akari_ra_deg,
            "akari_dec_deg": pair.akari_dec_deg,
            "iras_flux_60": pair.iras_flux_60,
            "iras_flux_100": pair.iras_flux_100,
            "akari_flux_65": pair.akari_flux_65,
            "akari_flux_90": pair.akari_flux_90,
            "separation_arcmin": motion.separation_arcmin,
            "rate_arcmin_yr": motion.proper_motion_arcmin_yr,
            "position_angle_deg": motion.position_angle_deg,
            "delta_ra_arcmin": motion.delta_ra_arcmin,
            "delta_dec_arcmin": motion.delta_dec_arcmin,
            "distance_au": implied_distance(motion.separation_arcmin, baseline),
            "distance_circular_au": implied_distance_circular(motion.proper_motion_arcmin_yr),
            "flux_ratio_60_90": pair.iras_flux_60 / pair.akari_flux_90,
        },
        "candidate_separation": motion.separation_arcmin,
        "candidate_distance": implied_distance(motion.separation_arcmin, baseline),
        "blackbody_ratio_60_90": temps.iter().map(|&t| flux_ratio_60_90(t)).collect::<Vec<_>>(),
        "temperatures_k": temps,
        "flux": {"distance_au": flux_d, "curves": flux},
        "post_cut_iras": PAPER_POST_CUT_IRAS,
        "post_cut_akari": PAPER_POST_CUT_AKARI,
        "chance_pairs": chance,
        "chance_probability": false_alarm_probability(chance),
        "published": {
            "pairs": PUBLISHED_PAIRS,
            "good": PUBLISHED_GOOD,
            "distance_au": [PUBLISHED_DISTANCE_AU.0, PUBLISHED_DISTANCE_AU.1],
            "mass_earth": [PUBLISHED_MASS_EARTH.0, PUBLISHED_MASS_EARTH.1],
            "baseline_years": PUBLISHED_BASELINE_YEARS,
        },
    })
}
