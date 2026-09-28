//! Film export for `p9-2026-alpha-slope`: the numbers its scene and ledger entry draw.

use p9_2026_alpha_slope::funnel::{
    ANOMALOUS, CANDIDATE_BINS, PASS_Q1_Q3, PASS_Q4_Q6, PLUTO_BINS, PLUTO_REFLECTED_RECOVERIES,
    PUBLISHED_FUNNEL, REFLECTED, SELF_LUMINOUS, SELF_LUMINOUS_ALL_PANSTARRS,
};
use p9_2026_alpha_slope::slope::{
    ALPHA_REFLECTED, ALPHA_SELF_LUMINOUS, ALPHA_TOLERANCE, PhotometryPoint, SlopeClass,
    classify_alpha, fit_alpha,
};
use p9_core::analysis::thermal::{C_LIGHT, reflected_flux_jy, thermal_flux_jy};
use serde_json::{Value, json};

/// Candidate bins quoted in the abstract and Table 2 of Eldadi & Loeb (2026).
const PUBLISHED_CANDIDATE_BINS: u32 = 8606;
/// Numbered TNOs the census covers (Table 2).
const PUBLISHED_TNOS: u32 = 913;
/// Pluto's heliocentric distance since its discovery: perihelion in 1989 and
/// the 1930 discovery distance (AU).
const PLUTO_RANGE_AU: (f64, f64) = (29.7, 41.0);
/// Number of synthetic epochs per fit.
const N_EPOCHS: usize = 24;
/// Body used for the synthetic photometry: radius (m), albedo, temperature (K).
const RADIUS_M: f64 = 1.0e6;
const ALBEDO: f64 = 0.1;
const TEMPERATURE_K: f64 = 40.0;

fn epochs(d_min: f64, d_max: f64) -> Vec<f64> {
    (0..N_EPOCHS)
        .map(|k| d_min + (d_max - d_min) * k as f64 / (N_EPOCHS - 1) as f64)
        .collect()
}

/// Fitted slope of reflected-light photometry between `d_min` and `d_max` when
/// the magnitude zero point drifts by `drift_mag` across the span.
fn alpha_with_drift(d_min: f64, d_max: f64, drift_mag: f64) -> f64 {
    let nu_v = C_LIGHT / 0.55e-6;
    let points: Vec<PhotometryPoint> = epochs(d_min, d_max)
        .into_iter()
        .enumerate()
        .map(|(k, d)| {
            let flux = reflected_flux_jy(ALBEDO, RADIUS_M, d, nu_v);
            let offset = drift_mag * k as f64 / (N_EPOCHS - 1) as f64;
            PhotometryPoint::from_magnitude(d, -2.5 * flux.log10() + offset)
        })
        .collect();
    fit_alpha(&points).map_or(f64::NAN, |fit| fit.alpha)
}

fn class_name(alpha: f64) -> &'static str {
    match classify_alpha(alpha) {
        SlopeClass::Reflected => "reflected",
        SlopeClass::SelfLuminous => "self-luminous",
        SlopeClass::Anomalous => "anomalous",
    }
}

pub fn export() -> Value {
    let nu_v = C_LIGHT / 0.55e-6;
    let nu_mm = 150.0e9;

    // 1. the two laws, and the slopes the crate's regression recovers
    let distances: Vec<f64> = (0..=45).map(|k| 30.0 + 2.0 * k as f64).collect();
    let reflected: Vec<f64> = distances
        .iter()
        .map(|&d| reflected_flux_jy(ALBEDO, RADIUS_M, d, nu_v))
        .collect();
    let thermal: Vec<f64> = distances
        .iter()
        .map(|&d| thermal_flux_jy(TEMPERATURE_K, RADIUS_M, d, nu_mm))
        .collect();
    let fit_of = |flux: &[f64]| {
        let points: Vec<PhotometryPoint> = distances
            .iter()
            .zip(flux)
            .map(|(&d, &f)| PhotometryPoint::from_flux(d, f))
            .collect();
        fit_alpha(&points).expect("slope fit")
    };
    let fit_reflected = fit_of(&reflected);
    let fit_thermal = fit_of(&thermal);

    // 2. why the archive struggles: a TNO barely changes distance, so a small
    //    calibration drift swamps the slope
    let drifts: Vec<f64> = (0..=40).map(|k| -1.0 + 0.05 * k as f64).collect();
    let pluto_alpha: Vec<f64> = drifts
        .iter()
        .map(|&m| alpha_with_drift(PLUTO_RANGE_AU.0, PLUTO_RANGE_AU.1, m))
        .collect();
    let pluto_clean = alpha_with_drift(PLUTO_RANGE_AU.0, PLUTO_RANGE_AU.1, 0.0);
    let per_tenth = alpha_with_drift(PLUTO_RANGE_AU.0, PLUTO_RANGE_AU.1, 0.1) - pluto_clean;
    // Drift that carries the fitted slope from the reflected to the
    // self-luminous window, by bisection (the slope rises with drift).
    let edge = ALPHA_SELF_LUMINOUS - ALPHA_TOLERANCE;
    let (mut lo, mut hi) = (-2.0_f64, 2.0_f64);
    let rising = alpha_with_drift(PLUTO_RANGE_AU.0, PLUTO_RANGE_AU.1, hi) > pluto_clean;
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        let above = alpha_with_drift(PLUTO_RANGE_AU.0, PLUTO_RANGE_AU.1, mid) > edge;
        if above == rising {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let flip = 0.5 * (lo + hi);

    let spans: Vec<f64> = (0..=30).map(|k| 1.02 + 0.02 * k as f64).collect();
    let span_error: Vec<f64> = spans
        .iter()
        .map(|&s| {
            (alpha_with_drift(40.0, 40.0 * s, 0.1) - alpha_with_drift(40.0, 40.0 * s, 0.0)).abs()
        })
        .collect();

    let funnel = PUBLISHED_FUNNEL;
    json!({
        "alpha_reflected": fit_reflected.alpha,
        "alpha_thermal": fit_thermal.alpha,
        "class_reflected": class_name(fit_reflected.alpha),
        "class_thermal": class_name(fit_thermal.alpha),
        "theory": {
            "reflected": ALPHA_REFLECTED,
            "self_luminous": ALPHA_SELF_LUMINOUS,
            "tolerance": ALPHA_TOLERANCE,
        },
        "laws": {
            "distance_au": distances,
            "reflected": reflected.iter().map(|f| f / reflected[0]).collect::<Vec<_>>(),
            "thermal": thermal.iter().map(|f| f / thermal[0]).collect::<Vec<_>>(),
        },
        "pluto": {
            "range_au": [PLUTO_RANGE_AU.0, PLUTO_RANGE_AU.1],
            "alpha_clean": pluto_clean,
            "drift_mag": drifts,
            "alpha": pluto_alpha,
            "alpha_per_tenth_mag": per_tenth,
            "drift_to_flip_mag": flip,
            "bins": PLUTO_BINS,
            "recoveries": PLUTO_REFLECTED_RECOVERIES,
        },
        "lever": {"span": spans, "alpha_error_per_tenth_mag": span_error},
        "funnel": {
            "candidate_bins": PUBLISHED_CANDIDATE_BINS,
            "crate_candidate_bins": CANDIDATE_BINS,
            "pass_q1_q3": PASS_Q1_Q3,
            "pass_q4_q6": PASS_Q4_Q6,
            "reflected": REFLECTED,
            "self_luminous": SELF_LUMINOUS,
            "anomalous": ANOMALOUS,
            "consistent": funnel.is_consistent(),
            "tnos": PUBLISHED_TNOS,
            "self_luminous_all_panstarrs": SELF_LUMINOUS_ALL_PANSTARRS,
        },
    })
}
