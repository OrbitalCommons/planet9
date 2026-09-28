//! Film export for `p9-2025-stacking`: the numbers its scene and ledger entry draw.

use p9_2021_ztf::detection_efficiency::detection_probability_for_orbit as ztf_probability;
use p9_2021_ztf::survey_model::ZtfSurvey;
use p9_2022_des::color_models::fiducial;
use p9_2022_des::survey_model::DesSurvey;
use p9_2024_panstarrs::detection_pipeline::{
    detection_probability_for_orbit, equatorial_position_deg,
};
use p9_2024_panstarrs::survey_model::Ps1Survey;
use p9_2025_stacking::published::{
    STACKED_DEPTH_REACH, TRIAL_ORBITS_ORDER, ZTF_BASELINE_YEARS, ZTF_SINGLE_DEPTH,
};
use p9_2025_stacking::significance::{
    asymptotic_threshold_sigma, look_elsewhere_threshold_sigma, net_depth_gain_mag,
    trials_penalty_mag,
};
use p9_core::analysis::stacking::matched_filter::{
    frames_to_reach_depth, stack_depth_gain_mag, stacked_limiting_mag,
};
use p9_core::analysis::stacking::orbit_metric::{
    n_trial_orbits, rate_resolution, snr_retained_exact, snr_retained_quadratic,
};
use p9_core::constants::YEAR_DAYS;
use p9_core::data::reference_population::{generate_reference_population, heliocentric_distance};
use p9_core::types::{OrbitalElements, P9Params};
use rand::SeedableRng;
use serde_json::{Value, json};

/// Images stacked per field: "thousands of images over six years".
const STACK_FRAMES: usize = 2000;
/// Point-spread function width (1 sigma, arcsec) and the span of sky rates
/// searched (arcsec/day) for the trial-orbit count.
const PSF_ARCSEC: f64 = 1.0;
const RATE_RANGE: f64 = 40.0;
/// Fraction of the stacked signal a trial orbit must retain.
const SNR_TOLERANCE: f64 = 0.9;
/// Single-trial threshold the trials penalty is measured against, and the
/// threshold the paper quotes for Planet Nine candidates.
const BASELINE_SIGMA: f64 = 5.0;
const P9_SIGMA: f64 = 10.0;
/// Global false-alarm probability held fixed over the whole search.
const GLOBAL_ALPHA: f64 = 0.01;
/// Rubin single-exposure depth the paper compares the stack with.
const RUBIN_SINGLE_DEPTH: f64 = 24.5;
/// Population scored for the Planet Nine reach (same draw as the searches).
const N_POP: usize = 3000;
const SEED: u64 = 2024;

/// Geringer-Sameth et al. (2025), Fig. 6: share of the reference population
/// already ruled out, and the share ZTF stacking would newly reach.
const PUBLISHED_ALREADY: f64 = 0.78;
const PUBLISHED_NEW: f64 = 0.18;

pub fn export() -> Value {
    let baseline_days = ZTF_BASELINE_YEARS * YEAR_DAYS;

    // 1. depth against the number of frames stacked
    let frames: Vec<f64> = (0..=60).map(|k| 10f64.powf(0.1 * k as f64)).collect();
    let depth: Vec<f64> = frames
        .iter()
        .map(|&n| stacked_limiting_mag(ZTF_SINGLE_DEPTH, n.round() as usize))
        .collect();

    // 2. how fast a wrong trial orbit loses the signal
    let cell = rate_resolution(PSF_ARCSEC, baseline_days, SNR_TOLERANCE);
    let rate_error: Vec<f64> = (0..=60).map(|k| 4.0 * cell * k as f64 / 60.0).collect();
    let retained_exact: Vec<f64> = rate_error
        .iter()
        .map(|&dv| snr_retained_exact(PSF_ARCSEC, baseline_days, dv, 2001))
        .collect();
    let retained_quadratic: Vec<f64> = rate_error
        .iter()
        .map(|&dv| snr_retained_quadratic(PSF_ARCSEC, baseline_days, dv))
        .collect();

    // 3. trial orbits against baseline, and what that many trials cost
    let baselines: Vec<f64> = (0..=40)
        .map(|k| 10f64.powf(-0.5 + (baseline_days.log10() + 0.5) * k as f64 / 40.0))
        .collect();
    let trials: Vec<f64> = baselines
        .iter()
        .map(|&t| n_trial_orbits(RATE_RANGE, PSF_ARCSEC, t, SNR_TOLERANCE))
        .collect();
    let n_trials = n_trial_orbits(RATE_RANGE, PSF_ARCSEC, baseline_days, SNR_TOLERANCE);
    let log_trials: Vec<f64> = (0..=48).map(|k| 0.25 * k as f64).collect();
    let threshold: Vec<f64> = log_trials
        .iter()
        .map(|&l| look_elsewhere_threshold_sigma(10f64.powf(l), GLOBAL_ALPHA))
        .collect();
    let z = look_elsewhere_threshold_sigma(n_trials, GLOBAL_ALPHA);
    let gain = stack_depth_gain_mag(STACK_FRAMES);
    let penalty = trials_penalty_mag(z, BASELINE_SIGMA);
    let net = net_depth_gain_mag(STACK_FRAMES, n_trials, GLOBAL_ALPHA, BASELINE_SIGMA);

    // 4. what that depth means for Planet Nine: the members no search has
    //    reached that a ZTF stack would see at the paper's 10-sigma threshold
    let p9_depth = ZTF_SINGLE_DEPTH + gain - trials_penalty_mag(P9_SIGMA, BASELINE_SIGMA);
    let ztf = ZtfSurvey::default();
    let des = DesSurvey::default();
    let des_color = fiducial();
    let ps1 = Ps1Survey::default();
    let mut rng = rand::rngs::StdRng::seed_from_u64(SEED);
    let population = generate_reference_population(N_POP, &mut rng);
    let (mut already, mut newly) = (0.0, 0.0);
    let mut v_mag = Vec::with_capacity(N_POP);
    let mut p_found = Vec::with_capacity(N_POP);
    let mut in_reach = Vec::with_capacity(N_POP);
    for obj in &population {
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
        let p = 1.0 - (1.0 - p_ztf) * (1.0 - p_des) * (1.0 - p_ps1);
        let (_, dec) = equatorial_position_deg(&elements);
        let reach = dec > ztf.dec_limit_deg && obj.v_magnitude < p9_depth;
        already += p;
        if reach {
            newly += 1.0 - p;
        }
        v_mag.push(obj.v_magnitude);
        p_found.push(p);
        in_reach.push(reach);
    }
    let n = N_POP as f64;

    json!({
        "single_depth": ZTF_SINGLE_DEPTH,
        "rubin_single_depth": RUBIN_SINGLE_DEPTH,
        "stack_frames": STACK_FRAMES,
        "baseline_years": ZTF_BASELINE_YEARS,
        "depth": {"frames": frames, "mag": depth},
        "frames_to_rubin": frames_to_reach_depth(ZTF_SINGLE_DEPTH, RUBIN_SINGLE_DEPTH),
        "frames_to_27": frames_to_reach_depth(ZTF_SINGLE_DEPTH, STACKED_DEPTH_REACH),
        "gain_mag": gain,
        "stacked_depth": stacked_limiting_mag(ZTF_SINGLE_DEPTH, STACK_FRAMES),
        "mismatch": {
            "psf_arcsec": PSF_ARCSEC,
            "tolerance": SNR_TOLERANCE,
            "cell_arcsec_day": cell,
            "cell_mas_day": 1000.0 * cell,
            "rate_error": rate_error,
            "retained_exact": retained_exact,
            "retained_quadratic": retained_quadratic,
        },
        "trials": {
            "rate_range": RATE_RANGE,
            "baseline_days": baselines,
            "n": trials,
            "log10_n": log_trials,
            "threshold_sigma": threshold,
        },
        "n_trials": n_trials,
        "log10_trials": n_trials.log10(),
        "threshold_sigma": z,
        "threshold_asymptotic": asymptotic_threshold_sigma(n_trials),
        "baseline_sigma": BASELINE_SIGMA,
        "global_alpha": GLOBAL_ALPHA,
        "penalty_mag": penalty,
        "net_gain_mag": net,
        "net_depth": ZTF_SINGLE_DEPTH + net,
        "p9": {
            "sigma": P9_SIGMA,
            "depth": p9_depth,
            "dec_limit_deg": ztf.dec_limit_deg,
            "already": already / n,
            "new": newly / n,
            "remaining": 1.0 - (already + newly) / n,
            "v_mag": v_mag,
            "p_found": p_found,
            "in_reach": in_reach,
        },
        "p9_new": newly / n,
        "published": {
            "trial_orbits": TRIAL_ORBITS_ORDER,
            "stacked_depth": STACKED_DEPTH_REACH,
            "already": PUBLISHED_ALREADY,
            "new": PUBLISHED_NEW,
        },
    })
}
