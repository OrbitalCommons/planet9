//! Headline reproduction of Eldadi & Loeb (2026), "Applying the Loeb–Turner
//! α-Slope Test to Archival Photometry of TNOs" (arXiv:2605.17098).
//!
//! Two claims are checked:
//!
//! 1. **Physics.** Synthetic photometry generated from p9-core's canonical flux
//!    laws recovers the right power-law index when regressed by [`fit_alpha`]:
//!    reflected sunlight (`reflected_flux_jy`, two inverse-square legs) gives
//!    α ≈ −4 and classifies [`SlopeClass::Reflected`]; thermal emission
//!    (`thermal_flux_jy`, one leg) gives α ≈ −2 and classifies
//!    [`SlopeClass::SelfLuminous`]. The α values are *computed by regression*,
//!    not hard-coded.
//!
//!    Computed values (this build, TNO distances 35–120 AU):
//!      reflected: α ≈ −4.04   thermal: α ≈ −2.00
//!
//! 2. **The funnel.** The published Q1–Q6 bin counts (8557 → 1089 → 186 →
//!    {53, 24, 109}) are internally consistent and match the paper, and Pluto's
//!    22 bins recover the reflected slope in none of them.

use p9_2026_alpha_slope::{
    classify_alpha, fit_alpha, PhotometryPoint, SlopeClass, ANOMALOUS, CANDIDATE_BINS, PASS_Q1_Q3,
    PASS_Q4_Q6, PLUTO_BINS, PLUTO_REFLECTED_RECOVERIES, PUBLISHED_FUNNEL, REFLECTED, SELF_LUMINOUS,
};
use p9_core::analysis::thermal::{reflected_flux_jy, thermal_flux_jy};

/// Heliocentric distances (AU) of the synthetic epochs — a TNO-scale lever arm.
const DISTANCES_AU: [f64; 6] = [35.0, 45.0, 60.0, 80.0, 100.0, 120.0];

#[test]
fn synthetic_reflected_recovers_alpha_minus_four() {
    // Reflected sunlight from a fixed body, observed near opposition, across a
    // range of heliocentric distances. p9-core's reflected_flux_jy carries both
    // inverse-square legs, so the flux falls ~d^-4.
    let albedo = 0.1;
    let radius_m = 1.0e6;
    let nu_hz = 2.997_924_58e8 / 0.55e-6; // V band
    let points: Vec<PhotometryPoint> = DISTANCES_AU
        .iter()
        .map(|&d| PhotometryPoint::from_flux(d, reflected_flux_jy(albedo, radius_m, d, nu_hz)))
        .collect();

    let fit = fit_alpha(&points).expect("reflected fit");
    assert!(
        (fit.alpha - (-4.0)).abs() < 0.1,
        "reflected α = {:.4}, want within 0.1 of -4 (R²={:.5})",
        fit.alpha,
        fit.r_squared
    );
    assert_eq!(classify_alpha(fit.alpha), SlopeClass::Reflected);
}

#[test]
fn synthetic_thermal_recovers_alpha_minus_two() {
    // A self-luminous body of fixed temperature and radius: only the body→Earth
    // inverse-square leg survives, so the flux falls ~d^-2.
    let t_k = 40.0;
    let radius_m = 1.0e6;
    let nu_hz = 150e9; // mm-band, deep in Rayleigh-Jeans for a 40 K body
    let points: Vec<PhotometryPoint> = DISTANCES_AU
        .iter()
        .map(|&d| PhotometryPoint::from_flux(d, thermal_flux_jy(t_k, radius_m, d, nu_hz)))
        .collect();

    let fit = fit_alpha(&points).expect("thermal fit");
    assert!(
        (fit.alpha - (-2.0)).abs() < 0.1,
        "thermal α = {:.4}, want within 0.1 of -2 (R²={:.5})",
        fit.alpha,
        fit.r_squared
    );
    assert_eq!(classify_alpha(fit.alpha), SlopeClass::SelfLuminous);
}

#[test]
fn magnitude_built_reflected_also_recovers_minus_four() {
    // The same reflected flux fed in as apparent magnitudes (f ∝ 10^-0.4m)
    // must recover the same slope — the magnitude path and the flux path agree.
    let albedo = 0.1;
    let radius_m = 1.0e6;
    let nu_hz = 2.997_924_58e8 / 0.55e-6;
    let points: Vec<PhotometryPoint> = DISTANCES_AU
        .iter()
        .map(|&d| {
            let f = reflected_flux_jy(albedo, radius_m, d, nu_hz);
            // Arbitrary zero point; cancels in the regression.
            let m = -2.5 * f.log10();
            PhotometryPoint::from_magnitude(d, m)
        })
        .collect();

    let fit = fit_alpha(&points).expect("reflected-from-mag fit");
    assert!(
        (fit.alpha - (-4.0)).abs() < 0.1,
        "reflected-from-magnitude α = {:.4}",
        fit.alpha
    );
    assert_eq!(classify_alpha(fit.alpha), SlopeClass::Reflected);
}

#[test]
fn funnel_counts_match_published_and_are_consistent() {
    assert!(
        PUBLISHED_FUNNEL.is_consistent(),
        "published funnel failed internal consistency"
    );
    assert_eq!(CANDIDATE_BINS, 8557);
    assert_eq!(PASS_Q1_Q3, 1089);
    assert_eq!(PASS_Q4_Q6, 186);
    assert_eq!(REFLECTED, 53);
    assert_eq!(SELF_LUMINOUS, 24);
    assert_eq!(ANOMALOUS, 109);
    // The three α-classes partition the Q4–Q6 survivors.
    assert_eq!(REFLECTED + SELF_LUMINOUS + ANOMALOUS, PASS_Q4_Q6);
    // Each cut is a subset of the previous stage.
    assert!(PASS_Q4_Q6 <= PASS_Q1_Q3 && PASS_Q1_Q3 <= CANDIDATE_BINS);
}

#[test]
fn pluto_recovers_no_reflected_slope() {
    // Even the brightest TNO: 22 bins, none recovering α ≈ -4 under a single
    // instrument+band.
    assert_eq!(PLUTO_BINS, 22);
    assert_eq!(PLUTO_REFLECTED_RECOVERIES, 0);
}
