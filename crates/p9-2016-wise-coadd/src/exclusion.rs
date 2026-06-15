//! Parameter-space exclusion from the coadd non-detection.
//!
//! Meisner et al. (2016) detect no moving Planet Nine in their deep W1 coadds.
//! As in the sibling 3π-search crate, the fraction of a reference P9 population
//! the coadd *would* have detected is the fraction the non-detection excludes.
//! The per-object detection probability is the product of
//!   * a logistic magnitude-efficiency centered on the **coadd** depth (the
//!     `coadd` √N-deepened limit), evaluated at the object's band magnitude
//!     (`band`), and
//!   * the galactic-plane footprint mask, reused verbatim from
//!     `p9-2018-wise-search` (`WiseSurvey`).
//!
//! We reuse the sibling's [`WiseSurvey`] by constructing it with the coadd
//! depth as its `w1_depth`; its logistic `detection_efficiency` and
//! `in_footprint`/`sky_coverage_fraction` then apply unchanged.

use serde::{Deserialize, Serialize};

use p9_2018_wise_search::sky::galactic_latitude_deg;
use p9_2018_wise_search::survey_model::WiseSurvey;
use p9_core::types::OrbitalElements;

use crate::band::{band_magnitude, WiseBand};
use crate::coadd::{coadd_depth, single_frame_depth};

/// One reference Planet Nine realization (orbit, mass, current distance).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ReferenceP9 {
    pub elements: OrbitalElements,
    pub mass_earth: f64,
    pub distance_au: f64,
}

/// Result of the coadd exclusion analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExclusionResult {
    /// Expected fraction of the reference population the coadd would detect.
    pub fraction_excluded: f64,
    /// Number of realizations evaluated.
    pub n_total: u32,
    /// Expected number detected (sum of per-object probabilities).
    pub expected_excluded: f64,
}

/// A `WiseSurvey` whose depth is the coadd depth in `band` after `n_frames`,
/// reusing the sibling's logistic/footprint machinery.
pub fn coadd_survey(band: WiseBand, n_frames: u32) -> WiseSurvey {
    WiseSurvey {
        w1_depth: coadd_depth(band, n_frames),
        ..WiseSurvey::default()
    }
}

/// A `WiseSurvey` at the single-frame depth in `band` (for the coadd-vs-single
/// contrast).
pub fn single_frame_survey(band: WiseBand) -> WiseSurvey {
    WiseSurvey {
        w1_depth: single_frame_depth(band),
        ..WiseSurvey::default()
    }
}

/// Expected exclusion fraction of `population` in `band` for a survey whose
/// depth is `survey.w1_depth` (use [`coadd_survey`] or [`single_frame_survey`]).
/// Deterministic: accumulates per-object probabilities, no RNG inside.
pub fn compute_exclusion(
    survey: &WiseSurvey,
    band: WiseBand,
    population: &[ReferenceP9],
) -> ExclusionResult {
    let n_total = population.len() as u32;
    let expected_excluded: f64 = population
        .iter()
        .map(|p| {
            let mag = band_magnitude(p.mass_earth, p.distance_au, band);
            let footprint = if survey.in_footprint(galactic_latitude_deg(&p.elements)) {
                1.0
            } else {
                0.0
            };
            survey.detection_efficiency(mag) * footprint
        })
        .sum();

    let fraction_excluded = if n_total > 0 {
        expected_excluded / n_total as f64
    } else {
        0.0
    };

    ExclusionResult {
        fraction_excluded,
        n_total,
        expected_excluded,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::band::{W1, W2};
    use p9_core::constants::DEG2RAD;

    /// A seeded synthetic reference population spread over mass/distance and
    /// orbital angles, built without external crates (self-contained LCG so the
    /// crate has no dev-dependency beyond approx). Distances span 250–900 AU,
    /// masses 5–15 M⊕.
    fn population(n: usize, seed: u64) -> Vec<ReferenceP9> {
        let mut state = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64) / ((1u64 << 53) as f64)
        };
        (0..n)
            .map(|_| {
                let mass = 5.0 + 10.0 * next();
                let distance = 250.0 + 650.0 * next();
                let i = (10.0 + 20.0 * next()) * DEG2RAD;
                let omega = 360.0 * next() * DEG2RAD;
                let omega_big = 360.0 * next() * DEG2RAD;
                let m = 360.0 * next() * DEG2RAD;
                ReferenceP9 {
                    elements: OrbitalElements {
                        a: distance,
                        e: 0.2,
                        i,
                        omega,
                        omega_big,
                        mean_anomaly: m,
                    },
                    mass_earth: mass,
                    distance_au: distance,
                }
            })
            .collect()
    }

    #[test]
    fn exclusion_fraction_is_a_probability() {
        let pop = population(3000, 2016);
        let survey = coadd_survey(W1, 100);
        let r = compute_exclusion(&survey, W1, &pop);
        assert!(
            r.fraction_excluded > 0.0 && r.fraction_excluded < 1.0,
            "fraction = {}",
            r.fraction_excluded
        );
    }

    #[test]
    fn exclusion_increases_with_assumed_mass() {
        let base = population(3000, 7);
        let survey = coadd_survey(W1, 100);
        let light: Vec<ReferenceP9> = base
            .iter()
            .map(|p| ReferenceP9 {
                mass_earth: 5.0,
                ..*p
            })
            .collect();
        let heavy: Vec<ReferenceP9> = base
            .iter()
            .map(|p| ReferenceP9 {
                mass_earth: 18.0,
                ..*p
            })
            .collect();
        let f_light = compute_exclusion(&survey, W1, &light).fraction_excluded;
        let f_heavy = compute_exclusion(&survey, W1, &heavy).fraction_excluded;
        assert!(
            f_heavy > f_light + 0.005,
            "heavy {f_heavy:.4} should exceed light {f_light:.4}"
        );
    }

    #[test]
    fn coadd_excludes_more_than_single_frame() {
        // The deeper coadd detects (excludes) a larger fraction than single
        // frames — the central result of the paper.
        let pop = population(3000, 99);
        let single = single_frame_survey(W1);
        let coadd = coadd_survey(W1, 100);
        let f_single = compute_exclusion(&single, W1, &pop).fraction_excluded;
        let f_coadd = compute_exclusion(&coadd, W1, &pop).fraction_excluded;
        assert!(
            f_coadd > f_single,
            "coadd {f_coadd:.4} should exceed single-frame {f_single:.4}"
        );
    }

    #[test]
    fn w2_excludes_more_than_w1_for_cold_population() {
        // For a cold reference population, W2 (more favorable thermally)
        // excludes a larger fraction than W1 at the same coadd depth number.
        let pop = population(3000, 31);
        let depth = 17.0;
        let survey = WiseSurvey {
            w1_depth: depth,
            ..WiseSurvey::default()
        };
        let f_w1 = compute_exclusion(&survey, W1, &pop).fraction_excluded;
        let f_w2 = compute_exclusion(&survey, W2, &pop).fraction_excluded;
        assert!(
            f_w2 > f_w1,
            "W2 {f_w2:.4} should exceed W1 {f_w1:.4} for a cold population"
        );
    }
}
