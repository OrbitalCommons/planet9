//! Detection analysis for the PS1 Planet Nine search.
//!
//! Maps synthetic Planet Nine orbits to *equatorial* sky positions (P2 fix:
//! the dec > -30 deg footprint cut was previously applied to ecliptic
//! latitude, misclassifying exactly the PS1/DES/ZTF unique-coverage
//! regions), evaluates the cadence/linking detection probability, and
//! accumulates expected counts deterministically. Monte-Carlo realizations
//! thread a caller-provided seed (P1 fix: previously unseeded
//! `rand::random`).

use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use p9_2022_des::sky::{apparent_position_deg, apparent_position_with_earth_deg};
use p9_core::coords::observer::{EarthProvider, EarthState, Timescale};
use p9_core::types::OrbitalElements;

use crate::survey_model::Ps1Survey;

/// Summary of PS1 detection results from the synthetic population analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ps1Result {
    /// Total synthetic orbits detected in PS1 pipeline
    pub n_detected: u32,
    /// Detections unique to PS1 (not recovered by ZTF or DES)
    pub n_unique: u32,
}

impl Ps1Result {
    /// Paper's reported values (Brown, Holman & Batygin 2024): of the
    /// 100,000-member reference population, 69,082 are detected nine or
    /// more times by PS1, of which 17,054 are unique to PS1 (17.1% of the
    /// population, the shared p9-core figure).
    pub fn paper_values() -> Self {
        Self {
            n_detected: 69_082,
            n_unique: 17_054,
        }
    }

    /// Fraction of PS1 detections that are unique to PS1.
    pub fn unique_fraction_of_detected(&self) -> f64 {
        self.n_unique as f64 / self.n_detected as f64
    }
}

/// Equatorial (RA, dec) in degrees of an orbit's current apparent sky
/// position (geocentric, at the orbit's stored epoch), via the shared
/// ecliptic -> equatorial conversion in `p9_2022_des::sky`.
pub fn equatorial_position_deg(elements: &OrbitalElements) -> (f64, f64) {
    let (ra, dec, _) = apparent_position_deg(elements, 0.0);
    (ra, dec)
}

/// Per-orbit PS1 detection probability: equatorial footprint membership
/// plus the magnitude-dependent cadence/linking model.
pub fn detection_probability_for_orbit(
    survey: &Ps1Survey,
    elements: &OrbitalElements,
    v_magnitude: f64,
) -> f64 {
    let (_, dec) = equatorial_position_deg(elements);
    survey.detection_probability(v_magnitude, dec)
}

/// Earth state at the middle of the PS1 observing-epoch grid (real
/// ephemeris epochs from `Ps1Survey::survey_start`/`epoch_offsets_days`
/// instead of raw year floats). Compute once and share across orbits.
pub fn mid_survey_earth_state(survey: &Ps1Survey, earth: &mut impl EarthProvider) -> EarthState {
    let ts = Timescale::default();
    let offsets = survey.epoch_offsets_days();
    let mid = offsets[offsets.len() / 2];
    earth.earth_state(&ts.tt_jd(survey.survey_start(&ts).tt() + mid, None))
}

/// Ephemeris-path counterpart of `equatorial_position_deg`: apparent
/// (light-time + aberration) equatorial position from a real Earth state
/// (see `mid_survey_earth_state`).
pub fn equatorial_position_with_earth_deg(
    elements: &OrbitalElements,
    earth_state: &EarthState,
) -> (f64, f64) {
    let (ra, dec, _) = apparent_position_with_earth_deg(elements, earth_state, 0.0);
    (ra, dec)
}

/// Ephemeris-path counterpart of `detection_probability_for_orbit`.
pub fn detection_probability_for_orbit_with_earth(
    survey: &Ps1Survey,
    elements: &OrbitalElements,
    v_magnitude: f64,
    earth_state: &EarthState,
) -> f64 {
    let (_, dec) = equatorial_position_with_earth_deg(elements, earth_state);
    survey.detection_probability(v_magnitude, dec)
}

/// Deterministic expected number of detections for a synthetic population
/// of `(elements, V magnitude)` pairs.
pub fn expected_detections(survey: &Ps1Survey, orbits: &[(OrbitalElements, f64)]) -> f64 {
    orbits
        .iter()
        .map(|(elements, mag)| detection_probability_for_orbit(survey, elements, *mag))
        .sum()
}

/// Seeded Monte-Carlo realization of the detection count (one Bernoulli per
/// orbit at its detection probability).
pub fn compute_detection(survey: &Ps1Survey, orbits: &[(OrbitalElements, f64)], seed: u64) -> u32 {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    orbits
        .iter()
        .filter(|(elements, mag)| {
            rng.gen::<f64>() < detection_probability_for_orbit(survey, elements, *mag)
        })
        .count() as u32
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    fn orbit(i_deg: f64, omega_big_deg: f64, m_deg: f64) -> OrbitalElements {
        OrbitalElements {
            a: 450.0,
            e: 0.2,
            i: i_deg * DEG2RAD,
            omega: 0.0,
            omega_big: omega_big_deg * DEG2RAD,
            mean_anomaly: m_deg * DEG2RAD,
        }
    }

    #[test]
    fn paper_values_consistent_with_shared_unique_fraction() {
        let result = Ps1Result::paper_values();
        assert_eq!(result.n_detected, 69_082);
        assert_eq!(result.n_unique, 17_054);
        // 17,054 / 100,000 is the shared p9-core PS1 unique fraction.
        let shared = p9_core::analysis::surveys::EXCLUSION_FRACTIONS
            .iter()
            .find(|s| s.survey == "PS1")
            .unwrap()
            .unique_fraction;
        assert!((result.n_unique as f64 / 100_000.0 - shared).abs() < 1e-3);
    }

    #[test]
    fn declination_is_not_ecliptic_latitude() {
        // P2 regression: an i = 0 orbit has ecliptic latitude ~0 everywhere,
        // but its declination swings by +/- the obliquity (~23.4 deg); near
        // ecliptic longitude ~270 deg it sits well south of -20 deg.
        let mut min_dec = f64::INFINITY;
        let mut max_dec = f64::NEG_INFINITY;
        for m in 0..72 {
            let (_, dec) = equatorial_position_deg(&orbit(0.0, 0.0, m as f64 * 5.0));
            min_dec = min_dec.min(dec);
            max_dec = max_dec.max(dec);
        }
        assert!(min_dec < -20.0, "min dec = {min_dec:.1}");
        assert!(max_dec > 20.0, "max dec = {max_dec:.1}");
    }

    #[test]
    fn bright_northern_objects_all_detected() {
        let survey = Ps1Survey::default();
        // Positions chosen on the northern declination side.
        let orbits: Vec<(OrbitalElements, f64)> = (0..100)
            .map(|k| (orbit(10.0, 0.0, 40.0 + k as f64 * 0.5), 18.0))
            .collect();
        let n = compute_detection(&survey, &orbits, 42);
        let expected = expected_detections(&survey, &orbits);
        assert!(n >= 90, "Expected most bright objects detected, got {n}");
        assert!((n as f64 - expected).abs() < 10.0);
    }

    #[test]
    fn faint_objects_all_missed() {
        let survey = Ps1Survey::default();
        let orbits: Vec<(OrbitalElements, f64)> = (0..50)
            .map(|k| (orbit(10.0, 0.0, k as f64), 25.0))
            .collect();
        assert_eq!(compute_detection(&survey, &orbits, 42), 0);
    }

    #[test]
    fn seeded_runs_reproduce() {
        let survey = Ps1Survey::default();
        let orbits: Vec<(OrbitalElements, f64)> = (0..200)
            .map(|k| (orbit(15.0, 30.0, k as f64), 21.4 + (k % 10) as f64 * 0.1))
            .collect();
        let a = compute_detection(&survey, &orbits, 7);
        let b = compute_detection(&survey, &orbits, 7);
        assert_eq!(a, b, "same seed must reproduce");
    }
}
