//! Parameter space exploration from Bailey+ (2016) Section 3.
//!
//! Surveys (a₉, e₉) space for each mass to find the inclination i₉
//! required to produce the observed 6° solar obliquity.

use p9_core::constants::*;

use crate::secular_hamiltonian::{integrate_obliquity, SecularParams, SpinOrbitState};
use crate::solar_model;

/// Result of a parameter survey point.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SurveyResult {
    pub mass_earth: f64,
    pub a9: f64,
    pub e9: f64,
    /// Required i₉ to produce ~6° obliquity (rad), or None if not found
    pub required_i9: Option<f64>,
    /// Final obliquity achieved (deg)
    pub final_obliquity_deg: f64,
}

/// Find the Planet Nine inclination that produces a target obliquity.
///
/// Uses bisection search over i₉ ∈ [5°, 50°].
pub fn find_required_inclination(
    mass_earth: f64,
    a9: f64,
    e9: f64,
    target_obliquity_deg: f64,
    tolerance_deg: f64,
) -> Option<(f64, f64)> {
    let q = a9 * (1.0 - e9);
    if q < 150.0 || q > 350.0 {
        return None;
    }

    let m9_solar = mass_earth * EARTH_MASS_SOLAR;
    let mut i_low = 5.0 * DEG2RAD;
    let mut i_high = 50.0 * DEG2RAD;

    let obl_low = run_obliquity(m9_solar, a9, e9, i_low);
    let obl_high = run_obliquity(m9_solar, a9, e9, i_high);

    let target = target_obliquity_deg;
    if (obl_low - target) * (obl_high - target) > 0.0 {
        if (obl_low - target).abs() < (obl_high - target).abs() {
            return Some((i_low, obl_low));
        } else {
            return Some((i_high, obl_high));
        }
    }

    for _ in 0..20 {
        let i_mid = 0.5 * (i_low + i_high);
        let obl_mid = run_obliquity(m9_solar, a9, e9, i_mid);

        if (obl_mid - target).abs() < tolerance_deg {
            return Some((i_mid, obl_mid));
        }

        if (obl_mid - target) * (obl_low - target) < 0.0 {
            i_high = i_mid;
        } else {
            i_low = i_mid;
        }
    }

    let i_mid = 0.5 * (i_low + i_high);
    let obl_mid = run_obliquity(m9_solar, a9, e9, i_mid);
    Some((i_mid, obl_mid))
}

/// Create an initial state for given P9 parameters.
fn make_initial_state(m9_solar: f64, a9: f64, e9: f64, i9: f64) -> SpinOrbitState {
    let l_gp_mag = solar_model::giant_planet_angular_momentum();
    let l_9_mag = solar_model::p9_angular_momentum(m9_solar, a9, e9);
    SpinOrbitState::from_inclinations(i9, std::f64::consts::PI, l_gp_mag, l_9_mag)
}

/// Run a single obliquity integration and return final obliquity in degrees.
fn run_obliquity(m9_solar: f64, a9: f64, e9: f64, i9: f64) -> f64 {
    let initial = make_initial_state(m9_solar, a9, e9, i9);

    let params = SecularParams {
        m9_solar,
        a9,
        e9,
        t_total: 4.5 * GYR_DAYS,
        dt: 5e4 * YEAR_DAYS,
    };

    let snapshots = integrate_obliquity(initial, &params, 4.5 * GYR_DAYS);
    let final_snap = snapshots.last().unwrap();
    final_snap.obliquity * RAD2DEG
}

/// Run the full parameter survey for a given mass.
pub fn parameter_survey(mass_earth: f64) -> Vec<SurveyResult> {
    let a_values: Vec<f64> = (4..=9).map(|i| i as f64 * 100.0).collect();
    let e_values: Vec<f64> = (3..=9).map(|i| i as f64 * 0.1).collect();

    let mut results = Vec::new();

    for &a9 in &a_values {
        for &e9 in &e_values {
            let q = a9 * (1.0 - e9);
            if q < 150.0 || q > 350.0 {
                continue;
            }

            let search = find_required_inclination(mass_earth, a9, e9, 6.0, 0.5);

            let (required_i9, final_obl) = match search {
                Some((i9, obl)) => (Some(i9), obl),
                None => (None, 0.0),
            };

            results.push(SurveyResult {
                mass_earth,
                a9,
                e9,
                required_i9,
                final_obliquity_deg: final_obl,
            });
        }
    }

    results
}

/// Quick survey for testing (fewer points, shorter integration).
pub fn quick_survey() -> Vec<SurveyResult> {
    let m9_solar = 10.0 * EARTH_MASS_SOLAR;

    let test_cases = [(700.0, 0.6), (500.0, 0.5)];

    test_cases
        .iter()
        .map(|&(a9, e9)| {
            let initial = make_initial_state(m9_solar, a9, e9, 20.0 * DEG2RAD);

            let params = SecularParams {
                m9_solar,
                a9,
                e9,
                t_total: 1e8 * YEAR_DAYS,
                dt: 1e5 * YEAR_DAYS,
            };

            let snapshots = integrate_obliquity(initial, &params, 1e8 * YEAR_DAYS);
            let final_snap = snapshots.last().unwrap();

            SurveyResult {
                mass_earth: 10.0,
                a9,
                e9,
                required_i9: Some(20.0 * DEG2RAD),
                final_obliquity_deg: final_snap.obliquity * RAD2DEG,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quick_survey_runs() {
        let results = quick_survey();
        assert_eq!(results.len(), 2);

        for r in &results {
            assert!(r.final_obliquity_deg >= 0.0);
            assert!(r.final_obliquity_deg.is_finite());
        }
    }

    #[test]
    fn test_run_obliquity_monotonic_in_inclination() {
        let m9_solar = 10.0 * EARTH_MASS_SOLAR;

        let obl_10 = run_obliquity(m9_solar, 700.0, 0.6, 10.0 * DEG2RAD);
        let obl_30 = run_obliquity(m9_solar, 700.0, 0.6, 30.0 * DEG2RAD);

        assert!(
            obl_30 > obl_10,
            "30° P9 inclination should produce more obliquity ({:.2}°) than 10° ({:.2}°)",
            obl_30,
            obl_10
        );
    }
}
