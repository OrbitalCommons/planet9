//! Multi-survey combined exclusion analysis.
//!
//! The rigorous combination is a *per-orbit OR* over one shared reference
//! population: each synthetic Planet Nine drawn from the Brown & Batygin
//! (2021) posterior is pushed through all three survey models (ZTF via
//! `p9_2021_ztf`, DES via `p9_2022_des`, PS1 here), and the combined
//! exclusion is the mean probability that *at least one* survey detects it.
//! Unique contributions follow the papers' ZTF -> DES -> PS1 ordering. The
//! published unique fractions live in `p9_core::analysis::surveys` and are
//! used as a cross-check.

use rand::SeedableRng;
use serde::{Deserialize, Serialize};

use p9_2021_orbit::reference_population::{generate_reference_population, heliocentric_distance};
use p9_2022_des::color_models as des_colors;
use p9_2022_des::survey_model::{DesBand, DesSurvey};
use p9_core::constants::DEG2RAD;
use p9_core::coords::observer::{EarthProvider, EarthState};
use p9_core::types::{OrbitalElements, P9Params};

use crate::detection_pipeline::{
    detection_probability_for_orbit, detection_probability_for_orbit_with_earth,
    mid_survey_earth_state,
};
use crate::survey_model::Ps1Survey;

/// Updated Planet Nine parameter estimates from the combined analysis
/// (Brown, Holman & Batygin 2024, after the 78% exclusion).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct UpdatedParameters {
    /// Semi-major axis in AU (median)
    pub a_median: f64,
    /// Semi-major axis upper uncertainty (+170 AU)
    pub a_upper: f64,
    /// Semi-major axis lower uncertainty (-120 AU)
    pub a_lower: f64,
    /// Mass in Earth masses (median)
    pub mass_earth_median: f64,
    /// Mass upper uncertainty (+2.6 ME)
    pub mass_upper: f64,
    /// Mass lower uncertainty (-1.7 ME)
    pub mass_lower: f64,
    /// Apparent V-band magnitude (median)
    pub v_mag_median: f64,
    /// V magnitude upper uncertainty (+1.1)
    pub v_mag_upper: f64,
    /// V magnitude lower uncertainty (-1.4)
    pub v_mag_lower: f64,
}

impl UpdatedParameters {
    /// Paper's updated parameter estimates after three-survey exclusion:
    /// a = 500 (+170/-120) AU, m = 6.6 (+2.6/-1.7) M_earth,
    /// V = 22.0 (+1.1/-1.4).
    pub fn paper_values() -> Self {
        Self {
            a_median: 500.0,
            a_upper: 170.0,
            a_lower: 120.0,
            mass_earth_median: 6.6,
            mass_upper: 2.6,
            mass_lower: 1.7,
            v_mag_median: 22.0,
            v_mag_upper: 1.1,
            v_mag_lower: 1.4,
        }
    }

    /// Convert median parameters to a P9Params for orbital computations.
    ///
    /// Sources for the angular/shape elements (the medians of the Brown &
    /// Batygin 2021 posterior, which the PS1 paper does not update):
    /// - e: from the post-ZTF/DES median perihelion q ~ 340 AU (Belyakov
    ///   et al. 2022) with this struct's a: e = 1 - 340/a (~0.32 at 500 AU)
    /// - i = 16 deg (BB21 median inclination)
    /// - omega = 150 deg, Omega = 97 deg (BB21 median argument of
    ///   perihelion and ascending node)
    /// - mean anomaly 0 (arbitrary epoch; the posterior is uninformative)
    pub fn to_p9_params(&self) -> P9Params {
        P9Params {
            mass_earth: self.mass_earth_median,
            a: self.a_median,
            e: 1.0 - 340.0 / self.a_median,
            i: 16.0 * DEG2RAD,
            omega: 150.0 * DEG2RAD,
            omega_big: 97.0 * DEG2RAD,
            mean_anomaly: 0.0,
        }
    }
}

/// Combined exclusion fractions from the three-survey analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CombinedExclusion {
    /// Fraction excluded by ZTF alone
    pub ztf_frac: f64,
    /// Additional unique fraction excluded by DES
    pub des_unique: f64,
    /// Additional unique fraction excluded by PS1
    pub ps1_unique: f64,
    /// Total combined exclusion
    pub combined: f64,
}

impl CombinedExclusion {
    /// Published exclusion fractions, from the shared p9-core survey table
    /// (ZTF 56.4%, DES unique 5.0%, PS1 unique 17.1%; the paper rounds the
    /// 78.5% sum to 78%).
    pub fn paper_values() -> Self {
        let table = &p9_core::analysis::surveys::EXCLUSION_FRACTIONS;
        let get = |name: &str| {
            table
                .iter()
                .find(|s| s.survey == name)
                .expect("survey in shared table")
                .unique_fraction
        };
        let (ztf, des, ps1) = (get("ZTF"), get("DES"), get("PS1"));
        Self {
            ztf_frac: ztf,
            des_unique: des,
            ps1_unique: ps1,
            combined: p9_core::analysis::surveys::combined_unique_exclusion(&[ztf, des, ps1]),
        }
    }
}

/// Combine *published unique* exclusion fractions (already disjoint by
/// construction) via the shared p9-core helper.
pub fn compute_combined(ztf_frac: f64, des_unique: f64, ps1_unique: f64) -> CombinedExclusion {
    let combined =
        p9_core::analysis::surveys::combined_unique_exclusion(&[ztf_frac, des_unique, ps1_unique]);
    CombinedExclusion {
        ztf_frac,
        des_unique,
        ps1_unique,
        combined,
    }
}

/// Per-orbit OR combination over a shared seeded reference population
/// (P4): draws `n` synthetic Planet Nines from the BB21 posterior emulator,
/// evaluates each survey's detection probability for every orbit, and
/// accumulates
///
///   ztf_frac   = <p_ZTF>
///   des_unique = <p_DES (1 - p_ZTF)>
///   ps1_unique = <p_PS1 (1 - p_ZTF)(1 - p_DES)>
///   combined   = <1 - (1-p_ZTF)(1-p_DES)(1-p_PS1)>
///
/// so `combined = ztf_frac + des_unique + ps1_unique` exactly, per-orbit
/// (not a declarative addition of stored numbers).
pub fn compute_combined_from_population(n: usize, seed: u64) -> CombinedExclusion {
    run_combined_from_population(n, seed, None)
}

/// Ephemeris-path variant of `compute_combined_from_population`: the DES
/// per-night apparent positions and the PS1 footprint position use real
/// Earth states (light-time + aberration) from `earth` instead of the
/// analytic circular-Earth fallback. ZTF's model has no per-epoch
/// geometry and is unchanged.
pub fn compute_combined_from_population_with_earth(
    n: usize,
    seed: u64,
    earth: &mut impl EarthProvider,
) -> CombinedExclusion {
    let des_nights = DesSurvey::default().observing_epochs_with_earth(earth);
    let ps1_state = mid_survey_earth_state(&Ps1Survey::default(), earth);
    run_combined_from_population(n, seed, Some((&des_nights, &ps1_state)))
}

fn run_combined_from_population(
    n: usize,
    seed: u64,
    earth_geometry: Option<(&[(DesBand, f64, EarthState)], &EarthState)>,
) -> CombinedExclusion {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let population = generate_reference_population(n, &mut rng);

    let ztf = p9_2021_ztf::survey_model::ZtfSurvey::default();
    let des = DesSurvey::default();
    let des_color = des_colors::fiducial();
    let ps1 = Ps1Survey::default();

    let mut sum_ztf = 0.0;
    let mut sum_des_unique = 0.0;
    let mut sum_ps1_unique = 0.0;
    let mut sum_combined = 0.0;

    for obj in &population {
        let elements = OrbitalElements {
            a: obj.a,
            e: obj.e,
            i: obj.i,
            omega: obj.omega,
            omega_big: obj.omega_big,
            mean_anomaly: obj.mean_anomaly,
        };

        let p_ztf = p9_2021_ztf::detection_efficiency::detection_probability_for_orbit(
            &ztf,
            &elements,
            obj.v_magnitude,
        );

        let r_au = heliocentric_distance(&P9Params {
            mass_earth: obj.mass,
            a: obj.a,
            e: obj.e,
            i: obj.i,
            omega: obj.omega,
            omega_big: obj.omega_big,
            mean_anomaly: obj.mean_anomaly,
        });
        let mags = des_color.band_magnitudes(obj.mass, r_au);
        let (p_des, p_ps1) = match earth_geometry {
            Some((des_nights, ps1_state)) => (
                des.detection_probability_for_orbit_with_earth(&elements, &mags, des_nights),
                detection_probability_for_orbit_with_earth(
                    &ps1,
                    &elements,
                    obj.v_magnitude,
                    ps1_state,
                ),
            ),
            None => (
                des.detection_probability_for_orbit(&elements, &mags),
                detection_probability_for_orbit(&ps1, &elements, obj.v_magnitude),
            ),
        };

        sum_ztf += p_ztf;
        sum_des_unique += p_des * (1.0 - p_ztf);
        sum_ps1_unique += p_ps1 * (1.0 - p_ztf) * (1.0 - p_des);
        sum_combined += 1.0 - (1.0 - p_ztf) * (1.0 - p_des) * (1.0 - p_ps1);
    }

    let nf = n as f64;
    CombinedExclusion {
        ztf_frac: sum_ztf / nf,
        des_unique: sum_des_unique / nf,
        ps1_unique: sum_ps1_unique / nf,
        combined: sum_combined / nf,
    }
}

/// Describe the remaining viable parameter space after combined exclusion.
///
/// Returns the fraction of the prior that has not been excluded by any survey.
pub fn remaining_parameter_space(exclusion: &CombinedExclusion) -> f64 {
    1.0 - exclusion.combined
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn paper_exclusion_values_from_shared_table() {
        let ex = CombinedExclusion::paper_values();
        assert_relative_eq!(ex.ztf_frac, 0.564, epsilon = 1e-10);
        assert_relative_eq!(ex.des_unique, 0.050, epsilon = 1e-10);
        assert_relative_eq!(ex.ps1_unique, 0.171, epsilon = 1e-10);
        assert_relative_eq!(ex.combined, 0.785, epsilon = 1e-10);
    }

    #[test]
    fn compute_combined_matches_paper() {
        let ex = compute_combined(0.564, 0.050, 0.171);
        assert_relative_eq!(ex.combined, 0.785, epsilon = 1e-10);
    }

    #[test]
    fn remaining_space_is_complement() {
        let ex = CombinedExclusion::paper_values();
        let remaining = remaining_parameter_space(&ex);
        assert!((remaining - 0.215).abs() < 1e-10);
    }

    #[test]
    fn combined_capped_at_unity() {
        let ex = compute_combined(0.8, 0.2, 0.3);
        assert_relative_eq!(ex.combined, 1.0, epsilon = 1e-10);
    }

    /// X1 regression: the per-orbit OR over the shared seeded reference
    /// population must land near the published fractions of the p9-core
    /// table (ZTF 0.564, DES unique 0.050, PS1 unique 0.171, combined
    /// 0.785).
    #[test]
    fn per_orbit_combination_near_published_fractions() {
        let ex = compute_combined_from_population(4000, 2024);
        let paper = CombinedExclusion::paper_values();
        assert!(
            (ex.ztf_frac - paper.ztf_frac).abs() < 0.08,
            "ZTF {:.3} vs published {:.3}",
            ex.ztf_frac,
            paper.ztf_frac
        );
        assert!(
            (ex.des_unique - paper.des_unique).abs() < 0.04,
            "DES unique {:.3} vs published {:.3}",
            ex.des_unique,
            paper.des_unique
        );
        assert!(
            (ex.ps1_unique - paper.ps1_unique).abs() < 0.08,
            "PS1 unique {:.3} vs published {:.3}",
            ex.ps1_unique,
            paper.ps1_unique
        );
        assert!(
            (ex.combined - paper.combined).abs() < 0.10,
            "combined {:.3} vs published {:.3}",
            ex.combined,
            paper.combined
        );
        // Internal consistency of the per-orbit decomposition.
        let sum = ex.ztf_frac + ex.des_unique + ex.ps1_unique;
        assert!((sum - ex.combined).abs() < 1e-9);
    }

    /// Ephemeris-path counterpart of the X1 regression: the per-orbit OR
    /// with real DE421 Earth geometry (DES per-night apparent positions,
    /// PS1 footprint position) must land near the same published
    /// fractions (the 0.564/0.050/0.171 pins). Kernel-gated: skips when
    /// the starfield cache lacks de421.bsp.
    #[test]
    fn per_orbit_combination_near_published_fractions_with_ephemeris() {
        let mut earth = match p9_core::coords::observer::EphemerisEarth::try_new() {
            Ok(e) => e,
            Err(_) => {
                eprintln!("skipping: no ephemeris kernel in starfield cache");
                return;
            }
        };
        let ex = compute_combined_from_population_with_earth(4000, 2024, &mut earth);
        let paper = CombinedExclusion::paper_values();
        assert!(
            (ex.ztf_frac - paper.ztf_frac).abs() < 0.08,
            "ZTF {:.3} vs published {:.3}",
            ex.ztf_frac,
            paper.ztf_frac
        );
        assert!(
            (ex.des_unique - paper.des_unique).abs() < 0.04,
            "DES unique {:.3} vs published {:.3}",
            ex.des_unique,
            paper.des_unique
        );
        assert!(
            (ex.ps1_unique - paper.ps1_unique).abs() < 0.08,
            "PS1 unique {:.3} vs published {:.3}",
            ex.ps1_unique,
            paper.ps1_unique
        );
        assert!(
            (ex.combined - paper.combined).abs() < 0.10,
            "combined {:.3} vs published {:.3}",
            ex.combined,
            paper.combined
        );
        // The analytic fallback gives nearly identical fractions: the
        // Earth-model difference is sub-0.3 deg of apparent position.
        let analytic = compute_combined_from_population(4000, 2024);
        assert!(
            (ex.combined - analytic.combined).abs() < 0.01,
            "ephemeris {:.3} vs analytic {:.3}",
            ex.combined,
            analytic.combined
        );
    }

    #[test]
    fn updated_parameters_match_paper() {
        let params = UpdatedParameters::paper_values();
        assert_relative_eq!(params.a_median, 500.0, epsilon = 1e-10);
        assert_relative_eq!(params.a_upper, 170.0, epsilon = 1e-10);
        assert_relative_eq!(params.a_lower, 120.0, epsilon = 1e-10);
        assert_relative_eq!(params.mass_earth_median, 6.6, epsilon = 1e-10);
        assert_relative_eq!(params.mass_upper, 2.6, epsilon = 1e-10);
        assert_relative_eq!(params.mass_lower, 1.7, epsilon = 1e-10);
        assert_relative_eq!(params.v_mag_median, 22.0, epsilon = 1e-10);
        assert_relative_eq!(params.v_mag_upper, 1.1, epsilon = 1e-10);
        assert_relative_eq!(params.v_mag_lower, 1.4, epsilon = 1e-10);
    }

    #[test]
    fn updated_params_to_p9_params() {
        let params = UpdatedParameters::paper_values();
        let p9 = params.to_p9_params();
        assert_relative_eq!(p9.mass_earth, 6.6, epsilon = 1e-10);
        assert_relative_eq!(p9.a, 500.0, epsilon = 1e-10);
        // e from the post-ZTF/DES median perihelion 340 AU.
        assert_relative_eq!(p9.e, 0.32, epsilon = 1e-10);
        // BB21 median node, not the previously hard-coded 100 deg.
        assert_relative_eq!(p9.omega_big, 97.0 * DEG2RAD, epsilon = 1e-10);
    }
}
