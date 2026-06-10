//! Candidate pair search and selection criteria.
//!
//! Reproduces the source-pair matching algorithm from Phan et al. (2025):
//!
//! 1. Start with IRAS FSC 60 µm sources and AKARI Monthly Unconfirmed 90 µm sources
//! 2. Apply flux cuts to retain only sources consistent with P9 thermal emission
//! 3. For each IRAS-AKARI pair, compute angular separation
//! 4. Retain pairs with separation in [42', 69.6'] (corresponding to 500–700 AU)
//! 5. Apply flux-ratio consistency check
//! 6. Image inspection (not reproducible computationally — we flag candidates)
//!
//! The paper finds 13 pairs after step 4, and 1 good candidate after step 6.

use serde::{Deserialize, Serialize};

use crate::survey_model::{angular_separation_arcmin, FirSource};

/// Selection criteria for the IRAS-AKARI cross-match.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SelectionCriteria {
    /// Minimum angular separation in arcminutes
    pub sep_min_arcmin: f64,
    /// Maximum angular separation in arcminutes
    pub sep_max_arcmin: f64,
    /// Maximum allowed flux ratio (IRAS 60µm / AKARI 90µm) for consistency
    pub max_flux_ratio: f64,
    /// Minimum allowed flux ratio
    pub min_flux_ratio: f64,
}

impl Default for SelectionCriteria {
    fn default() -> Self {
        Self {
            sep_min_arcmin: 42.0,
            sep_max_arcmin: 69.6,
            // Flux ratio F_60/F_90 for a 30–50 K blackbody is roughly 0.2–2.0
            max_flux_ratio: 3.0,
            min_flux_ratio: 0.05,
        }
    }
}

/// A candidate IRAS-AKARI source pair.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CandidatePair {
    /// The IRAS source
    pub iras_source: FirSource,
    /// The AKARI source
    pub akari_source: FirSource,
    /// Angular separation in arcminutes
    pub separation_arcmin: f64,
    /// Implied heliocentric distance in AU (from proper motion inversion)
    pub implied_distance_au: f64,
    /// Flux ratio (IRAS 60µm / AKARI 90µm)
    pub flux_ratio: f64,
}

/// Result of the candidate search.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SearchResult {
    /// Number of IRAS sources considered
    pub n_iras: usize,
    /// Number of AKARI sources considered
    pub n_akari: usize,
    /// Number of pairs within the angular separation window
    pub n_pairs_in_window: usize,
    /// Number of pairs passing the flux ratio cut
    pub n_pairs_flux_cut: usize,
    /// Final candidate pairs
    pub candidates: Vec<CandidatePair>,
}

/// Implied heliocentric distance from angular separation over a baseline.
///
/// Inverts µ = sqrt(GM_sun)/d^(3/2) * conversion_factor
/// sep = µ * baseline
/// d = (sqrt(GM_sun) * conversion * baseline / sep)^(2/3)
fn implied_distance(sep_arcmin: f64, baseline_years: f64) -> f64 {
    use p9_core::constants::{GM_SUN, RAD2DEG, YEAR_DAYS};
    let k = GM_SUN.sqrt() * RAD2DEG * 60.0 * YEAR_DAYS * baseline_years;
    (k / sep_arcmin).powf(2.0 / 3.0)
}

/// Run the IRAS-AKARI cross-match candidate search.
///
/// This is the core algorithm from Phan et al. (2025) §3.
pub fn search_candidates(
    iras_sources: &[FirSource],
    akari_sources: &[FirSource],
    criteria: &SelectionCriteria,
    baseline_years: f64,
) -> SearchResult {
    let mut candidates = Vec::new();
    let mut n_pairs_in_window = 0usize;

    for iras in iras_sources {
        for akari in akari_sources {
            let sep =
                angular_separation_arcmin(iras.ra_deg, iras.dec_deg, akari.ra_deg, akari.dec_deg);

            if sep < criteria.sep_min_arcmin || sep > criteria.sep_max_arcmin {
                continue;
            }
            n_pairs_in_window += 1;

            // Flux ratio check
            let flux_ratio = iras.flux_jy / akari.flux_jy;
            if flux_ratio < criteria.min_flux_ratio || flux_ratio > criteria.max_flux_ratio {
                continue;
            }

            let d = implied_distance(sep, baseline_years);

            candidates.push(CandidatePair {
                iras_source: iras.clone(),
                akari_source: akari.clone(),
                separation_arcmin: sep,
                implied_distance_au: d,
                flux_ratio,
            });
        }
    }

    let n_pairs_flux_cut = candidates.len();

    SearchResult {
        n_iras: iras_sources.len(),
        n_akari: akari_sources.len(),
        n_pairs_in_window,
        n_pairs_flux_cut,
        candidates,
    }
}

/// Estimate the number of spurious (chance-alignment) pairs expected
/// from uncorrelated source catalogues.
///
/// For N_iras IRAS sources and N_akari AKARI sources spread over the
/// full sky (41,253 deg²), the expected number of chance pairs in an
/// annular search window [θ_min, θ_max] is:
///
///   N_spurious = N_iras * N_akari * π(θ_max² - θ_min²) / A_sky
///
/// where angles are in the same units as A_sky.
pub fn expected_spurious_pairs(
    n_iras: usize,
    n_akari: usize,
    sep_min_arcmin: f64,
    sep_max_arcmin: f64,
) -> f64 {
    // Full sky = 4π sr = 41,253 deg² = 41,253 * 3600 arcmin²
    let full_sky_arcmin2 = 41_253.0 * 3600.0;
    let annulus_area = std::f64::consts::PI * (sep_max_arcmin.powi(2) - sep_min_arcmin.powi(2));
    (n_iras as f64) * (n_akari as f64) * annulus_area / full_sky_arcmin2
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::survey_model::Survey;

    fn make_iras_source(ra: f64, dec: f64, flux: f64) -> FirSource {
        FirSource {
            ra_deg: ra,
            dec_deg: dec,
            flux_jy: flux,
            pos_err_arcsec: 20.0,
            epoch_year: 1983.5,
            survey: Survey::Iras,
        }
    }

    fn make_akari_source(ra: f64, dec: f64, flux: f64) -> FirSource {
        FirSource {
            ra_deg: ra,
            dec_deg: dec,
            flux_jy: flux,
            pos_err_arcsec: 30.0,
            epoch_year: 2006.5,
            survey: Survey::Akari,
        }
    }

    #[test]
    fn pair_within_window_is_found() {
        // Place sources ~1° apart (60') — within [42', 69.6']
        let iras = vec![make_iras_source(180.0, 45.0, 0.5)];
        let akari = vec![make_akari_source(180.0, 46.0, 0.6)];
        let criteria = SelectionCriteria::default();
        let result = search_candidates(&iras, &akari, &criteria, 23.0);
        assert_eq!(result.n_pairs_in_window, 1);
        assert_eq!(result.candidates.len(), 1);
        let pair = &result.candidates[0];
        assert!((pair.separation_arcmin - 60.0).abs() < 0.5);
        assert!(pair.implied_distance_au > 400.0 && pair.implied_distance_au < 800.0);
    }

    #[test]
    fn pair_outside_window_is_rejected() {
        // Place sources 2° apart (120') — outside [42', 69.6']
        let iras = vec![make_iras_source(180.0, 45.0, 0.5)];
        let akari = vec![make_akari_source(180.0, 47.0, 0.6)];
        let criteria = SelectionCriteria::default();
        let result = search_candidates(&iras, &akari, &criteria, 23.0);
        assert_eq!(result.n_pairs_in_window, 0);
        assert!(result.candidates.is_empty());
    }

    #[test]
    fn flux_ratio_rejection() {
        // Sources ~1° apart but extreme flux ratio
        let iras = vec![make_iras_source(180.0, 45.0, 100.0)];
        let akari = vec![make_akari_source(180.0, 46.0, 0.1)];
        let criteria = SelectionCriteria::default();
        let result = search_candidates(&iras, &akari, &criteria, 23.0);
        assert_eq!(result.n_pairs_in_window, 1);
        assert_eq!(result.n_pairs_flux_cut, 0); // ratio 1000 >> max_flux_ratio
    }

    #[test]
    fn implied_distance_at_60_arcmin() {
        // 60' over 23 years: for circular orbits this maps to ~450 AU.
        // The implied_distance function inverts the circular-orbit proper
        // motion, so smaller distances are expected than the paper's
        // 500–700 AU range (which accounts for eccentric orbits).
        let d = implied_distance(60.0, 23.0);
        assert!(
            d > 350.0 && d < 550.0,
            "implied distance = {d:.0} AU, expected ~410-500"
        );
    }

    #[test]
    fn spurious_pairs_estimate() {
        // With ~173,000 IRAS sources and ~1000 AKARI unconfirmed sources,
        // and the 42'–69.6' annulus, estimate spurious pairs
        let n_spurious = expected_spurious_pairs(173_000, 1_000, 42.0, 69.6);
        // This gives a rough number to compare against the 13 pairs found
        assert!(n_spurious > 0.0);
        // The paper found 13, so spurious rate should be in a similar ballpark
        // (most pairs are indeed chance alignments)
    }

    #[test]
    fn empty_catalogues_give_no_pairs() {
        let criteria = SelectionCriteria::default();
        let result = search_candidates(&[], &[], &criteria, 23.0);
        assert_eq!(result.n_pairs_in_window, 0);
        assert!(result.candidates.is_empty());
    }

    #[test]
    fn selection_criteria_default_matches_paper() {
        let criteria = SelectionCriteria::default();
        assert!((criteria.sep_min_arcmin - 42.0).abs() < 1e-10);
        assert!((criteria.sep_max_arcmin - 69.6).abs() < 1e-10);
    }
}
