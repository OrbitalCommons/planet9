//! Candidate pair search and selection criteria.
//!
//! Reproduces the source-pair matching algorithm from Phan et al. (2025):
//!
//! 1. Start with IRAS FSC 60 µm sources and AKARI Monthly Unconfirmed 90 µm sources
//! 2. Apply flux cuts to retain only sources consistent with P9 thermal emission
//! 3. For each IRAS-AKARI pair, compute angular separation
//! 4. Retain pairs with separation in [42', 69.6'] (corresponding to 500–700 AU)
//! 5. Apply the paper's colour-consistency check between the two surveys
//! 6. Image inspection (not reproducible computationally — we flag candidates)
//!
//! The paper finds 22 raw pairs in the separation window and 13 after the
//! colour-consistency cut (its Fig. 4), with 1 good candidate after image
//! inspection. Two pipelines are provided: [`search_candidates_paper`] — the
//! paper's published cuts verbatim (within-survey colour rejection,
//! cross-survey inconsistency, colour consistency, upper flux limit,
//! galactic plane/bulge exclusion, un-inflated window) — and the labelled
//! thermal-window variant [`search_candidates`] retained as an alternative
//! (blackbody-derived flux-ratio window, 3σ-inflated separation).

use p9_core::coords::candidate_pair::implied_distance;
use p9_core::units::{arcseconds, au, Angle, Length};
use serde::{Deserialize, Serialize};

use crate::survey_model::{angular_separation_arcmin, FirSource, Survey};
use crate::thermal_model::flux_ratio_60_90;

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
    /// Separation window inflation in units of the pair's combined 1σ
    /// positional uncertainty (the catalogued 20"/30" errors were
    /// previously carried but never used in the match).
    pub pos_err_n_sigma: f64,
}

impl SelectionCriteria {
    /// Minimum angular separation as a dimension-checked [`Angle`].
    pub fn sep_min(&self) -> Angle {
        arcseconds(self.sep_min_arcmin * 60.0)
    }

    /// Maximum angular separation as a dimension-checked [`Angle`].
    pub fn sep_max(&self) -> Angle {
        arcseconds(self.sep_max_arcmin * 60.0)
    }
}

impl SelectionCriteria {
    /// Criteria with the flux-ratio window derived from the blackbody
    /// thermal model over `[t_min_k, t_max_k]`, widened by a fractional
    /// `margin` for flux-calibration error and emissivity deviations.
    ///
    /// For 30–50 K and margin 0.5 this gives ≈ [0.16, 0.99] — physically
    /// motivated, unlike the previous ad-hoc [0.05, 3.0].
    pub fn from_thermal_model(t_min_k: f64, t_max_k: f64, margin: f64) -> Self {
        let r_lo = flux_ratio_60_90(t_min_k);
        let r_hi = flux_ratio_60_90(t_max_k);
        let (r_lo, r_hi) = (r_lo.min(r_hi), r_lo.max(r_hi));
        Self {
            sep_min_arcmin: 42.0,
            sep_max_arcmin: 69.6,
            min_flux_ratio: r_lo / (1.0 + margin),
            max_flux_ratio: r_hi * (1.0 + margin),
            pos_err_n_sigma: 3.0,
        }
    }
}

impl Default for SelectionCriteria {
    /// Paper window 42'–69.6'; flux-ratio bounds derived from the 30–50 K
    /// thermal model with a 50% margin; 3σ positional-error inflation.
    fn default() -> Self {
        Self::from_thermal_model(30.0, 50.0, 0.5)
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
    /// Combined 1σ positional uncertainty of the pair (arcmin)
    pub combined_pos_err_arcmin: f64,
    /// Implied heliocentric distance in AU (geocentric maximum-distance
    /// inversion, see `proper_motion::implied_distance`)
    pub implied_distance_au: f64,
    /// Flux ratio (IRAS 60µm / AKARI 90µm)
    pub flux_ratio: f64,
}

impl CandidatePair {
    /// Angular separation as a dimension-checked [`Angle`].
    pub fn separation(&self) -> Angle {
        arcseconds(self.separation_arcmin * 60.0)
    }

    /// Combined 1σ positional uncertainty as a dimension-checked [`Angle`].
    pub fn combined_pos_err(&self) -> Angle {
        arcseconds(self.combined_pos_err_arcmin * 60.0)
    }

    /// Implied heliocentric distance as a dimension-checked [`Length`].
    pub fn implied_distance(&self) -> Length {
        au(self.implied_distance_au)
    }
}

/// Result of the candidate search.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SearchResult {
    /// Number of IRAS sources considered
    pub n_iras: usize,
    /// Number of AKARI sources considered
    pub n_akari: usize,
    /// Number of pairs within the (uncertainty-inflated) separation window
    pub n_pairs_in_window: usize,
    /// Number of pairs passing the flux ratio cut
    pub n_pairs_flux_cut: usize,
    /// Final candidate pairs
    pub candidates: Vec<CandidatePair>,
}

/// Combined 1σ positional uncertainty of a pair in arcminutes.
fn combined_pos_err_arcmin(iras: &FirSource, akari: &FirSource) -> f64 {
    (iras.pos_err_arcsec.powi(2) + akari.pos_err_arcsec.powi(2)).sqrt() / 60.0
}

/// Run the IRAS-AKARI cross-match candidate search.
///
/// This is the core algorithm from Phan et al. (2025) §3. The separation
/// window is inflated per-pair by `pos_err_n_sigma` × the combined
/// positional uncertainty, so genuine pairs scattered just outside the
/// nominal window by measurement error are not lost.
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

            let sigma = combined_pos_err_arcmin(iras, akari);
            let pad = criteria.pos_err_n_sigma * sigma;
            if sep < criteria.sep_min_arcmin - pad || sep > criteria.sep_max_arcmin + pad {
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
                combined_pos_err_arcmin: sigma,
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

/// The paper's published per-source and per-pair cuts (Phan et al. 2025
/// §3.3 / Fig. 4), verbatim:
///
/// per source —
///  - flux-quality rejection (quality 1 in any of 60/65/90/100 µm)
///  - galactic plane/bulge exclusion: reject |b| < 10° or within 27.5° of
///    the galactic centre
///  - within-survey colour rejection (T ≳ 100 K): reject if
///    S₁₀₀/S₆₀ > 3 or S₆₀/S₁₀₀ > 3 (IRAS), S₉₀/S₆₅ > 3 or S₆₅/S₉₀ > 3
///    (AKARI)
///  - upper flux limit: reject fluxes above the 500 AU expectation
///
/// per pair —
///  - separation ∈ [42′, 69.6′], un-inflated
///  - cross-survey inconsistency: reject if S₆₀/S₆₅ or S₆₅/S₆₀ > 4, or
///    S₉₀/S₁₀₀ or S₁₀₀/S₉₀ > 4
///  - colour consistency: reject when the IRAS colour (S₁₀₀/S₆₀) and AKARI
///    colour (S₉₀/S₆₅) differ by more than a factor 3.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PaperCuts {
    /// Within-survey colour rejection threshold (3).
    pub within_survey_ratio_max: f64,
    /// Cross-survey flux inconsistency threshold (4).
    pub cross_survey_ratio_max: f64,
    /// Colour-consistency threshold between the surveys (3).
    pub colour_ratio_max: f64,
    /// Upper flux limits (Jy) at the 500 AU expectation: (IRAS 60 µm,
    /// AKARI 90 µm).
    pub upper_flux_limit_jy: (f64, f64),
    /// Galactic plane half-width (deg).
    pub plane_b_min_deg: f64,
    /// Bulge exclusion radius about the galactic centre (deg).
    pub bulge_radius_deg: f64,
    /// Separation window (arcmin), un-inflated.
    pub sep_window_arcmin: (f64, f64),
}

impl Default for PaperCuts {
    fn default() -> Self {
        // Upper limits: the warm/heavy corner of the paper's grid at 500 AU
        // (10 M⊕, 50 K) — fluxes above this cannot be the sought planet.
        let warm = crate::thermal_model::P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 500.0,
            t_eff: 50.0,
        };
        Self {
            within_survey_ratio_max: 3.0,
            cross_survey_ratio_max: 4.0,
            colour_ratio_max: 3.0,
            upper_flux_limit_jy: (warm.flux_density_jy(60.0e-6), warm.flux_density_jy(90.0e-6)),
            plane_b_min_deg: 10.0,
            bulge_radius_deg: 27.5,
            sep_window_arcmin: (42.0, 69.6),
        }
    }
}

/// Galactic latitude (deg) and angular distance to the galactic centre (deg)
/// for an equatorial direction, via the shared p9-core frame matrices.
fn galactic_b_and_gc_distance(ra_deg: f64, dec_deg: f64) -> (f64, f64) {
    use p9_core::coords::sky::equatorial_to_galactic_matrix;
    let (ra, dec) = (ra_deg.to_radians(), dec_deg.to_radians());
    let eq = nalgebra::Vector3::new(dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin());
    let gal = equatorial_to_galactic_matrix() * eq;
    let b = (gal.z / gal.norm()).asin().to_degrees();
    let gc_cos = gal.x / gal.norm(); // galactic centre is +x (l = 0, b = 0)
    (b, gc_cos.clamp(-1.0, 1.0).acos().to_degrees())
}

/// Does a source survive the paper's per-source cuts?
pub fn passes_source_cuts(src: &FirSource, cuts: &PaperCuts) -> bool {
    if !src.flux_quality_ok {
        return false;
    }
    let (b, gc_dist) = galactic_b_and_gc_distance(src.ra_deg, src.dec_deg);
    if b.abs() < cuts.plane_b_min_deg || gc_dist < cuts.bulge_radius_deg {
        return false;
    }
    // Within-survey colour rejection (skipped when no secondary flux is
    // catalogued, as in the sparse MUSL rows).
    if let Some(sec) = src.flux_secondary_jy {
        let r = src.flux_jy / sec;
        if r > cuts.within_survey_ratio_max || 1.0 / r > cuts.within_survey_ratio_max {
            return false;
        }
    }
    let limit = match src.survey {
        Survey::Iras => cuts.upper_flux_limit_jy.0,
        Survey::Akari => cuts.upper_flux_limit_jy.1,
    };
    src.flux_jy <= limit
}

/// Does a (source-cut-surviving) pair pass the paper's per-pair cuts?
/// Returns the separation (arcmin) when it does.
pub fn passes_pair_cuts(iras: &FirSource, akari: &FirSource, cuts: &PaperCuts) -> Option<f64> {
    let sep = angular_separation_arcmin(iras.ra_deg, iras.dec_deg, akari.ra_deg, akari.dec_deg);
    let (lo, hi) = cuts.sep_window_arcmin;
    if sep < lo || sep > hi {
        return None;
    }
    // Cross-survey flux inconsistency (needs the secondary bands).
    if let (Some(s100), Some(s65)) = (iras.flux_secondary_jy, akari.flux_secondary_jy) {
        let (s60, s90) = (iras.flux_jy, akari.flux_jy);
        for r in [s60 / s65, s65 / s60, s90 / s100, s100 / s90] {
            if r > cuts.cross_survey_ratio_max {
                return None;
            }
        }
        // Colour consistency between the surveys.
        let colour_iras = s100 / s60;
        let colour_akari = s90 / s65;
        let c = colour_iras / colour_akari;
        if c > cuts.colour_ratio_max || 1.0 / c > cuts.colour_ratio_max {
            return None;
        }
    }
    Some(sep)
}

/// The paper's published pipeline end to end: per-source cuts, un-inflated
/// separation window, cross-survey and colour-consistency cuts.
pub fn search_candidates_paper(
    iras_sources: &[FirSource],
    akari_sources: &[FirSource],
    cuts: &PaperCuts,
    baseline_years: f64,
) -> SearchResult {
    let iras_kept: Vec<&FirSource> = iras_sources
        .iter()
        .filter(|s| passes_source_cuts(s, cuts))
        .collect();
    let akari_kept: Vec<&FirSource> = akari_sources
        .iter()
        .filter(|s| passes_source_cuts(s, cuts))
        .collect();

    let mut candidates = Vec::new();
    let mut n_pairs_in_window = 0usize;
    let (lo, hi) = cuts.sep_window_arcmin;
    for iras in &iras_kept {
        for akari in &akari_kept {
            let sep =
                angular_separation_arcmin(iras.ra_deg, iras.dec_deg, akari.ra_deg, akari.dec_deg);
            if sep >= lo && sep <= hi {
                n_pairs_in_window += 1;
            }
            if let Some(sep) = passes_pair_cuts(iras, akari, cuts) {
                candidates.push(CandidatePair {
                    iras_source: (*iras).clone(),
                    akari_source: (*akari).clone(),
                    separation_arcmin: sep,
                    combined_pos_err_arcmin: combined_pos_err_arcmin(iras, akari),
                    implied_distance_au: implied_distance(sep, baseline_years),
                    flux_ratio: iras.flux_jy / akari.flux_jy,
                });
            }
        }
    }

    let n_pairs_flux_cut = candidates.len();
    SearchResult {
        n_iras: iras_kept.len(),
        n_akari: akari_kept.len(),
        n_pairs_in_window,
        n_pairs_flux_cut,
        candidates,
    }
}

/// Galactic-latitude-dependent source-density model for far-infrared
/// catalogues.
///
/// IRAS 60 µm (and AKARI 90 µm) source counts are strongly concentrated
/// toward the galactic plane (young stellar objects, HII regions, cirrus
/// knots) on top of a roughly isotropic extragalactic floor. We model the
/// relative density as
///
///   w(b) = f_disk · exp(−|sin b| / s) / Z + (1 − f_disk),
///
/// with Z = s(1 − e^{−1/s}) normalizing the disk term to unit sky average,
/// so ⟨w⟩ = 1 and N_total is preserved. Defaults: f_disk = 0.5 of sources
/// in the disk component, scale s = 0.15 (|b| scale height ≈ 9°) — a
/// simple documented stand-in for the FSC count maps.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct GalacticDensityModel {
    /// Fraction of sources in the plane-concentrated component
    pub disk_fraction: f64,
    /// Exponential scale in sin|b|
    pub sin_b_scale: f64,
}

impl Default for GalacticDensityModel {
    fn default() -> Self {
        Self {
            disk_fraction: 0.5,
            sin_b_scale: 0.15,
        }
    }
}

impl GalacticDensityModel {
    /// Relative density w(b) at galactic latitude `b_deg` (sky average 1).
    pub fn relative_density(&self, b_deg: f64) -> f64 {
        let s = self.sin_b_scale;
        let z = s * (1.0 - (-1.0 / s).exp());
        let u = b_deg.to_radians().sin().abs();
        self.disk_fraction * (-u / s).exp() / z + (1.0 - self.disk_fraction)
    }

    /// Chance-pair enhancement factor ⟨w²⟩ ≥ 1: both catalogues share the
    /// latitude concentration, so chance pairs are more likely than the
    /// uniform-sky estimate by the sky average of w_iras·w_akari
    /// (= ⟨w²⟩ for a shared model). Computed by numerical integration over
    /// sin b.
    pub fn pair_enhancement(&self) -> f64 {
        let n = 2000;
        let mut sum = 0.0;
        for k in 0..n {
            // uniform in sin b ∈ [0, 1] — the area-weighted measure
            let u = (k as f64 + 0.5) / n as f64;
            let b = u.asin().to_degrees();
            sum += self.relative_density(b).powi(2);
        }
        sum / n as f64
    }
}

/// Fraction of the sky RETAINED after the paper's galactic cuts (|b| ≥ 10°
/// plus the 27.5°-radius bulge exclusion, Phan et al. Fig. 4): ≈ 0.79.
pub const RETAINED_SKY_FRACTION: f64 = 0.79;

/// Paper post-cut catalogue counts (Phan et al. Fig. 4): IRAS sources and
/// AKARI-MUSL sources surviving all flux/colour/plane cuts, i.e. the
/// populations the pairing actually runs on.
pub const PAPER_POST_CUT_IRAS: usize = 2479;
pub const PAPER_POST_CUT_AKARI: usize = 108;

/// Expected number of chance-alignment pairs in the separation annulus,
/// conditioned on the *post-cut* source counts over the *retained* sky:
///
///   λ = N_iras · N_akari · π(θ_max² − θ_min²) / A_retained.
///
/// The paper removes |b| < 10° and the bulge BEFORE pairing, so the post-cut
/// populations live only on the retained ~79% of the sky where the disk
/// density term is absent — the correct area factor is 1/A_retained with NO
/// galactic pair enhancement. (A previous version divided by the full sky
/// and multiplied by ⟨w²⟩ ≈ 1.59 — modeling plane-concentrated catalogues
/// that the paper's cuts had already excluded — which inflated λ to ~19 and
/// was most of the documented λ-vs-13 residual. With the Fig. 4 counts this
/// form gives λ ≈ 22, matching the paper's 22 raw window pairs.)
pub fn expected_chance_pairs(
    n_iras_post_cut: usize,
    n_akari_post_cut: usize,
    sep_min_arcmin: f64,
    sep_max_arcmin: f64,
    _density: &GalacticDensityModel,
) -> f64 {
    // Full sky = 41,253 deg² = 41,253 × 3600 arcmin²
    let retained_arcmin2 = 41_253.0 * 3600.0 * RETAINED_SKY_FRACTION;
    let annulus_area = std::f64::consts::PI * (sep_max_arcmin.powi(2) - sep_min_arcmin.powi(2));
    (n_iras_post_cut as f64) * (n_akari_post_cut as f64) * annulus_area / retained_arcmin2
}

/// Poisson false-alarm probability of finding at least one surviving pair
/// by chance: P(N ≥ 1) = 1 − e^{−λ}. This is the key statistic for the
/// single candidate the paper retains after image inspection.
pub fn false_alarm_probability(expected_pairs: f64) -> f64 {
    1.0 - (-expected_pairs).exp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::survey_model::Survey;
    use approx::assert_relative_eq;
    use rand::{Rng, SeedableRng};
    use uom::si::angle::minute;
    use uom::si::length::astronomical_unit;

    #[test]
    fn typed_accessors_match_f64_sources() {
        let criteria = SelectionCriteria::default();
        assert_relative_eq!(
            criteria.sep_min().get::<minute>(),
            criteria.sep_min_arcmin,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            criteria.sep_max().get::<minute>(),
            criteria.sep_max_arcmin,
            max_relative = 1e-12
        );

        let iras = vec![make_iras_source(180.0, 45.0, 0.5)];
        let akari = vec![make_akari_source(180.0, 46.0, 0.6)];
        let pair = &search_candidates(&iras, &akari, &criteria, 23.0).candidates[0];
        assert_relative_eq!(
            pair.separation().get::<minute>(),
            pair.separation_arcmin,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            pair.combined_pos_err().get::<minute>(),
            pair.combined_pos_err_arcmin,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            pair.implied_distance().get::<astronomical_unit>(),
            pair.implied_distance_au,
            epsilon = 1e-9
        );
    }

    fn make_iras_source(ra: f64, dec: f64, flux: f64) -> FirSource {
        // Cold-consistent secondary flux (S100 ≈ 2.5·S60 for ~40 K).
        FirSource {
            flux_secondary_jy: Some(flux * 2.5),
            flux_quality_ok: true,
            ra_deg: ra,
            dec_deg: dec,
            flux_jy: flux,
            pos_err_arcsec: 20.0,
            epoch_year: 1983.5,
            survey: Survey::Iras,
        }
    }

    fn make_akari_source(ra: f64, dec: f64, flux: f64) -> FirSource {
        // Cold-consistent secondary flux (S65 ≈ 0.4·S90 for ~40 K).
        FirSource {
            flux_secondary_jy: Some(flux * 0.4),
            flux_quality_ok: true,
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
        // Place sources ~1° apart (60') — within [42', 69.6'].
        // Flux ratio 0.5/0.6 ≈ 0.83 sits inside the thermal window.
        let iras = vec![make_iras_source(180.0, 45.0, 0.5)];
        let akari = vec![make_akari_source(180.0, 46.0, 0.6)];
        let criteria = SelectionCriteria::default();
        let result = search_candidates(&iras, &akari, &criteria, 23.0);
        assert_eq!(result.n_pairs_in_window, 1);
        assert_eq!(result.candidates.len(), 1);
        let pair = &result.candidates[0];
        assert!((pair.separation_arcmin - 60.0).abs() < 0.5);
        // Geocentric maximum-distance inversion: 60' → ~560 AU.
        assert!(
            pair.implied_distance_au > 500.0 && pair.implied_distance_au < 630.0,
            "implied d = {:.0} AU",
            pair.implied_distance_au
        );
    }

    #[test]
    fn pair_outside_window_is_rejected() {
        // Place sources 2° apart (120') — outside even the inflated window
        let iras = vec![make_iras_source(180.0, 45.0, 0.5)];
        let akari = vec![make_akari_source(180.0, 47.0, 0.6)];
        let criteria = SelectionCriteria::default();
        let result = search_candidates(&iras, &akari, &criteria, 23.0);
        assert_eq!(result.n_pairs_in_window, 0);
        assert!(result.candidates.is_empty());
    }

    #[test]
    fn positional_uncertainty_inflates_window() {
        // Combined σ = √(20² + 30²)/60 = 0.6'; 3σ = 1.8'. A pair at
        // 70.5' (0.9' beyond the nominal 69.6') must be kept, one at
        // 75' rejected.
        let criteria = SelectionCriteria::default();
        let iras = vec![make_iras_source(180.0, 45.0, 0.5)];
        let akari_in = vec![make_akari_source(180.0, 45.0 + 70.5 / 60.0, 0.6)];
        let akari_out = vec![make_akari_source(180.0, 45.0 + 75.0 / 60.0, 0.6)];
        assert_eq!(
            search_candidates(&iras, &akari_in, &criteria, 23.0).n_pairs_in_window,
            1
        );
        assert_eq!(
            search_candidates(&iras, &akari_out, &criteria, 23.0).n_pairs_in_window,
            0
        );
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
    fn flux_ratio_window_derived_from_thermal_model() {
        let criteria = SelectionCriteria::default();
        // 30–50 K blackbody ratio range 0.23–0.66, ±50% margin
        assert!(
            (criteria.min_flux_ratio - 0.233 / 1.5).abs() < 0.02,
            "min ratio = {:.3}",
            criteria.min_flux_ratio
        );
        assert!(
            (criteria.max_flux_ratio - 0.66 * 1.5).abs() < 0.05,
            "max ratio = {:.3}",
            criteria.max_flux_ratio
        );
        // The paper's good candidate (0.24 Jy / 0.27 Jy = 0.89) passes;
        // the old ad-hoc window [0.05, 3.0] does not bound it usefully.
        let candidate_ratio = 0.24 / 0.27;
        assert!(
            candidate_ratio > criteria.min_flux_ratio && candidate_ratio < criteria.max_flux_ratio
        );
    }

    #[test]
    fn chance_pairs_conditioned_on_post_cut_counts() {
        let density = GalacticDensityModel::default();
        // Post-flux-cut populations of order 10³ (IRAS) × 10² (AKARI)
        // give an expectation of the same order as the paper's 13 pairs.
        let lambda = expected_chance_pairs(1_800, 100, 42.0, 69.6, &density);
        assert!(
            lambda > 5.0 && lambda < 50.0,
            "λ = {lambda:.1}, expected same order as the paper's 13 pairs"
        );

        // The uniform full-catalogue estimate is ~3 orders of magnitude
        // larger — the error this conditioning fixes.
        let lambda_full = expected_chance_pairs(173_000, 1_000, 42.0, 69.6, &density);
        assert!(
            lambda_full / lambda > 100.0,
            "full-catalogue λ = {lambda_full:.0} should dwarf post-cut λ = {lambda:.1}"
        );
    }

    #[test]
    fn galactic_density_model_normalized() {
        let density = GalacticDensityModel::default();
        // Sky average of w(b) is 1 (area-weighted, uniform in sin b).
        let n = 4000;
        let mean: f64 = (0..n)
            .map(|k| {
                let u = (k as f64 + 0.5) / n as f64;
                density.relative_density(u.asin().to_degrees())
            })
            .sum::<f64>()
            / n as f64;
        assert!((mean - 1.0).abs() < 0.01, "⟨w⟩ = {mean:.3}");
        // Plane denser than pole; enhancement ⟨w²⟩ > 1.
        assert!(density.relative_density(0.0) > density.relative_density(90.0));
        let enh = density.pair_enhancement();
        assert!(enh > 1.0 && enh < 5.0, "⟨w²⟩ = {enh:.2}");
    }

    #[test]
    fn false_alarm_probability_limits() {
        assert!((false_alarm_probability(0.01) - 0.00995).abs() < 1e-4);
        assert!(false_alarm_probability(20.0) > 0.999);
        assert!(false_alarm_probability(0.0).abs() < 1e-12);
    }

    #[test]
    fn empty_catalogues_give_no_pairs() {
        let criteria = SelectionCriteria::default();
        let result = search_candidates(&[], &[], &criteria, 23.0);
        assert_eq!(result.n_pairs_in_window, 0);
        assert!(result.candidates.is_empty());
    }

    #[test]
    fn selection_criteria_default_matches_paper_window() {
        let criteria = SelectionCriteria::default();
        assert!((criteria.sep_min_arcmin - 42.0).abs() < 1e-10);
        assert!((criteria.sep_max_arcmin - 69.6).abs() < 1e-10);
    }

    /// Seeded end-to-end smoke test: synthetic catalogues with one
    /// injected moving source; the pipeline must recover the injected
    /// pair and the chance-pair count must be of the expected order.
    #[test]
    fn end_to_end_pipeline_recovers_injected_candidate() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(20250610);

        // Background populations in a 30°×30° patch.
        let patch = |rng: &mut rand::rngs::StdRng| {
            (
                150.0 + rng.gen::<f64>() * 30.0,
                30.0 + rng.gen::<f64>() * 30.0,
            )
        };
        let mut iras: Vec<FirSource> = (0..200)
            .map(|_| {
                let (ra, dec) = patch(&mut rng);
                make_iras_source(ra, dec, 0.2 + rng.gen::<f64>() * 0.4)
            })
            .collect();
        let mut akari: Vec<FirSource> = (0..50)
            .map(|_| {
                let (ra, dec) = patch(&mut rng);
                make_akari_source(ra, dec, 0.55 + rng.gen::<f64>() * 0.4)
            })
            .collect();

        // Injected Planet Nine: moved 50' south between epochs, with a
        // thermally consistent flux ratio (0.30/0.50 = 0.6).
        iras.push(make_iras_source(165.0, 45.0, 0.30));
        akari.push(make_akari_source(165.0, 45.0 - 50.0 / 60.0, 0.50));

        let criteria = SelectionCriteria::default();
        let result = search_candidates(&iras, &akari, &criteria, 23.0);

        let injected = result.candidates.iter().find(|c| {
            (c.iras_source.ra_deg - 165.0).abs() < 1e-9 && (c.separation_arcmin - 50.0).abs() < 0.5
        });
        let injected = injected.expect("injected pair must be recovered");
        assert!(
            injected.implied_distance_au > 450.0 && injected.implied_distance_au < 800.0,
            "implied d = {:.0} AU",
            injected.implied_distance_au
        );

        // The chance-pair count in the patch should be of the same order
        // as the Poisson expectation scaled to the patch area (the patch
        // is 900 deg² of the 41,253 deg² sky; densities are patch-local
        // so we compare directly against the patch-scaled estimate).
        let patch_area_arcmin2 = 30.0 * 30.0 * 3600.0;
        let annulus = std::f64::consts::PI
            * (criteria.sep_max_arcmin.powi(2) - criteria.sep_min_arcmin.powi(2));
        let lambda_patch = 200.0 * 50.0 * annulus / patch_area_arcmin2;
        let n_chance = result.n_pairs_in_window as f64 - 1.0; // minus injected
        assert!(
            n_chance > 0.2 * lambda_patch && n_chance < 5.0 * lambda_patch,
            "chance pairs {n_chance} vs expectation {lambda_patch:.1}"
        );
    }

    #[test]
    fn lambda_with_paper_counts_matches_raw_window_pairs() {
        // With the paper's Fig. 4 post-cut counts over the retained sky,
        // the expected chance pairs in the 42'-69.6' window is ~22 —
        // matching the paper's 22 raw window pairs (13 remain after their
        // colour-consistency cut).
        let density = GalacticDensityModel::default();
        let lambda = expected_chance_pairs(
            PAPER_POST_CUT_IRAS,
            PAPER_POST_CUT_AKARI,
            42.0,
            69.6,
            &density,
        );
        assert!((18.0..26.0).contains(&lambda), "lambda = {lambda:.1}");
    }

    #[test]
    fn paper_cuts_exercise_each_criterion() {
        let cuts = PaperCuts::default();
        // A cold, off-plane, quality source passes. RA/Dec ≈ (30°, -40°) is
        // far from the plane and bulge.
        let good = make_iras_source(30.0, -40.0, 0.3);
        assert!(passes_source_cuts(&good, &cuts));

        // Hot within-survey colour (S100/S60 < 1/3) rejected.
        let mut hot = good.clone();
        hot.flux_secondary_jy = Some(good.flux_jy / 4.0);
        assert!(!passes_source_cuts(&hot, &cuts));

        // Bad flux quality rejected.
        let mut badq = good.clone();
        badq.flux_quality_ok = false;
        assert!(!passes_source_cuts(&badq, &cuts));

        // Galactic-plane source rejected (galactic centre direction:
        // RA ≈ 266.4°, Dec ≈ -29°).
        let plane = make_iras_source(266.4, -29.0, 0.3);
        assert!(!passes_source_cuts(&plane, &cuts));

        // Over-bright source (above the 500 AU expectation) rejected.
        let mut bright = good.clone();
        bright.flux_jy = cuts.upper_flux_limit_jy.0 * 2.0;
        bright.flux_secondary_jy = Some(bright.flux_jy * 2.5);
        assert!(!passes_source_cuts(&bright, &cuts));

        // Pair: consistent cold pair at 50' passes; colour-inconsistent
        // partner fails the factor-3 colour cut.
        let akari = make_akari_source(30.0, -40.0 + 50.0 / 60.0, 0.25);
        assert!(passes_pair_cuts(&good, &akari, &cuts).is_some());
        let mut off_colour = akari.clone();
        // AKARI colour S90/S65 = 1/0.4 = 2.5 nominally; push it past 3x the
        // IRAS colour (2.5) by shrinking S65.
        off_colour.flux_secondary_jy = Some(akari.flux_jy / 20.0);
        assert!(passes_pair_cuts(&good, &off_colour, &cuts).is_none());

        // End to end: the good pair is found by the paper pipeline.
        let res =
            search_candidates_paper(&[good.clone(), hot, plane], &[akari.clone()], &cuts, 23.5);
        assert_eq!(res.n_iras, 1, "only the good IRAS source survives");
        assert_eq!(res.candidates.len(), 1);
    }
}
