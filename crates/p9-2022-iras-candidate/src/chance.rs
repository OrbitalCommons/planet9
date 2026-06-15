//! Chance-alignment probability for the single-HCON association.
//!
//! Rowan-Robinson's signature is a pair (or triplet) of single-HCON Reject
//! File sources with separations 2–35 arcmin — the body's motion between
//! the IRAS HCON passes. The natural question is how often such a pairing
//! arises *by chance* given the IRAS source surface density. For a Poisson
//! field of surface density σ, the expected number of neighbours of a given
//! source within an annulus `[r_min, r_max]` (arcmin) is
//!
//!   μ = σ · π · (r_max² − r_min²),
//!
//! and the probability of *at least one* chance neighbour is `1 − e^{−μ}`.
//! Multiplying μ by the number of search sources gives the expected number
//! of chance associations across the whole search — the paper's "several
//! hundred candidate associations" — which is precisely why a *single*
//! survivor of the stricter HCON detect/non-detect pattern is presented as
//! tentative rather than as a detection.
//!
//! # The search is over subsets, not the whole FSC
//!
//! The "several hundred" figure is *not* every FSC source tested against
//! every other FSC source — that all-pairs estimate runs to ~10⁵–10⁶ (an
//! upper bound, pinned in [`tests::all_pairs_upper_bound_is_large`]).
//! Rowan-Robinson searched the much smaller pool of *unidentified* 60 µm
//! sources ([`N_UNIDENTIFIED_60UM`]) for a neighbour drawn from the sparse
//! single-HCON Reject File ([`N_REJECT_SINGLE_HCON`]). With those two
//! populations the raw Poisson estimate is ~10³ geometric coincidences (the
//! paper's "several hundred" after Scanpi vetting thinned them) rather than
//! the ~10⁵–10⁶ all-pairs figure. The exact subset sizes are not tabulated
//! in the paper, so the
//! constants below are order-of-magnitude reference values labelled as such;
//! the reproduction therefore pins the *order of magnitude* ("several
//! hundred"), not a precise count.

/// IRAS Faint Source Catalogue 60 µm source count (≈ 173 000 sources).
pub const IRAS_FSC_N_60UM: f64 = 173_000.0;

/// Approximate number of *unidentified* 60 µm sources searched (the pool of
/// candidate primaries, after removing 2MASS galaxies, Galactic sources and
/// cirrus). Order-of-magnitude reference value — not tabulated in the paper.
pub const N_UNIDENTIFIED_60UM: f64 = 3_000.0;

/// Approximate number of single-HCON sources in the IRAS Reject File acting
/// as the secondary pool whose surface density sets the chance rate.
/// Order-of-magnitude reference value — not tabulated in the paper.
pub const N_REJECT_SINGLE_HCON: f64 = 50_000.0;

/// Surface density (sources per arcmin²) of a catalogue of `n_sources`
/// spread over `sky_coverage` of the full sphere.
pub fn surface_density_per_arcmin2(n_sources: f64, sky_coverage: f64) -> f64 {
    let full_sphere_sr = 4.0 * std::f64::consts::PI;
    // Steradians → arcmin²: 1 sr = (180/π · 60)² arcmin².
    let arcmin2_per_sr = (180.0 / std::f64::consts::PI * 60.0).powi(2);
    let sky_arcmin2 = full_sphere_sr * sky_coverage * arcmin2_per_sr;
    n_sources / sky_arcmin2
}

/// Expected number of chance neighbours of a source within the annulus
/// `[r_min_arcmin, r_max_arcmin]`, given a Poisson field of surface density
/// `sigma_per_arcmin2`.
pub fn expected_neighbours(sigma_per_arcmin2: f64, r_min_arcmin: f64, r_max_arcmin: f64) -> f64 {
    let area = std::f64::consts::PI * (r_max_arcmin.powi(2) - r_min_arcmin.powi(2));
    sigma_per_arcmin2 * area
}

/// Probability of at least one chance neighbour in the annulus: `1 − e^{−μ}`
/// with μ = [`expected_neighbours`].
pub fn chance_alignment_probability(
    sigma_per_arcmin2: f64,
    r_min_arcmin: f64,
    r_max_arcmin: f64,
) -> f64 {
    let mu = expected_neighbours(sigma_per_arcmin2, r_min_arcmin, r_max_arcmin);
    1.0 - (-mu).exp()
}

/// Expected number of chance associations across an entire search of
/// `n_search_sources`, each tested for a neighbour in the annulus — the
/// quantity that produces the paper's "several hundred" candidate
/// associations.
pub fn expected_chance_associations(
    n_search_sources: f64,
    sigma_per_arcmin2: f64,
    r_min_arcmin: f64,
    r_max_arcmin: f64,
) -> f64 {
    n_search_sources * expected_neighbours(sigma_per_arcmin2, r_min_arcmin, r_max_arcmin)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sky_coverage() -> f64 {
        p9_2025_iras_akari::survey_model::IrasSurvey::default().sky_coverage
    }

    fn iras_sigma() -> f64 {
        surface_density_per_arcmin2(IRAS_FSC_N_60UM, sky_coverage())
    }

    fn reject_sigma() -> f64 {
        surface_density_per_arcmin2(N_REJECT_SINGLE_HCON, sky_coverage())
    }

    #[test]
    fn surface_density_is_about_1e_minus_3_per_arcmin2() {
        // 173k sources over ~96% of the sky ≈ 1.2e-3 per arcmin².
        let sigma = iras_sigma();
        assert!(
            (1.0e-3..1.4e-3).contains(&sigma),
            "σ = {sigma:e} per arcmin²"
        );
    }

    #[test]
    fn single_pair_chance_is_small() {
        // Within a tight 2 arcmin radius, a chance neighbour is rare (~1.5%):
        // it is the *number of trials*, not any single pairing, that inflates
        // the candidate count.
        let p = chance_alignment_probability(iras_sigma(), 0.0, 2.0);
        assert!(p < 0.05, "single-pair chance probability = {p:.4}");
    }

    #[test]
    fn wide_annulus_more_likely_than_narrow() {
        let sigma = iras_sigma();
        let narrow = chance_alignment_probability(sigma, 2.0, 10.0);
        let wide = chance_alignment_probability(sigma, 2.0, 35.0);
        assert!(wide > narrow);
    }

    #[test]
    fn probability_in_unit_interval() {
        let p = chance_alignment_probability(iras_sigma(), 2.0, 35.0);
        assert!((0.0..=1.0).contains(&p));
    }

    #[test]
    fn search_yields_hundreds_to_thousands_of_chance_associations() {
        // The real search: ~3000 unidentified 60 µm primaries, each tested
        // for a single-HCON Reject-File neighbour (sparser population) in the
        // 2–35 arcmin window. The raw geometric-coincidence count is ~10³
        // (here ≈ 4×10³ from the order-of-magnitude reference populations);
        // the paper's "several hundred candidate associations" are these
        // *after* Scanpi vetting thinned them. We pin the order of magnitude
        // — hundreds to a few thousand — not a precise count, since the
        // subset sizes are not tabulated in the paper. This is the residual
        // documented in REPRODUCTION_NOTES.md.
        let n = expected_chance_associations(N_UNIDENTIFIED_60UM, reject_sigma(), 2.0, 35.0);
        assert!(
            (1.0e2..1.0e4).contains(&n),
            "expected chance associations = {n:.0}"
        );
    }

    #[test]
    fn all_pairs_upper_bound_is_large() {
        // Testing every FSC source against the full FSC density is a loose
        // upper bound that runs to ~10⁵–10⁶: it is the restriction to the
        // unidentified primaries and the sparse Reject secondaries (above)
        // that brings the count down to the paper's several hundred.
        let n = expected_chance_associations(IRAS_FSC_N_60UM, iras_sigma(), 2.0, 35.0);
        assert!(n > 1.0e5, "all-pairs upper bound = {n:.0}");
    }

    #[test]
    fn expected_neighbours_scales_with_annulus_area() {
        let sigma = iras_sigma();
        let mu_small = expected_neighbours(sigma, 0.0, 10.0);
        let mu_big = expected_neighbours(sigma, 0.0, 20.0);
        // Area ratio is (20/10)² = 4.
        assert!((mu_big / mu_small - 4.0).abs() < 1e-9);
    }
}
