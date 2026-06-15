//! Cross-paper survey of the Planet Nine prediction.
//!
//! This crate is a "survey paper" in code: it gathers the Planet Nine orbit
//! solutions reproduced elsewhere in the workspace (the Caltech 2016–2021
//! lineage from `p9_core::types::P9Params`, plus the independent Siraj, Chyba
//! & Tremaine 2024 solution from `p9-2024-siraj-orbit`), and for each one
//! computes — purely from p9-core primitives — where on the sky the planet
//! most probably is *now* (dwell-time-weighted over the unknown orbital phase
//! and the solution's element uncertainties), how bright it would be
//! (reflected-light V), and whether a given survey could catch it.
//!
//! The `survey` binary writes the whole result as a strongly-typed JSON
//! [`schema::SurveyDataset`]; `scripts/plot_survey.py` reads that exact schema
//! and renders the comparison heatmap. The chart is therefore a direct
//! rendering of numerical outputs of the reproduced papers — nothing is drawn
//! that the Rust side did not compute.

pub mod ephemeris;
pub mod parallax;
pub mod plan;
pub mod refine;
pub mod schema;
pub mod skymap;
pub mod studies;
pub mod telescope;

use schema::{Constraints, SkyGrid, SurveyDataset, TelescopeResult, SCHEMA_VERSION};

/// The fixed RA/Dec grid used for every probability map: full RA, ±60° Dec, 3°
/// cells (P9's declination never leaves this band for any surveyed solution).
pub fn default_grid() -> SkyGrid {
    SkyGrid {
        ra_min_deg: 0.0,
        ra_max_deg: 360.0,
        n_ra: 120,
        dec_min_deg: -60.0,
        dec_max_deg: 60.0,
        n_dec: 40,
    }
}

/// Run the full survey: every study against every telescope.
pub fn run_survey(n_samples: usize, seed: u64) -> SurveyDataset {
    let grid = default_grid();
    let solutions = studies::catalog();
    let telescopes = telescope::catalog();

    let mut studies_out = Vec::with_capacity(solutions.len());
    let mut per_telescope: Vec<Vec<schema::StudyDetection>> =
        vec![Vec::with_capacity(solutions.len()); telescopes.len()];

    for (idx, sol) in solutions.iter().enumerate() {
        let (result, detections) = skymap::run_study(
            sol,
            &grid,
            &telescopes,
            n_samples,
            seed.wrapping_add(idx as u64),
        );
        for (k, d) in detections.into_iter().enumerate() {
            per_telescope[k].push(d);
        }
        studies_out.push(result);
    }

    let telescopes_out = telescopes
        .into_iter()
        .zip(per_telescope)
        .map(|(telescope, per_study)| TelescopeResult {
            telescope,
            per_study,
        })
        .collect();

    let (nu_lo, nu_hi) = ephemeris::favored_interval_deg();
    let constraints = Constraints {
        rubin_dec_max_deg: refine::RUBIN_DEC_MAX_DEG,
        favored_nu_lo_deg: nu_lo,
        favored_nu_hi_deg: nu_hi,
        // The favored-ν zone belongs to the orbit Fienga et al. actually fit
        // (the Brown & Batygin geometry), so map it through that orbit.
        favored_arc: ephemeris::favored_arc(&p9_2016_cassini_ranging::brown_batygin_orbit(), 64),
    };

    SurveyDataset {
        schema_version: SCHEMA_VERSION,
        generated_by: "p9-survey".to_string(),
        samples_per_study: n_samples,
        rng_seed: seed,
        grid,
        studies: studies_out,
        telescopes: telescopes_out,
        overlays: skymap::overlays(360),
        constraints,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_survey_is_well_formed() {
        let d = run_survey(20_000, 2026);
        assert_eq!(d.schema_version, SCHEMA_VERSION);
        assert_eq!(d.studies.len(), 5);
        assert_eq!(d.telescopes.len(), 2);
        // every telescope carries a detection entry for every study
        for t in &d.telescopes {
            assert_eq!(t.per_study.len(), d.studies.len());
        }
        // round-trips through JSON
        let json = serde_json::to_string(&d).unwrap();
        let back: SurveyDataset = serde_json::from_str(&json).unwrap();
        assert_eq!(back.studies.len(), d.studies.len());
    }
}
