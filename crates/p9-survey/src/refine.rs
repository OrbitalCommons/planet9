//! Refine the Planet Nine position prior by removing sky that *other* assets
//! have already effectively covered — so JBT/SPENCER is pointed only where it
//! adds something new.
//!
//! Two cuts:
//!
//! 1. **Cede Rubin's sky.** Rubin/LSST will survey the south to r≈24.5 (deeper
//!    than any plausible P9). Assume it gets anything it can reach, so drop all
//!    cells with Dec ≤ its northern limit (~+12°).
//!
//! 2. **Survive prior optical surveys.** ZTF (r≈20.5) and Pan-STARRS (r≈21.5)
//!    have already imaged the northern sky; DES (r≈23.8) the south. If a survey
//!    covered a cell to a depth fainter than the source would be there, a P9
//!    would already have been found — so that cell survives only with
//!    probability `1 − coverage`. The product over surveys is the chance P9 is
//!    *still* there and *un-found*.
//!
//! The crucial consequence: a **bright** P9 (the close 2021/2024 solutions,
//! V≈19) in the north is mostly excluded by ZTF/PS1 already, so JBT adds little
//! there. JBT's unique value is a **faint/distant** P9 (V≳22, the wide 2016-era
//! solutions) in the north — fainter than ZTF/PS1 and out of Rubin's reach: a
//! region nothing else covers deep enough.
//!
//! Survey footprints/depths reuse `p9_core::analysis::surveys`; the refined map
//! is intentionally *not* renormalized — its sum is the residual probability
//! mass that is genuinely JBT's to chase.

use crate::schema::{Footprint, SkyGrid};
use p9_core::analysis::surveys::limiting_magnitude;
use p9_core::constants::DEG2RAD;
use p9_core::coords::sky::equatorial_to_galactic_matrix;

/// Rubin/LSST northern survey limit (deg); cells at or below are ceded.
pub const RUBIN_DEC_MAX_DEG: f64 = 12.0;

/// A survey that has already observed part of the sky, with the depth it
/// reached and how completely it covered its footprint.
#[derive(Debug, Clone)]
pub struct PriorSurvey {
    pub name: &'static str,
    pub footprint: Footprint,
    pub limiting_mag: f64,
}

/// The optical surveys whose non-detections constrain P9's current sky cell.
/// (WISE is omitted: a cold P9's thermal flux at W1 is negligible — see
/// `p9-2018-wise-search` — so it excludes essentially nothing here.)
pub fn prior_surveys() -> Vec<PriorSurvey> {
    vec![
        PriorSurvey {
            name: "ZTF",
            footprint: Footprint {
                dec_min_deg: -30.0,
                dec_max_deg: 90.0,
                galactic_lat_min_deg: 10.0,
                coverage_fraction: 0.70,
            },
            limiting_mag: limiting_magnitude("ZTF").unwrap_or(20.5),
        },
        PriorSurvey {
            name: "PS1 3pi",
            footprint: Footprint {
                dec_min_deg: -30.0,
                dec_max_deg: 90.0,
                galactic_lat_min_deg: 10.0,
                coverage_fraction: 0.75,
            },
            limiting_mag: limiting_magnitude("PS1 3pi").unwrap_or(21.5),
        },
        PriorSurvey {
            name: "DES",
            footprint: Footprint {
                dec_min_deg: -70.0,
                dec_max_deg: 5.0,
                galactic_lat_min_deg: 0.0,
                coverage_fraction: 0.30,
            },
            limiting_mag: limiting_magnitude("DES").unwrap_or(23.8),
        },
    ]
}

/// Galactic latitude (deg) of an equatorial sky position.
pub fn galactic_latitude_deg(ra_deg: f64, dec_deg: f64) -> f64 {
    let (ra, dec) = (ra_deg * DEG2RAD, dec_deg * DEG2RAD);
    let v = nalgebra::Vector3::new(dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin());
    let g = equatorial_to_galactic_matrix() * v;
    (g.z / g.norm()).asin().to_degrees()
}

/// Probability a source of magnitude `v` at (`ra`,`dec`) is still un-found by
/// the prior surveys: ∏ (1 − coverage) over surveys that cover the cell to a
/// depth fainter than `v`.
pub fn survival_probability(v: f64, ra_deg: f64, dec_deg: f64, surveys: &[PriorSurvey]) -> f64 {
    let gal_b = galactic_latitude_deg(ra_deg, dec_deg);
    let mut s = 1.0;
    for sv in surveys {
        let covered = sv.footprint.accepts(dec_deg, gal_b);
        let detectable = v <= sv.limiting_mag;
        if covered && detectable {
            s *= 1.0 - sv.footprint.coverage_fraction;
        }
    }
    s
}

/// Refine a position-probability grid for a source of representative magnitude
/// `source_mag`: cede Rubin's sky (Dec ≤ `RUBIN_DEC_MAX_DEG`) and weight each
/// surviving cell by `survival_probability`. Returns an *un-normalized* map
/// whose sum is the residual probability mass left for JBT.
pub fn refine(prob: &[f64], grid: &SkyGrid, source_mag: f64, surveys: &[PriorSurvey]) -> Vec<f64> {
    prob.iter()
        .enumerate()
        .map(|(idx, &p)| {
            let dec = grid.dec_center(idx / grid.n_ra);
            if dec <= RUBIN_DEC_MAX_DEG {
                return 0.0; // ceded to Rubin
            }
            let ra = grid.ra_min_deg + (idx % grid.n_ra) as f64 * grid.dra() + grid.dra() / 2.0;
            p * survival_probability(source_mag, ra, dec, surveys)
        })
        .collect()
}

/// Total residual probability mass after refinement (the size of JBT's niche).
pub fn residual_mass(
    prob: &[f64],
    grid: &SkyGrid,
    source_mag: f64,
    surveys: &[PriorSurvey],
) -> f64 {
    refine(prob, grid, source_mag, surveys).iter().sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{default_grid, run_survey};

    #[test]
    fn faint_source_survives_shallow_surveys() {
        let s = prior_surveys();
        // A V=22.7 source in the northern ZTF/PS1 sky: fainter than both
        // (20.5, 21.5) → not excluded by them → survives fully.
        let surv = survival_probability(22.7, 60.0, 20.0, &s);
        assert!(surv > 0.99, "faint north survival = {surv}");
        // A bright V=19.6 source in the same cell IS within ZTF+PS1 depth →
        // heavily suppressed.
        let bright = survival_probability(19.6, 60.0, 20.0, &s);
        assert!(bright < 0.1, "bright north survival = {bright}");
    }

    #[test]
    fn rubin_cede_zeros_the_south() {
        let d = run_survey(40_000, 7);
        let grid = default_grid();
        let s = prior_surveys();
        let refined = refine(&d.studies[3].prob, &grid, 19.6, &s);
        for (idx, &p) in refined.iter().enumerate() {
            if grid.dec_center(idx / grid.n_ra) <= RUBIN_DEC_MAX_DEG {
                assert_eq!(p, 0.0);
            }
        }
    }

    #[test]
    fn distant_solution_has_larger_residual_niche() {
        let d = run_survey(60_000, 11);
        let grid = default_grid();
        let s = prior_surveys();
        let bright = d
            .studies
            .iter()
            .find(|x| x.solution.name.contains("Siraj"))
            .unwrap();
        let faint = d
            .studies
            .iter()
            .find(|x| x.solution.name.contains("2016 Batygin & Brown (nominal)"))
            .unwrap();
        let r_bright = residual_mass(&bright.prob, &grid, bright.v_median, &s);
        let r_faint = residual_mass(&faint.prob, &grid, faint.v_median, &s);
        // The faint/distant solution leaves JBT a bigger un-covered niche than
        // the bright close one (which Rubin + ZTF + PS1 mostly mop up).
        assert!(
            r_faint > r_bright,
            "faint residual {r_faint} should exceed bright {r_bright}"
        );
    }
}
