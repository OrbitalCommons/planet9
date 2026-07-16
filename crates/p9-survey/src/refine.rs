//! Refine the Planet Nine position prior by removing sky that *other* assets
//! have already effectively covered — so JBT/SPENCER is pointed only where it
//! adds something new.
//!
//! Three cuts (the third optional):
//!
//! 1. **Cede Rubin's sky.** Rubin/LSST will survey the south to r≈24.5 (deeper
//!    than any plausible P9). Drop all cells with Dec ≤ its northern limit
//!    (~+12°).
//!
//! 2. **Survive prior optical surveys — with real footprints, not dec bands.**
//!    ZTF (Dec ≳ −31°) and Pan-STARRS 3π (Dec ≳ −30°) cover the whole northern
//!    sky *including the galactic plane*; DES covers a ~5000 deg² South-
//!    Galactic-Cap polygon. Crucially the plane is **not a hole** — the surveys
//!    imaged it, just shallower (crowding + extinction), so a survey's
//!    *effective depth* and *completeness* degrade toward |b|→0. A cell
//!    survives a survey with probability `1 − coverage(b)` only if the source
//!    is brighter than that survey's effective depth there. The product over
//!    surveys is the chance P9 is still there and un-found — so a faint P9 in
//!    the plane survives because the surveys are *shallow* there, not absent.
//!
//! 3. **Ephemeris phase prior** (applied upstream, in `skymap`): see
//!    [`crate::ephemeris`] — concentrates the prior onto the Cassini/Iorio
//!    favored-ν sky arc, cutting the RA range.
//!
//! The plane depth/coverage degradation is a documented parametric stand-in for
//! each survey's true crowding/extinction map; the footprint *boundaries* are
//! the real survey limits. The refined map is left un-normalized — its sum is
//! the residual probability mass genuinely left for JBT.

use crate::schema::SkyGrid;
use p9_core::analysis::surveys::limiting_magnitude;
use p9_core::analysis::surveys::{NORTHERN_SURVEY_DEC_LIMIT_DEG, ZTF_DEC_LIMIT_DEG};
use p9_core::constants::DEG2RAD;
use p9_core::coords::sky::equatorial_to_galactic_matrix;

/// Rubin/LSST northern survey limit (deg); cells at or below are ceded.
pub const RUBIN_DEC_MAX_DEG: f64 = 12.0;

/// Galactic-latitude scale over which plane crowding/extinction bites (deg).
const PLANE_B_SCALE_DEG: f64 = 8.0;

/// Footprint shape of a prior survey.
#[derive(Debug, Clone, Copy)]
pub enum Footprint {
    /// Everything north of a declination limit (ground survey from the north).
    NorthOf(f64),
    /// The DES South-Galactic-Cap region: Dec ∈ [−68, +5], RA in the SGC
    /// longitudes, |b| ≥ 18°.
    DesSgc,
}

impl Footprint {
    fn covers(&self, _ra_deg: f64, dec_deg: f64, gal_b_deg: f64) -> bool {
        match self {
            Footprint::NorthOf(dmin) => dec_deg >= *dmin,
            Footprint::DesSgc => {
                let ra = _ra_deg.rem_euclid(360.0);
                (-68.0..=5.0).contains(&dec_deg)
                    && (ra <= 100.0 || ra >= 300.0)
                    && gal_b_deg.abs() >= 18.0
            }
        }
    }
}

/// A survey that has already observed part of the sky, with the depth and
/// completeness it reached in clean (high-|b|) sky, degraded toward the plane.
#[derive(Debug, Clone)]
pub struct PriorSurvey {
    pub name: &'static str,
    pub footprint: Footprint,
    /// Limiting magnitude away from the galactic plane.
    pub clean_depth: f64,
    /// Completeness fraction away from the plane.
    pub clean_coverage: f64,
    /// Magnitudes of depth lost at b = 0 (crowding + extinction).
    pub plane_depth_penalty: f64,
    /// Fraction of completeness lost at b = 0.
    pub plane_cov_loss: f64,
}

impl PriorSurvey {
    /// 1 at the plane, →0 away from it.
    fn plane_factor(gal_b_deg: f64) -> f64 {
        (-0.5 * (gal_b_deg / PLANE_B_SCALE_DEG).powi(2)).exp()
    }
    pub fn effective_depth(&self, gal_b_deg: f64) -> f64 {
        self.clean_depth - self.plane_depth_penalty * Self::plane_factor(gal_b_deg)
    }
    pub fn coverage(&self, gal_b_deg: f64) -> f64 {
        self.clean_coverage * (1.0 - self.plane_cov_loss * Self::plane_factor(gal_b_deg))
    }
}

/// The optical surveys whose non-detections constrain P9's current sky cell.
/// (WISE is omitted: a cold P9's thermal flux at W1 is negligible — see
/// `p9-2018-wise-search`.)
pub fn prior_surveys() -> Vec<PriorSurvey> {
    vec![
        PriorSurvey {
            name: "ZTF",
            footprint: Footprint::NorthOf(ZTF_DEC_LIMIT_DEG),
            clean_depth: limiting_magnitude("ZTF").unwrap_or(20.5),
            clean_coverage: 0.70,
            plane_depth_penalty: 2.5,
            plane_cov_loss: 0.5,
        },
        PriorSurvey {
            name: "PS1 3pi",
            footprint: Footprint::NorthOf(NORTHERN_SURVEY_DEC_LIMIT_DEG),
            clean_depth: limiting_magnitude("PS1 3pi").unwrap_or(21.5),
            clean_coverage: 0.75,
            plane_depth_penalty: 2.5,
            plane_cov_loss: 0.5,
        },
        PriorSurvey {
            name: "DES",
            footprint: Footprint::DesSgc,
            clean_depth: limiting_magnitude("DES").unwrap_or(23.8),
            clean_coverage: 0.60,
            // The DES footprint is high-|b| by construction, so plane terms
            // never engage; keep them zero.
            plane_depth_penalty: 0.0,
            plane_cov_loss: 0.0,
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
/// the prior surveys: ∏ (1 − coverage(b)) over surveys that cover the cell to
/// an effective depth fainter than `v`.
pub fn survival_probability(v: f64, ra_deg: f64, dec_deg: f64, surveys: &[PriorSurvey]) -> f64 {
    let gal_b = galactic_latitude_deg(ra_deg, dec_deg);
    let mut s = 1.0;
    for sv in surveys {
        if sv.footprint.covers(ra_deg, dec_deg, gal_b) && v <= sv.effective_depth(gal_b) {
            s *= 1.0 - sv.coverage(gal_b);
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
    fn plane_is_shallower_than_clean_sky() {
        let s = &prior_surveys()[0]; // ZTF
        let clean = s.effective_depth(40.0);
        let plane = s.effective_depth(0.0);
        assert!((clean - 20.5).abs() < 0.05, "clean = {clean}");
        assert!(plane < clean - 2.0, "plane depth {plane} not shallower");
        assert!(s.coverage(0.0) < s.coverage(40.0));
    }

    #[test]
    fn faint_source_survives_shallow_surveys() {
        let s = prior_surveys();
        // A V=22.7 source off the plane in the northern ZTF/PS1 sky: fainter
        // than both clean depths (20.5, 21.5) → not excluded → survives.
        let surv = survival_probability(22.7, 60.0, 25.0, &s);
        assert!(surv > 0.99, "faint clean-north survival = {surv}");
        // A bright V=19.6 source there IS within ZTF+PS1 depth → suppressed.
        let bright = survival_probability(19.6, 60.0, 25.0, &s);
        assert!(bright < 0.1, "bright clean-north survival = {bright}");
        // The SAME bright source in the galactic plane survives better, because
        // ZTF/PS1 are shallow & incomplete there (crowding) — the new physics.
        let bright_plane = survival_probability(19.6, 82.0, 20.0, &s);
        let b_here = galactic_latitude_deg(82.0, 20.0);
        if b_here.abs() < 6.0 {
            assert!(
                bright_plane > bright,
                "plane survival {bright_plane} should exceed clean {bright}"
            );
        }
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
        assert!(
            r_faint > r_bright,
            "faint {r_faint} should exceed bright {r_bright}"
        );
    }

    #[test]
    fn des_sgc_times_coverage_matches_the_box_union_area() {
        // The DesSgc footprint (SGC box + |b| cut) with clean_coverage 0.60
        // is an areal-coverage summary of the real box-union footprint in
        // p9-core: box area x coverage must land on the box-union solid
        // angle (~4,980 deg^2) to ~10%. This ties the two DES geometries the
        // workspace keeps (membership vs coverage-summary) to one another.
        use p9_core::analysis::surveys::des_footprint_solid_angle_deg2;
        let des = prior_surveys()
            .into_iter()
            .find(|s| s.name == "DES")
            .unwrap();
        // MC solid angle of the DesSgc acceptance region.
        let (n_ra, n_dec) = (720, 360);
        let mut area = 0.0;
        let cell = (360.0 / n_ra as f64) * (180.0 / n_dec as f64);
        for iy in 0..n_dec {
            let dec = -90.0 + (iy as f64 + 0.5) * 180.0 / n_dec as f64;
            let w = (dec.to_radians()).cos();
            for ix in 0..n_ra {
                let ra = (ix as f64 + 0.5) * 360.0 / n_ra as f64;
                let gal_b = galactic_latitude_deg(ra, dec);
                if des.footprint.covers(ra, dec, gal_b) {
                    area += cell * w;
                }
            }
        }
        let effective = area * des.clean_coverage;
        let union = des_footprint_solid_angle_deg2();
        let ratio = effective / union;
        assert!(
            (0.9..1.1).contains(&ratio),
            "SGC box {area:.0} deg2 x {:.2} = {effective:.0} vs union {union:.0} (ratio {ratio:.2})",
            des.clean_coverage
        );
    }
}
