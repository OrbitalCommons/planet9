//! The "what has been observed" side: the wide optical surveys that have
//! actually hunted for Planet Nine, each as a footprint (declination band +
//! galactic-plane cut + areal coverage) and a limiting depth.
//!
//! Depths come from the single-source-of-truth table
//! [`p9_core::analysis::surveys::SURVEY_DEPTHS`]; the δ > −30° northern cut is
//! [`p9_core::analysis::surveys::NORTHERN_SURVEY_DEC_LIMIT_DEG`]. Footprint
//! geometry and areal coverage are the literature values cited per entry.
//!
//! Combining them gives, per sky direction, the *deepest* search that has
//! reached it — the survey-coverage hull. For a predicted Planet Nine
//! magnitude `v` at a direction, the probability the body would already have
//! been caught is a per-survey OR:
//!
//! ```text
//!   P_excl(v) = 1 − Π_s ( 1 − cov_s · 1[depth_s ≥ v] )
//! ```
//!
//! over the surveys `s` whose footprint contains that direction.

use p9_2020_tess_shiftstack::published::{SENSITIVITY_DISTANCE_AU, SENSITIVITY_V_LIMIT};
use p9_core::analysis::surveys::{
    limiting_magnitude, Footprint, NORTHERN_SURVEY_DEC_LIMIT_DEG, ZTF_DEC_LIMIT_DEG,
};

/// One archival wide-field optical search for Planet Nine.
#[derive(Debug, Clone)]
pub struct RealSurvey {
    pub name: &'static str,
    /// Limiting magnitude (r/V; the two are within the bookkeeping error of
    /// this map for a neutral-ish giant planet).
    pub depth: f64,
    pub footprint: Footprint,
    /// Optional right-ascension EXCLUSION interval (lo, hi) in degrees: the
    /// survey covers RA ≤ lo or RA ≥ hi (a wrap-around acceptance window).
    /// `None` = all RA. Needed because the core `Footprint` carries no RA
    /// term — without it DES's 5,000 deg² was painted across the whole
    /// δ ∈ [−65, +5] band (~14,000 deg² it never imaged).
    pub ra_exclusion_deg: Option<(f64, f64)>,
    /// Maximum heliocentric distance (au) to which this search's method is
    /// valid, or `None` for stationary-depth catalog searches. Shift-stack
    /// searches only probe the trial tracks they ran: TESS's grid stops at
    /// 150 au, so its non-detection says nothing about a body at 400 au even
    /// if that body is brighter than the stacked depth.
    pub max_distance_au: Option<f64>,
    pub reference: &'static str,
}

impl RealSurvey {
    /// Does this survey's geometry contain the direction (ra, dec, b)?
    pub fn accepts(&self, ra_deg: f64, dec_deg: f64, gal_b_deg: f64) -> bool {
        if !self.footprint.accepts(dec_deg, gal_b_deg) {
            return false;
        }
        match self.ra_exclusion_deg {
            Some((lo, hi)) => {
                let ra = ra_deg.rem_euclid(360.0);
                ra <= lo || ra >= hi
            }
            None => true,
        }
    }

    /// Does this survey constrain a body at direction (ra, dec, b) and
    /// heliocentric distance `dist_au`? Geometry plus the method's distance
    /// validity; a NaN distance fails any cap (unknown distance is treated as
    /// hull-typical, i.e. far).
    pub fn constrains(&self, ra_deg: f64, dec_deg: f64, gal_b_deg: f64, dist_au: f64) -> bool {
        self.accepts(ra_deg, dec_deg, gal_b_deg)
            && self.max_distance_au.is_none_or(|cap| dist_au <= cap)
    }
}

/// Look up a depth from the shared table, panicking if the name is missing so a
/// typo can never silently zero out a survey.
fn depth(name: &str) -> f64 {
    limiting_magnitude(name).unwrap_or_else(|| panic!("no depth for survey {name:?}"))
}

/// The wide optical surveys whose non-detections constrain Planet Nine.
///
/// These are the optical reflected-light searches; the thermal/IR channels
/// (WISE, AKARI) constrain a different (size, distance) slice and are tracked
/// by `p9-viability`, so they are deliberately out of this optical hull.
///
/// TESS (Rice & Laughlin 2020) is included with a DISTANCE cap: a shift-stack
/// search only constrains bodies whose motion matches a searched trial track,
/// and that search's track grid stops at heliocentric distances d ≤ 150 au
/// (`p9_2020_tess_shiftstack::published::SENSITIVITY_DISTANCE_AU`). Most of
/// the posterior ensemble sits far beyond that, but the low-a/high-e tails
/// (e.g. Siraj 2024 draws near perihelion) do dip inside — so TESS is a real
/// exclusion channel for exactly those draws and must not claim the rest.
/// Listing it uncapped at its stacked depth would falsely mark directions
/// excluded where no track at Planet Nine's sky rate was ever searched.
pub fn real_surveys() -> Vec<RealSurvey> {
    vec![
        // Catalina Real-Time Transient Survey: shallow but very wide, both
        // hemispheres away from the plane.
        RealSurvey {
            name: "CRTS",
            depth: depth("CRTS"),
            footprint: Footprint {
                dec_min_deg: -40.0,
                dec_max_deg: 70.0,
                galactic_lat_min_deg: 10.0,
                coverage_fraction: 0.60,
            },
            ra_exclusion_deg: None,
            max_distance_au: None,
            reference: "Drake et al. (2009); footprint approx.",
        },
        // ZTF (Palomar): δ > −30°, used by Brown & Batygin (2022) for a P9
        // shift-stack non-detection.
        // Plane treatment: ZTF images |b| < 10° but tracklet linking there
        // is crowding-limited; the hull takes the CONSERVATIVE hard hole
        // (unexcluded), while p9-survey::refine models it as
        // shallower-but-imaged. Same survey, two deliberate fidelities.
        RealSurvey {
            name: "ZTF",
            depth: depth("ZTF"),
            footprint: Footprint {
                dec_min_deg: ZTF_DEC_LIMIT_DEG,
                dec_max_deg: 90.0,
                galactic_lat_min_deg: 10.0,
                coverage_fraction: 0.70,
            },
            ra_exclusion_deg: None,
            max_distance_au: None,
            reference: "Brown & Batygin (2022), ZTF non-detection",
        },
        // Pan-STARRS1 3π: δ > −30°, near-all-northern-sky to r ≈ 21.5.
        RealSurvey {
            name: "PS1 3pi",
            depth: depth("PS1 3pi"),
            footprint: Footprint {
                dec_min_deg: NORTHERN_SURVEY_DEC_LIMIT_DEG,
                dec_max_deg: 90.0,
                galactic_lat_min_deg: 0.0,
                coverage_fraction: 0.75,
            },
            ra_exclusion_deg: None,
            max_distance_au: None,
            reference: "Chambers et al. (2016); Belyakov et al. (2024)",
        },
        // Dark Energy Survey: deep but small southern footprint (~5000 deg²
        // of the ~8,400 deg² south-galactic-cap box), matching the DesSgc
        // geometry in p9-survey::refine (RA ≤ 100° or ≥ 300°, δ ∈ [−68, +5],
        // |b| ≥ 18°). The old entry had no RA term, so cell_coverage reported
        // DES as the deepest survey at EVERY RA in the declination band and
        // its areal coverage was diluted to 0.13 to compensate.
        RealSurvey {
            name: "DES",
            depth: depth("DES"),
            footprint: Footprint {
                dec_min_deg: -68.0,
                dec_max_deg: 5.0,
                galactic_lat_min_deg: 18.0,
                coverage_fraction: 0.60,
            },
            ra_exclusion_deg: Some((100.0, 300.0)),
            max_distance_au: None,
            reference: "Bernardinelli et al. (2022); SGC box as p9-survey::refine",
        },
        // TESS shift-stack (Rice & Laughlin 2020): the blind search's
        // demonstrated V < 21 sensitivity applies only along the trial
        // tracks it ran, d <= 150 au — hence the distance cap. The blind
        // sectors (18-19) are two ~2,300 deg^2 northern strips; the coarse
        // dec band + areal fraction below carries that ~4,600 deg^2.
        RealSurvey {
            name: "TESS",
            depth: SENSITIVITY_V_LIMIT,
            footprint: Footprint {
                dec_min_deg: 0.0,
                dec_max_deg: 90.0,
                galactic_lat_min_deg: 0.0,
                coverage_fraction: 0.20,
            },
            ra_exclusion_deg: None,
            max_distance_au: Some(SENSITIVITY_DISTANCE_AU),
            reference: "Rice & Laughlin (2020), PSJ 1, 81; sectors 18-19",
        },
    ]
}

/// Per-direction coverage summary at one sky cell.
pub struct CellCoverage {
    /// Deepest limiting magnitude reaching this direction (−∞ proxy: f64::NAN
    /// when no survey covers it).
    pub best_depth: Option<f64>,
    /// Name of the survey supplying `best_depth`.
    pub best_survey: Option<&'static str>,
}

/// Deepest survey reaching a direction for a body at heliocentric distance
/// `dist_au` (areal coverage of that survey is returned alongside, not
/// applied). The distance gates the method-validity caps: pass the cell's
/// posterior mean distance (a NaN — empty cell — fails any cap).
pub fn cell_coverage(
    surveys: &[RealSurvey],
    ra_deg: f64,
    dec_deg: f64,
    gal_b_deg: f64,
    dist_au: f64,
) -> CellCoverage {
    let mut best: Option<(&RealSurvey,)> = None;
    for s in surveys {
        if s.constrains(ra_deg, dec_deg, gal_b_deg, dist_au) {
            match best {
                Some((b,)) if b.depth >= s.depth => {}
                _ => best = Some((s,)),
            }
        }
    }
    match best {
        Some((s,)) => CellCoverage {
            best_depth: Some(s.depth),
            best_survey: Some(s.name),
        },
        None => CellCoverage {
            best_depth: None,
            best_survey: None,
        },
    }
}

/// Probability a body of apparent magnitude `v` at this direction and
/// heliocentric distance `dist_au` would already have been caught: the
/// per-survey OR of (covers this direction and distance) × (deep enough) ×
/// (areal coverage).
pub fn exclusion_probability(
    surveys: &[RealSurvey],
    ra_deg: f64,
    dec_deg: f64,
    gal_b_deg: f64,
    dist_au: f64,
    v: f64,
) -> f64 {
    let mut survive = 1.0_f64;
    for s in surveys {
        if s.constrains(ra_deg, dec_deg, gal_b_deg, dist_au) && s.depth >= v {
            survive *= 1.0 - s.footprint.coverage_fraction;
        }
    }
    1.0 - survive
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn depths_come_from_the_shared_table() {
        let s = real_surveys();
        let ztf = s.iter().find(|x| x.name == "ZTF").unwrap();
        assert_eq!(ztf.depth, limiting_magnitude("ZTF").unwrap());
        let des = s.iter().find(|x| x.name == "DES").unwrap();
        assert_eq!(des.depth, 23.8);
    }

    #[test]
    fn deepest_survey_wins_the_hull() {
        // A southern, off-plane direction inside the DES SGC box (RA 30°):
        // DES (23.8) is deepest.
        let c = cell_coverage(&real_surveys(), 30.0, -40.0, 40.0, 400.0);
        assert_eq!(c.best_survey, Some("DES"));
        assert_eq!(c.best_depth, Some(23.8));
    }

    #[test]
    fn des_does_not_cover_outside_its_ra_window() {
        // Same declination band but RA 200° — the ~14,000 deg² the old
        // RA-less footprint wrongly painted as DES-deep. CRTS (19.5) is the
        // only survey there.
        let c = cell_coverage(&real_surveys(), 200.0, -40.0, 40.0, 400.0);
        assert_eq!(c.best_survey, Some("CRTS"));
        // And the RA wrap: RA 330° is back inside the SGC box.
        let c2 = cell_coverage(&real_surveys(), 330.0, -40.0, 40.0, 400.0);
        assert_eq!(c2.best_survey, Some("DES"));
    }

    #[test]
    fn northern_direction_falls_to_ps1() {
        // Far north, on no southern survey: PS1 (21.5) is the deepest all-sky
        // northern reach (deeper than ZTF 20.5).
        let c = cell_coverage(&real_surveys(), 30.0, 60.0, 40.0, 400.0);
        assert_eq!(c.best_survey, Some("PS1 3pi"));
    }

    #[test]
    fn faint_planet_is_less_excluded_than_bright() {
        let s = real_surveys();
        // A southern, off-plane direction in the SGC box: CRTS + DES cover it.
        let p_bright = exclusion_probability(&s, 30.0, -40.0, 40.0, 400.0, 19.0); // CRTS + DES
        let p_faint = exclusion_probability(&s, 30.0, -40.0, 40.0, 400.0, 23.0); // only DES reaches
        assert!(p_bright > p_faint, "{p_bright} !> {p_faint}");
        // A planet fainter than every survey here cannot be excluded.
        let p_too_faint = exclusion_probability(&s, 30.0, -40.0, 40.0, 400.0, 25.0);
        assert_eq!(p_too_faint, 0.0);
        // Outside the RA window, the faint planet is not excluded at all.
        let p_faint_out = exclusion_probability(&s, 200.0, -40.0, 40.0, 400.0, 23.0);
        assert_eq!(p_faint_out, 0.0);
    }

    #[test]
    fn tess_constrains_only_its_searched_distances() {
        let s = real_surveys();
        let tess = s.iter().find(|x| x.name == "TESS").unwrap();
        assert_eq!(tess.max_distance_au, Some(SENSITIVITY_DISTANCE_AU));
        assert_eq!(tess.depth, SENSITIVITY_V_LIMIT);
        // A northern direction, body inside the searched track grid: TESS
        // constrains it...
        assert!(tess.constrains(30.0, 40.0, 30.0, 100.0));
        // ...but the same direction at hull-typical distance is NOT
        // constrained (no trial track was run there), and an unknown (NaN)
        // distance is treated as far.
        assert!(!tess.constrains(30.0, 40.0, 30.0, 400.0));
        assert!(!tess.constrains(30.0, 40.0, 30.0, f64::NAN));

        // The distance cap flows through exclusion: at V=20.7 north of the
        // CRTS band, PS1 (21.5, cov 0.75) reaches at any distance and ZTF
        // (20.5) at none — TESS's 20% only ORs in when the body is inside
        // its searched tracks.
        let p_close = exclusion_probability(&s, 30.0, 80.0, 30.0, 100.0, 20.7);
        let p_far = exclusion_probability(&s, 30.0, 80.0, 30.0, 400.0, 20.7);
        assert_eq!(p_far, 0.75);
        assert!((p_close - (1.0 - 0.25 * 0.80)).abs() < 1e-12, "{p_close}");

        // And through the hull map: at 100 au TESS's V<21 beats nothing
        // deeper than PS1 21.5 — PS1 still wins where both apply, so the cap
        // only ever ADDS coverage, never displaces a deeper survey.
        let c = cell_coverage(&s, 30.0, 60.0, 40.0, 100.0);
        assert_eq!(c.best_survey, Some("PS1 3pi"));
    }
}
