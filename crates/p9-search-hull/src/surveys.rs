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

use p9_core::analysis::surveys::{limiting_magnitude, Footprint, NORTHERN_SURVEY_DEC_LIMIT_DEG};

/// One archival wide-field optical search for Planet Nine.
#[derive(Debug, Clone)]
pub struct RealSurvey {
    pub name: &'static str,
    /// Limiting magnitude (r/V; the two are within the bookkeeping error of
    /// this map for a neutral-ish giant planet).
    pub depth: f64,
    pub footprint: Footprint,
    pub reference: &'static str,
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
            reference: "Drake et al. (2009); footprint approx.",
        },
        // ZTF (Palomar): δ > −30°, used by Brown & Batygin (2022) for a P9
        // shift-stack non-detection.
        RealSurvey {
            name: "ZTF",
            depth: depth("ZTF"),
            footprint: Footprint {
                dec_min_deg: NORTHERN_SURVEY_DEC_LIMIT_DEG,
                dec_max_deg: 90.0,
                galactic_lat_min_deg: 10.0,
                coverage_fraction: 0.70,
            },
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
            reference: "Chambers et al. (2016); Belyakov et al. (2024)",
        },
        // Dark Energy Survey: deep but small southern footprint (~5000 deg²),
        // off the galactic plane.
        RealSurvey {
            name: "DES",
            depth: depth("DES"),
            footprint: Footprint {
                dec_min_deg: -65.0,
                dec_max_deg: 5.0,
                galactic_lat_min_deg: 10.0,
                coverage_fraction: 0.13,
            },
            reference: "Bernardinelli et al. (2022)",
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

/// Deepest survey reaching a direction (geometry only — the areal coverage of
/// that survey is returned alongside, not applied).
pub fn cell_coverage(surveys: &[RealSurvey], dec_deg: f64, gal_b_deg: f64) -> CellCoverage {
    let mut best: Option<(&RealSurvey,)> = None;
    for s in surveys {
        if s.footprint.accepts(dec_deg, gal_b_deg) {
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

/// Probability a body of apparent magnitude `v` at this direction would already
/// have been caught: the per-survey OR of (covered this direction) × (deep
/// enough) × (areal coverage).
pub fn exclusion_probability(surveys: &[RealSurvey], dec_deg: f64, gal_b_deg: f64, v: f64) -> f64 {
    let mut survive = 1.0_f64;
    for s in surveys {
        if s.footprint.accepts(dec_deg, gal_b_deg) && s.depth >= v {
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
        // A southern, off-plane direction inside DES: DES (23.8) is deepest.
        let c = cell_coverage(&real_surveys(), -40.0, 40.0);
        assert_eq!(c.best_survey, Some("DES"));
        assert_eq!(c.best_depth, Some(23.8));
    }

    #[test]
    fn northern_direction_falls_to_ps1() {
        // Far north, on no southern survey: PS1 (21.5) is the deepest all-sky
        // northern reach (deeper than ZTF 20.5).
        let c = cell_coverage(&real_surveys(), 60.0, 40.0);
        assert_eq!(c.best_survey, Some("PS1 3pi"));
    }

    #[test]
    fn faint_planet_is_less_excluded_than_bright() {
        let s = real_surveys();
        // A southern, off-plane direction: CRTS (19.5) and DES (23.8) cover it.
        let p_bright = exclusion_probability(&s, -40.0, 40.0, 19.0); // CRTS + DES
        let p_faint = exclusion_probability(&s, -40.0, 40.0, 23.0); // only DES reaches
        assert!(p_bright > p_faint, "{p_bright} !> {p_faint}");
        // A planet fainter than every survey here cannot be excluded.
        let p_too_faint = exclusion_probability(&s, -40.0, 40.0, 25.0);
        assert_eq!(p_too_faint, 0.0);
    }
}
