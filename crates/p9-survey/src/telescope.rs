//! Telescope / survey catalog for the detectability assessment.
//!
//! Detection here is a reflected-light optical/NIR model: a sample is
//! "detected" if it lands in the imaged footprint *and* is brighter than the
//! limiting magnitude. Limiting depths are single-visit-to-stacked values for
//! a slow-moving (shift-and-stack-friendly) target; treat the magnitude as a
//! V-band reflected-light proxy (NIR instruments carry a documented color
//! caveat). Thermal-IR detectability of a *cold* Planet Nine is separately
//! shown to be negligible in `p9-2018-wise-search`, so this survey stays in
//! the reflected-light regime.

use crate::schema::{Footprint, Telescope};

/// The surveyed instruments: one ground baseline, two space assets.
pub fn catalog() -> Vec<Telescope> {
    vec![
        Telescope {
            name: "Rubin / LSST (ground baseline)".to_string(),
            band: "r (stacked)".to_string(),
            limiting_mag: 24.5,
            footprint: Footprint {
                dec_min_deg: -75.0,
                dec_max_deg: 12.0,
                galactic_lat_min_deg: 10.0,
                coverage_fraction: 0.85,
            },
            space_based: false,
            note: "8.4 m, ~18,000 deg² southern main survey; deep 10-yr stack \
                   reaches r≈27 but plane-avoiding. Ground comparison point."
                .to_string(),
        },
        Telescope {
            name: "Roman (space, NIR, pointed)".to_string(),
            band: "F146 (NIR proxy)".to_string(),
            limiting_mag: 26.0,
            footprint: Footprint {
                dec_min_deg: -90.0,
                dec_max_deg: 90.0,
                galactic_lat_min_deg: 0.0,
                coverage_fraction: 0.02,
            },
            space_based: true,
            note: "2.4 m, very deep but a tiny 0.28 deg² field — superb depth, \
                   negligible blind-survey sky coverage for a wide P9 search."
                .to_string(),
        },
        Telescope {
            name: "Proposed all-sky space survey".to_string(),
            band: "optical/NIR (proxy)".to_string(),
            limiting_mag: 26.0,
            footprint: Footprint {
                dec_min_deg: -90.0,
                dec_max_deg: 90.0,
                galactic_lat_min_deg: 5.0,
                coverage_fraction: 0.90,
            },
            space_based: true,
            note: "PLACEHOLDER specs for the upcoming mission — deep, near-all-\
                   sky from space (only the galactic plane still hurts). Tune \
                   limiting_mag / coverage / footprint to the real design."
                .to_string(),
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn footprint_dec_band() {
        let t = &catalog()[0]; // Rubin
        assert!(t.footprint.accepts(-30.0, 40.0)); // in band, off plane
        assert!(!t.footprint.accepts(60.0, 40.0)); // too far north
        assert!(!t.footprint.accepts(-30.0, 3.0)); // on galactic plane
    }

    #[test]
    fn space_survey_is_all_sky() {
        let t = catalog()
            .into_iter()
            .find(|t| t.name.contains("Proposed"))
            .unwrap();
        assert!(t.space_based);
        assert!(t.footprint.accepts(80.0, 30.0)); // northern, allowed from space
        assert!(t.footprint.coverage_fraction > 0.5);
    }
}
