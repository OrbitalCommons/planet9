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
//!
//! JBT 0.5 m + SPENCER parameters are derived from the `focalplane` simulator
//! (CosmicFrontierLabs): the default "Cosmic Frontier JBT .5m" telescope
//! (0.485 m aperture, f/12.3 → 5.97 m focal length) with the SPENCER imager
//! (4× Sony IMX455, 9568×6380 px at 3.76 µm, read noise 1.58 e⁻, dark
//! 0.0046 e⁻/s/px). That gives a 0.130″/px plate scale, ~0.080 deg² per chip
//! → ~0.32 deg² instantaneous field. A photon-budget SNR=5 calculation against
//! the space/zodiacal sky (V≈22.5/arcsec²) yields V≈23.8 in a single 300 s
//! visit, ~25 for an ~1 h per-field stack, and ~26 for a deep multi-epoch
//! shift-stack. So the limiter is not depth (P9 is V≈19–23) but the small
//! field: tiling the multi-thousand-deg² P9 prior is the constraint.

use crate::schema::Telescope;
use p9_core::analysis::surveys::Footprint;

/// JBT 0.5 m + SPENCER derived optics (see module docs / `focalplane`).
pub const JBT_APERTURE_M: f64 = 0.485;
pub const JBT_F_NUMBER: f64 = 12.3;
pub const JBT_PLATE_SCALE_ARCSEC_PER_PX: f64 = 0.130;
/// Instantaneous field of the 4-chip SPENCER array (deg²).
pub const SPENCER_FIELD_DEG2: f64 = 0.32;

/// The surveyed instruments: the Rubin/LSST ground baseline and JBT/SPENCER.
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
            name: "JBT 0.5 m + SPENCER (space)".to_string(),
            band: "broadband visible (V proxy)".to_string(),
            // Operational per-field depth (a few stacked 300 s visits). Single
            // 300 s visit ≈ 23.8; deep multi-epoch shift-stack ≈ 26.
            limiting_mag: 24.5,
            footprint: Footprint {
                dec_min_deg: -90.0,
                dec_max_deg: 90.0,
                galactic_lat_min_deg: 10.0,
                // The 0.32 deg² field must tile the multi-thousand-deg² P9
                // prior. Feasible (~months at ~0.3 deg²/visit) for a dedicated
                // campaign over the high-probability region; this fraction is
                // the dominant assumption — the limiter is sky coverage, not
                // depth. Galactic-plane crowding excluded (|b| > 10°).
                coverage_fraction: 0.70,
            },
            space_based: true,
            note: "Cosmic Frontier JBT 0.5 m (0.485 m, f/12.3) + SPENCER imager \
                   (4× Sony IMX455, 0.130″/px, ~0.32 deg² field), from the \
                   focalplane simulator. Depth is ample (V≈24–26 stacked vs P9 \
                   V≈19–23); detection is gated by how much of the prior region \
                   the small field can tile (coverage_fraction)."
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
    fn jbt_is_space_all_sky() {
        let t = catalog()
            .into_iter()
            .find(|t| t.name.contains("JBT"))
            .unwrap();
        assert!(t.space_based);
        assert!(t.footprint.accepts(80.0, 30.0)); // northern, allowed from space
        assert!(!t.footprint.accepts(80.0, 5.0)); // galactic plane excluded
                                                  // Deep enough that every surveyed P9 solution (V ≈ 19–23) is bright
                                                  // enough; coverage of the small field is the limiter.
        assert!(t.limiting_mag > 23.0);
    }
}
