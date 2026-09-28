//! Distant objects found for free: how many trans-Neptunian objects a tile
//! yields at a given depth, how many of them are the distant ones worth a
//! paper, and how many Rubin will not already have.
//!
//! * **Sky density.** Cumulative counts on the ecliptic follow a broken power
//!   law, Σ(< m_R) = 10^{α₁(m_R − 23.4)} per deg² up to the break at
//!   m_R = 24.5 and slope α₂ beyond it (α₁ = 0.75, α₂ = 0.35; the range
//!   bracketed by Fraser et al. 2014 and Bernstein et al. 2004).
//! * **Latitude.** A cold component (σ = 2.5°, 35% of the on-ecliptic
//!   density) and a hot one (σ = 15°) after Brown (2001).
//! * **Distant fraction.** About 2% of a flux-limited TNO sample lies beyond
//!   60 AU (the OSSOS/DES detected-distance distributions); these are the
//!   detached and sednoid candidates.
//!
//! These are survey-planning numbers, good to a factor of ~2.

use p9_core::analysis::photometry::SOLAR_V_MINUS_R;

pub const ALPHA_BRIGHT: f64 = 0.75;
pub const ALPHA_FAINT: f64 = 0.35;
pub const R_UNIT_DENSITY: f64 = 23.4;
pub const R_BREAK: f64 = 24.5;
pub const COLD_FRACTION: f64 = 0.35;
pub const COLD_SIGMA_DEG: f64 = 2.5;
pub const HOT_SIGMA_DEG: f64 = 15.0;
pub const DISTANT_FRACTION: f64 = 0.02;
/// Typical TNO colour redder than the Sun (V − R ≈ 0.6).
pub const TNO_V_MINUS_R: f64 = SOLAR_V_MINUS_R + 0.24;

/// TNOs per deg² on the ecliptic brighter than `v_depth`.
pub fn ecliptic_density(v_depth: f64) -> f64 {
    let m = v_depth - TNO_V_MINUS_R;
    if m <= R_BREAK {
        10f64.powf(ALPHA_BRIGHT * (m - R_UNIT_DENSITY))
    } else {
        10f64.powf(ALPHA_BRIGHT * (R_BREAK - R_UNIT_DENSITY) + ALPHA_FAINT * (m - R_BREAK))
    }
}

/// Density relative to the ecliptic at ecliptic latitude `beta_deg`.
pub fn latitude_profile(beta_deg: f64) -> f64 {
    let g = |s: f64| (-0.5 * (beta_deg / s).powi(2)).exp();
    COLD_FRACTION * g(COLD_SIGMA_DEG) + (1.0 - COLD_FRACTION) * g(HOT_SIGMA_DEG)
}

/// Expected TNOs in `area_deg2` at latitude `beta_deg` to depth `v_depth`.
pub fn expected_tnos(area_deg2: f64, beta_deg: f64, v_depth: f64) -> f64 {
    area_deg2 * ecliptic_density(v_depth) * latitude_profile(beta_deg)
}

/// Of those, the ones beyond 60 AU.
pub fn expected_distant(area_deg2: f64, beta_deg: f64, v_depth: f64) -> f64 {
    DISTANT_FRACTION * expected_tnos(area_deg2, beta_deg, v_depth)
}

/// Rubin's single-visit depth for a TNO-coloured source, in V.
pub const RUBIN_TNO_DEPTH_V: f64 = 24.5 + TNO_V_MINUS_R;
/// Fraction of its footprint Rubin images well enough to link movers.
pub const RUBIN_COVERAGE: f64 = 0.85;

/// Fraction of the TNOs brighter than `v_depth` in a tile that Rubin also
/// finds, given whether the tile is in its footprint and how much of the
/// tile survives Rubin's crowding.
pub fn rubin_share(v_depth: f64, in_footprint: bool, rubin_crowding: f64) -> f64 {
    if !in_footprint {
        return 0.0;
    }
    let reach = ecliptic_density(v_depth.min(RUBIN_TNO_DEPTH_V)) / ecliptic_density(v_depth);
    RUBIN_COVERAGE * rubin_crowding * reach
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rubin_takes_everything_it_can_reach() {
        assert_eq!(rubin_share(23.0, false, 1.0), 0.0);
        assert!((rubin_share(23.0, true, 1.0) - RUBIN_COVERAGE).abs() < 1e-12);
        assert!(rubin_share(26.5, true, 1.0) < 0.6 * RUBIN_COVERAGE);
        assert!(rubin_share(23.0, true, 0.1) < 0.1);
    }

    #[test]
    fn a_few_per_square_degree_at_rubin_depth() {
        let d = ecliptic_density(24.5);
        assert!((1.0..6.0).contains(&d), "{d} per deg2");
        assert!(ecliptic_density(26.0) > d);
        // The faint slope is shallower than the bright one.
        let bright = ecliptic_density(24.0) / ecliptic_density(23.0);
        let faint = ecliptic_density(27.0) / ecliptic_density(26.0);
        assert!(faint < bright);
    }

    #[test]
    fn the_belt_is_a_belt() {
        assert!((latitude_profile(0.0) - 1.0).abs() < 1e-12);
        assert!(latitude_profile(10.0) < 0.6);
        assert!(latitude_profile(40.0) < 0.03);
    }
}
