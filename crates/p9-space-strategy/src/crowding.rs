//! Stellar crowding: how much of a field is lost under stars, for any
//! telescope, as a function of where it points in the Galaxy.
//!
//! A moving-object search loses the sky within roughly one seeing disc of
//! every star it cannot subtract cleanly. With `N` stars per unit area
//! brighter than the survey depth and a masked area `A` per star, the
//! surviving fraction of a field is `exp(−N·A)` (Poisson voids). The same
//! expression, with each telescope's own resolution, describes why the
//! ground surveys go blind in the Galactic plane and why a sharp space
//! telescope does not.
//!
//! Star counts come from the Bahcall & Soneira (1980, ApJS 44, 73; Appendix
//! B) fitting formula for the integrated V-band counts of their disc +
//! spheroid Galaxy model. The formula has no interstellar extinction and is
//! stated for |b| ≳ 20°; it is continued to the plane here with the latitude
//! floored at [`B_FLOOR_DEG`] and the disc projection factor floored at
//! [`PROJECTION_FLOOR`], which keeps it finite toward the inner Galaxy. At
//! low latitude it therefore over-counts (dust hides stars), so the
//! completeness it implies is conservative for every telescope.

use p9_core::constants::DEG2RAD;

/// Latitude floor for the continued formula (deg).
pub const B_FLOOR_DEG: f64 = 2.0;
/// Floor on sin b·(1 − μ cot b cos l), the disc projection factor.
pub const PROJECTION_FLOOR: f64 = 0.02;

/// Integrated star counts N(< V) per deg² toward Galactic (l, b) in degrees.
pub fn stars_per_deg2(v_limit: f64, gal_l_deg: f64, gal_b_deg: f64) -> f64 {
    let m = v_limit;
    let b = gal_b_deg.abs().max(B_FLOOR_DEG) * DEG2RAD;
    let l = gal_l_deg * DEG2RAD;

    // Disc component.
    let (c1, alpha, beta, delta, m_star) = (925.0, -0.132, 0.035, 3.0, 15.75);
    let mu = if m <= 12.0 {
        0.03
    } else if m <= 20.0 {
        0.0075 * (m - 12.0) + 0.03
    } else {
        0.09
    };
    let gamma = if m <= 12.0 {
        0.36
    } else if m <= 20.0 {
        0.04 * (12.0 - m) + 0.36
    } else {
        0.04
    };
    let projection = (b.sin() * (1.0 - mu * l.cos() / b.tan())).max(PROJECTION_FLOOR);
    let disc = c1 * 10f64.powf(beta * (m - m_star))
        / (1.0 + 10f64.powf(alpha * (m - m_star))).powf(delta)
        / projection.powf(3.0 - 5.0 * gamma);

    // Spheroid component.
    let (c2, kappa, eta, lambda, m_dagger) = (1050.0, -0.180, 0.087, 2.50, 17.5);
    let sigma = 1.45 - 0.20 * b.cos() * l.cos();
    let spheroid = c2 * 10f64.powf(eta * (m - m_dagger))
        / (1.0 + 10f64.powf(kappa * (m - m_dagger))).powf(lambda)
        / (1.0 - b.cos() * l.cos()).max(PROJECTION_FLOOR).powf(sigma);

    disc + spheroid
}

/// Fraction of a field that survives crowding for a telescope that loses a
/// disc of radius `mask_radius_arcsec` around every star brighter than
/// `v_limit`.
pub fn completeness(v_limit: f64, mask_radius_arcsec: f64, gal_l_deg: f64, gal_b_deg: f64) -> f64 {
    let n = stars_per_deg2(v_limit, gal_l_deg, gal_b_deg);
    let mask_deg2 = std::f64::consts::PI * (mask_radius_arcsec / 3600.0).powi(2);
    (-n * mask_deg2).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_match_the_high_latitude_anchors() {
        // Galactic pole to V = 24: a few thousand stars per deg² (deep HST
        // fields find ~1–3 stars per arcmin² to V ~ 25).
        let pole = stars_per_deg2(24.0, 0.0, 90.0);
        assert!((3_000.0..8_000.0).contains(&pole), "pole = {pole}");
        // The classic V < 21 pole count is ~2,000 per deg².
        let pole21 = stars_per_deg2(21.0, 0.0, 90.0);
        assert!((1_200.0..3_500.0).contains(&pole21), "pole(21) = {pole21}");
    }

    #[test]
    fn counts_rise_toward_the_plane_and_the_centre() {
        let high = stars_per_deg2(24.0, 180.0, 40.0);
        let mid = stars_per_deg2(24.0, 180.0, 10.0);
        let plane = stars_per_deg2(24.0, 180.0, 0.0);
        assert!(high < mid && mid < plane);
        assert!(stars_per_deg2(24.0, 0.0, 10.0) > 3.0 * mid);
        assert!(plane.is_finite() && stars_per_deg2(24.0, 0.0, 0.0).is_finite());
    }

    #[test]
    fn a_sharp_telescope_keeps_the_anticentre_plane() {
        // 0.6" mask (space) vs 3" mask (2" ground seeing) at the anticentre.
        let space = completeness(24.5, 0.6, 180.0, 0.0);
        let ground = completeness(20.5, 3.0, 180.0, 0.0);
        assert!(
            space > 0.75,
            "space completeness at the anticentre = {space}"
        );
        assert!(ground < space);
        // Toward the inner Galaxy nothing survives from the ground.
        assert!(completeness(21.5, 2.0, 10.0, 0.0) < 0.2);
    }
}
