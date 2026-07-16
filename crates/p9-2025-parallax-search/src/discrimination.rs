//! Separating a bound distant planet from the stellar background.
//!
//! Two numbers make the parallax method work:
//!
//! 1. **Contrast.** A 500–700 AU body has a reflex parallax of order arcminutes
//!    ([`crate::parallax`]); a background reference star at hundreds of pc has a
//!    parallax of a few milliarcsec. The ratio is `10³–10⁴`, so P9 moves between
//!    epochs by a margin no field star can fake.
//!
//! 2. **Required astrometric precision.** To claim a shift at signal-to-noise
//!    `snr`, each epoch's position must be measured to better than
//!    `amplitude / snr`. Because the amplitude falls as `1/d`, the *tolerable*
//!    per-epoch error also falls as `1/d`: detecting P9 farther out demands ever
//!    finer astrometry. We report it both as the maximum tolerable error (arcsec,
//!    falling with `d`) and as the demanded precision in pixels for a given plate
//!    scale.

use crate::parallax::{epoch_parallax_arcsec, half_parallax_arcsec};
use p9_core::units::{arcseconds, au, Angle, Length};

/// Contrast ratio between a Planet Nine at `distance_au` and a background star of
/// parallax `star_parallax_arcsec`, using the **half-angle** (annual) parallax of
/// each. P9's reflex parallax is `1 AU / d`; the star's is its true parallax.
/// The ratio is dimensionless and `≫ 1`.
pub fn contrast_ratio(distance_au: f64, star_parallax_arcsec: f64) -> f64 {
    half_parallax_arcsec(distance_au) / star_parallax_arcsec
}

/// Maximum tolerable per-epoch astrometric error (arcsec) to detect the reflex
/// parallax of a body at `distance_au` at the requested `snr`, over a baseline of
/// `epoch_separation_days`. A measured shift `A` detected at `snr` needs
/// per-epoch error `σ ≤ A / snr`. Falls as `1/d`.
pub fn max_tolerable_error_arcsec(distance_au: f64, snr: f64, epoch_separation_days: f64) -> f64 {
    epoch_parallax_arcsec(distance_au, epoch_separation_days) / snr
}

/// Maximum tolerable per-epoch astrometric error as a dimension-checked [`Angle`].
pub fn max_tolerable_error(distance_au: f64, snr: f64, epoch_separation_days: f64) -> Angle {
    arcseconds(max_tolerable_error_arcsec(
        distance_au,
        snr,
        epoch_separation_days,
    ))
}

/// The same demanded precision expressed in detector pixels for a given plate
/// scale (arcsec/pixel). Smaller is harder; the value falls as `1/d`.
pub fn required_precision_pixels(
    distance_au: f64,
    snr: f64,
    epoch_separation_days: f64,
    plate_scale_arcsec_per_px: f64,
) -> f64 {
    max_tolerable_error_arcsec(distance_au, snr, epoch_separation_days) / plate_scale_arcsec_per_px
}

/// Largest heliocentric distance (AU) at which a body's reflex parallax is still
/// detectable at `snr` given a per-epoch astrometric error `sigma_arcsec` over a
/// baseline of `epoch_separation_days`. Solve `A(d) / snr = sigma` for `d` using
/// `A(d) = baseline_arcsec / d`.
pub fn max_detectable_distance_au(sigma_arcsec: f64, snr: f64, epoch_separation_days: f64) -> f64 {
    // A(d) = epoch_parallax_arcsec(d) = K / d, so K = A(d0) * d0 for any d0.
    let d0 = 500.0;
    let k = epoch_parallax_arcsec(d0, epoch_separation_days) * d0;
    k / (snr * sigma_arcsec)
}

/// Largest detectable heliocentric distance as a dimension-checked [`Length`].
pub fn max_detectable_distance(sigma_arcsec: f64, snr: f64, epoch_separation_days: f64) -> Length {
    au(max_detectable_distance_au(
        sigma_arcsec,
        snr,
        epoch_separation_days,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parallax::YEAR_DAYS;
    use crate::published;
    use approx::assert_relative_eq;
    use uom::si::angle::second as arcsecond;
    use uom::si::length::astronomical_unit;

    #[test]
    fn typed_siblings_match_f64_sources() {
        let base = YEAR_DAYS / 2.0;
        assert_relative_eq!(
            max_tolerable_error(700.0, 5.0, base).get::<arcsecond>(),
            max_tolerable_error_arcsec(700.0, 5.0, base),
            max_relative = 1e-12
        );
        assert_relative_eq!(
            max_detectable_distance(1.0, 5.0, base).get::<astronomical_unit>(),
            max_detectable_distance_au(1.0, 5.0, base),
            max_relative = 1e-12
        );
    }

    #[test]
    fn contrast_against_field_star_is_thousands() {
        // 460 AU P9 vs a 1 mas field star: ~7.5' / 1 mas.
        let r = contrast_ratio(
            published::REFERENCE_DISTANCE_AU,
            published::TYPICAL_FIELD_STAR_PARALLAX_ARCSEC,
        );
        assert!(r > 1.0e3, "contrast = {r:.0}, expected >1e3");
        // 460 AU half-parallax 448" / 0.001" = ~4.5e5.
        assert_relative_eq!(r, 4.484e5, epsilon = 1.0e3);
    }

    #[test]
    fn contrast_exceeds_1e3_to_1e4_across_prior_and_realistic_stars() {
        // Even against an unusually nearby (10 pc, 0.1") reference star, a P9 at
        // the far edge of the prior still beats it by >10^3.
        let nearby_star = 0.1; // 10 pc, parallax 0.1 arcsec
        let r = contrast_ratio(published::P9_DISTANCE_MAX_AU, nearby_star);
        assert!(r > 1.0e3, "even vs a 10 pc star, contrast = {r:.0}");
    }

    #[test]
    fn required_precision_grows_with_distance() {
        // Tolerable per-epoch error shrinks as 1/d: detecting a farther P9 needs
        // finer astrometry.
        let snr = 5.0;
        let base = YEAR_DAYS / 2.0;
        let e500 = max_tolerable_error_arcsec(500.0, snr, base);
        let e700 = max_tolerable_error_arcsec(700.0, snr, base);
        assert!(e700 < e500, "precision must tighten with distance");
        // Exactly inverse-linear in d.
        assert_relative_eq!(e500 / e700, 700.0 / 500.0, epsilon = 1e-9);
    }

    #[test]
    fn required_precision_is_loose_for_p9_full_baseline() {
        // Over a 6-month baseline at SNR 5, a 700 AU P9 still allows ~1' of
        // per-epoch error -- far coarser than routine sub-arcsec astrometry.
        let sigma = max_tolerable_error_arcsec(700.0, 5.0, YEAR_DAYS / 2.0);
        assert!(sigma > 60.0, "tolerable error = {sigma:.1} arcsec");
    }

    #[test]
    fn max_detectable_distance_scales_inversely_with_precision() {
        // Halving the achievable error doubles the reach.
        let snr = 5.0;
        let base = YEAR_DAYS / 2.0;
        let d_coarse = max_detectable_distance_au(1.0, snr, base);
        let d_fine = max_detectable_distance_au(0.5, snr, base);
        assert_relative_eq!(d_fine / d_coarse, 2.0, epsilon = 1e-9);
        // With 1" per-epoch error and SNR 5 over 6 months, reach is well beyond
        // the P9 prior (thousands of AU), confirming depth (not astrometry) is
        // the real limiter -- consistent with the paper's magnitude-limited
        // conclusion.
        assert!(d_coarse > 700.0);
    }

    #[test]
    fn required_precision_pixels_at_jbt_scale() {
        // At the 0.130"/px JBT/SPENCER scale, a 700 AU P9 over 6 months at SNR 5
        // tolerates hundreds of pixels of per-epoch error -- trivially met.
        let px = required_precision_pixels(700.0, 5.0, YEAR_DAYS / 2.0, 0.130);
        assert!(px > 100.0, "tolerable error = {px:.0} px");
    }
}
