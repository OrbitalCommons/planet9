//! TNO colours (paper Section 5, after Pfalzner, Wagner & Gibbon 2025).
//!
//! The protoplanetary disc is postulated to carry a radial colour gradient:
//! very red bodies in the cold-belt region, progressively more neutral
//! surfaces further out, modelled as a linear spectral-slope gradient from
//! 30 AU to 150 AU (their Figure 8a). Because the flyby delivers the
//! high-inclination and high-eccentricity orbits from the *outer* disc, the
//! model predicts a scarcity of very red objects at i > 21° (Marsset et al.
//! 2019) and e > 0.42 (Ali-Dib et al. 2021) — and the 2026 paper finds that
//! this survives 4.5 Gyr of evolution.

use p9_core::constants::RAD2DEG;
use p9_core::types::OrbitalElements;

/// Inner and outer ends of the colour gradient (AU).
pub const GRADIENT_RANGE_AU: (f64, f64) = (30.0, 150.0);
/// Spectral slope (%/100 nm) at the inner end: very red.
pub const SLOPE_INNER: f64 = 40.0;
/// Spectral slope at the outer end: neutral/grey.
pub const SLOPE_OUTER: f64 = 5.0;
/// "Very red" threshold on the spectral slope.
pub const VERY_RED_SLOPE: f64 = 25.0;
/// Inclination boundary of the observed colour dichotomy (deg).
pub const I_SPLIT_DEG: f64 = 21.0;
/// Eccentricity boundary of the observed colour dichotomy.
pub const E_SPLIT: f64 = 0.42;

/// Provisional spectral slope of a body formed at heliocentric radius
/// `r0_au`: linear between the gradient ends, clamped outside.
pub fn spectral_slope(r0_au: f64) -> f64 {
    let (r_in, r_out) = GRADIENT_RANGE_AU;
    let t = ((r0_au - r_in) / (r_out - r_in)).clamp(0.0, 1.0);
    SLOPE_INNER + t * (SLOPE_OUTER - SLOPE_INNER)
}

pub fn is_very_red(r0_au: f64) -> bool {
    spectral_slope(r0_au) > VERY_RED_SLOPE
}

/// Fraction of very red bodies among the `(r0, elements)` pairs passing
/// `pred`; `None` if none pass.
pub fn very_red_fraction<'a>(
    bodies: impl IntoIterator<Item = (f64, &'a OrbitalElements)>,
    pred: impl Fn(&OrbitalElements) -> bool,
) -> Option<f64> {
    let (mut n, mut red) = (0usize, 0usize);
    for (r0, el) in bodies {
        if pred(el) {
            n += 1;
            if is_very_red(r0) {
                red += 1;
            }
        }
    }
    (n > 0).then(|| red as f64 / n as f64)
}

/// The paper's Figure 8b/c bins: `(low-i, high-i, low-e, high-e)` very-red
/// fractions over bodies with 5° < i < 50°.
pub fn colour_dichotomy<'a>(
    bodies: impl IntoIterator<Item = (f64, &'a OrbitalElements)> + Clone,
) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
    let i_deg = |el: &OrbitalElements| el.i * RAD2DEG;
    (
        very_red_fraction(bodies.clone(), |el| (5.0..I_SPLIT_DEG).contains(&i_deg(el))),
        very_red_fraction(bodies.clone(), |el| {
            (I_SPLIT_DEG..50.0).contains(&i_deg(el))
        }),
        very_red_fraction(bodies.clone(), |el| el.e < E_SPLIT),
        very_red_fraction(bodies, |el| el.e >= E_SPLIT),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gradient_is_monotone_red_to_grey() {
        assert!(is_very_red(35.0));
        assert!(!is_very_red(140.0));
        assert!(spectral_slope(30.0) > spectral_slope(90.0));
        assert!(spectral_slope(90.0) > spectral_slope(150.0));
        assert_eq!(spectral_slope(10.0), SLOPE_INNER);
        assert_eq!(spectral_slope(500.0), SLOPE_OUTER);
    }
}
