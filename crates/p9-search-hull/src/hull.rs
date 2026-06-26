//! The parameter-space reach hull: the (heliocentric distance × apparent V)
//! plane, with the boundary that the wide surveys have already searched.
//!
//! A reflected-light planet of fixed size gets fainter with distance as
//! `V(d) = H + 5 log10(d (d-1))` (opposition), so a survey of depth `m`
//! reaches out to the distance where `V(d) = m`. The deepest all-sky depth
//! therefore defines a horizontal cut in this plane: everything brighter
//! (nearer) than that line has been searched; the posterior cloud that sits
//! *below* the line (fainter / farther) is what remains. All magnitudes come
//! from [`p9_core::analysis::photometry`]; the inversion reuses the monotone
//! bisection in [`p9_core::analysis::thermal::max_detectable_distance`].

use p9_core::analysis::photometry::{planet_apparent_magnitude, ALBEDO_NEPTUNE};
use p9_core::analysis::thermal::max_detectable_distance;

/// Distance bracket for the magnitude curve and all reach inversions (AU).
pub const D_MIN: f64 = 80.0;
pub const D_MAX: f64 = 2_000.0;

/// Fiducial Planet Nine mass for the parameter-space panel (Earth masses) —
/// the Brown & Batygin 2021 MCMC median.
pub const FIDUCIAL_MASS_EARTH: f64 = 6.2;

/// Apparent reflected V at distance `d_au` for the fiducial planet.
pub fn fiducial_v(d_au: f64) -> f64 {
    planet_apparent_magnitude(FIDUCIAL_MASS_EARTH, ALBEDO_NEPTUNE, d_au)
}

/// Max heliocentric distance (AU) at which a planet of `mass_earth` would be
/// brighter than `depth`; `None` if it is undetectable across the bracket.
pub fn reach_distance(mass_earth: f64, depth: f64) -> Option<f64> {
    max_detectable_distance(D_MIN, D_MAX, depth, |d| {
        planet_apparent_magnitude(mass_earth, ALBEDO_NEPTUNE, d)
    })
}

/// The V(d) curve for the fiducial planet, sampled on a log-spaced distance
/// grid of `n` points across the bracket.
pub fn magnitude_curve(n: usize) -> (Vec<f64>, Vec<f64>) {
    let lo = D_MIN.ln();
    let hi = D_MAX.ln();
    let mut ds = Vec::with_capacity(n);
    let mut vs = Vec::with_capacity(n);
    for k in 0..n {
        let f = k as f64 / (n - 1) as f64;
        let d = (lo + f * (hi - lo)).exp();
        ds.push(d);
        vs.push(fiducial_v(d));
    }
    (ds, vs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn brightness_falls_with_distance() {
        assert!(fiducial_v(300.0) < fiducial_v(700.0));
    }

    #[test]
    fn deeper_surveys_reach_farther() {
        let shallow = reach_distance(FIDUCIAL_MASS_EARTH, 20.5).unwrap();
        let deep = reach_distance(FIDUCIAL_MASS_EARTH, 24.5).unwrap();
        assert!(deep > shallow, "deep {deep} !> shallow {shallow}");
    }

    #[test]
    fn reach_is_consistent_with_the_curve() {
        // V at the reach distance should equal the depth we inverted for.
        let depth = 22.0;
        let d = reach_distance(FIDUCIAL_MASS_EARTH, depth).unwrap();
        assert!((fiducial_v(d) - depth).abs() < 0.05);
    }

    #[test]
    fn curve_is_monotone_and_spans_the_bracket() {
        let (ds, vs) = magnitude_curve(50);
        assert_eq!(ds.len(), 50);
        assert!((ds[0] - D_MIN).abs() < 1e-6);
        assert!((ds[49] - D_MAX).abs() < 1e-6);
        for w in vs.windows(2) {
            assert!(w[1] >= w[0]);
        }
    }
}
