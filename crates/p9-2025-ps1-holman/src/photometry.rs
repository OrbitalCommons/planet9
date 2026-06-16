//! Reflected-light photometry for the Holman PS1 search.
//!
//! Reuses `p9_core::analysis::photometry` entirely (Neptune-anchored
//! mass-radius relation, reflected-sunlight distance law, Neptune albedo
//! 0.41). The crate-specific quantity is the *maximum heliocentric
//! distance* at which a P9 of a given mass is still brighter than the
//! shift-and-stack effective depth — i.e. the outer edge of the search's
//! sensitivity, which grows with assumed mass.

use p9_core::analysis::photometry::{planet_apparent_magnitude, ALBEDO_NEPTUNE};
use p9_core::units::{au, Length};

use crate::survey_model::Ps1StackSurvey;

/// Apparent r magnitude (Neptune albedo, at opposition) of a P9 of
/// `mass_earth` at heliocentric distance `r_au`. Thin wrapper over the
/// shared `p9-core` reflected-sunlight law.
pub fn apparent_r_magnitude(mass_earth: f64, r_au: f64) -> f64 {
    planet_apparent_magnitude(mass_earth, ALBEDO_NEPTUNE, r_au)
}

/// Maximum heliocentric distance (AU) at which a P9 of `mass_earth` is still
/// at or brighter than the survey's effective (shift-and-stack) depth, i.e.
/// `apparent_r_magnitude(mass, r) == survey.effective_depth()`.
///
/// Apparent magnitude rises monotonically with distance (reflected light
/// dims as r^4 at opposition: m ~ H + 5 log10(r(r-1))), so the equation has
/// a unique root, found by bisection on [`R_MIN_AU`, `R_MAX_AU`]. Returns
/// `R_MAX_AU` if even there the object is detectable, `R_MIN_AU` if it is
/// never detectable.
pub fn max_detectable_distance(survey: &Ps1StackSurvey, mass_earth: f64) -> f64 {
    let target = survey.effective_depth();
    let mut lo = R_MIN_AU;
    let mut hi = R_MAX_AU;
    // At lo the object is brightest (smallest m); detectable means m <= target.
    if apparent_r_magnitude(mass_earth, lo) > target {
        return R_MIN_AU; // not detectable anywhere in range
    }
    if apparent_r_magnitude(mass_earth, hi) <= target {
        return R_MAX_AU; // detectable everywhere in range
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if apparent_r_magnitude(mass_earth, mid) <= target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// [`max_detectable_distance`] as a dimension-checked [`Length`].
pub fn max_detectable_distance_typed(survey: &Ps1StackSurvey, mass_earth: f64) -> Length {
    au(max_detectable_distance(survey, mass_earth))
}

/// Inner bound of the bisection bracket (AU). Objects closer than the giant
/// planets are outside the P9 search regime.
pub const R_MIN_AU: f64 = 30.0;

/// Outer bound of the bisection bracket (AU): well beyond any plausible P9
/// aphelion.
pub const R_MAX_AU: f64 = 5000.0;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn magnitude_matches_core_law() {
        // The wrapper must be exactly the shared p9-core law (no local copy).
        assert_relative_eq!(
            apparent_r_magnitude(6.0, 500.0),
            planet_apparent_magnitude(6.0, ALBEDO_NEPTUNE, 500.0),
            epsilon = 1e-12
        );
    }

    #[test]
    fn max_distance_is_finite_and_at_threshold() {
        let s = Ps1StackSurvey::default();
        let r = max_detectable_distance(&s, 6.0);
        // PINNED: a finite AU value strictly inside the bracket.
        assert!(r.is_finite());
        assert!(r > R_MIN_AU && r < R_MAX_AU, "r_max = {r:.1}");
        // At the root the apparent magnitude equals the effective depth.
        assert!(
            (apparent_r_magnitude(6.0, r) - s.effective_depth()).abs() < 1e-3,
            "m(r_max) = {:.3} vs depth {:.3}",
            apparent_r_magnitude(6.0, r),
            s.effective_depth()
        );
    }

    #[test]
    fn typed_distance_matches_f64() {
        use uom::si::length::astronomical_unit;
        let s = Ps1StackSurvey::default();
        assert_relative_eq!(
            max_detectable_distance_typed(&s, 6.0).get::<astronomical_unit>(),
            max_detectable_distance(&s, 6.0),
            epsilon = 1e-12
        );
    }

    #[test]
    fn max_distance_increases_with_mass() {
        let s = Ps1StackSurvey::default();
        let r5 = max_detectable_distance(&s, 5.0);
        let r10 = max_detectable_distance(&s, 10.0);
        let r20 = max_detectable_distance(&s, 20.0);
        // PINNED: more massive (larger, brighter) P9 detectable farther out.
        assert!(r10 > r5, "r10 {r10:.1} should exceed r5 {r5:.1}");
        assert!(r20 > r10, "r20 {r20:.1} should exceed r10 {r10:.1}");
    }

    #[test]
    fn deeper_stack_reaches_farther() {
        // The shift-and-stack gain extends the reach versus single-epoch.
        let stack = Ps1StackSurvey::default();
        let single = Ps1StackSurvey {
            stack_depth_gain: 0.0,
            ..Ps1StackSurvey::default()
        };
        let r_stack = max_detectable_distance(&stack, 6.0);
        let r_single = max_detectable_distance(&single, 6.0);
        assert!(
            r_stack > r_single,
            "stack reach {r_stack:.1} should exceed single-epoch {r_single:.1}"
        );
    }
}
