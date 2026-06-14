//! W1 detectability of a Planet Nine and the maximum heliocentric distance at
//! which a planet of a given mass crosses the WISE W1 depth.
//!
//! The detectability question reduces to: at what heliocentric distance does
//! the modelled W1 magnitude (thermal + reflected, `thermal_model`) equal the
//! survey depth? Because W1 magnitude increases (faints) monotonically with
//! distance, there is a single finite cross-over distance for any mass whose
//! near-Sun magnitude is brighter than the depth — the maximum detectable
//! distance. It is found by bisection.

use crate::survey_model::WiseSurvey;
use crate::thermal_model::{w1_magnitude, P9Thermal};

/// Whether a P9 of `mass_earth` at `distance_au` is brighter than the WISE W1
/// depth (i.e. its modelled W1 magnitude is at or below `survey.w1_depth`).
pub fn is_detectable_w1(survey: &WiseSurvey, mass_earth: f64, distance_au: f64) -> bool {
    w1_magnitude(mass_earth, distance_au) <= survey.w1_depth
}

/// Maximum heliocentric distance (AU) at which a P9 of `mass_earth` reaches
/// the WISE W1 depth, found by bisection on the monotone magnitude–distance
/// relation. Returns `None` if even the nearest distance considered
/// (`d_min`) is already fainter than the depth (the planet is undetectable at
/// any distance in the search range).
pub fn max_detectable_distance(survey: &WiseSurvey, mass_earth: f64) -> Option<f64> {
    let d_min = 50.0_f64;
    let d_max = 50_000.0_f64;

    // Brightness decreases monotonically with distance: detectable at d_min,
    // not at d_max, with a single crossing in between.
    if !is_detectable_w1(survey, mass_earth, d_min) {
        return None;
    }
    if is_detectable_w1(survey, mass_earth, d_max) {
        // Detectable even at the far edge of the search range; report the edge.
        return Some(d_max);
    }

    let mut lo = d_min;
    let mut hi = d_max;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if is_detectable_w1(survey, mass_earth, mid) {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-4 {
            break;
        }
    }
    Some(0.5 * (lo + hi))
}

/// Per-position W1 detection probability: the magnitude-efficiency at the
/// planet's modelled W1 magnitude, gated by the galactic-plane footprint mask.
pub fn detection_probability(survey: &WiseSurvey, p9: &P9Thermal, galactic_lat_deg: f64) -> f64 {
    let footprint = if survey.in_footprint(galactic_lat_deg) {
        1.0
    } else {
        0.0
    };
    survey.detection_efficiency(p9.w1_magnitude()) * footprint
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nearby_massive_planet_is_detectable() {
        let survey = WiseSurvey::default();
        // A 15 M⊕ planet at 180 AU is brighter than the W1 depth in reflected
        // light (the cold-body W1 signal; see thermal_model docs).
        assert!(is_detectable_w1(&survey, 15.0, 180.0));
    }

    #[test]
    fn distant_low_mass_planet_is_not_detectable() {
        let survey = WiseSurvey::default();
        // A 5 M⊕ planet at 1000 AU is far below the W1 depth.
        assert!(!is_detectable_w1(&survey, 5.0, 1000.0));
    }

    #[test]
    fn max_distance_is_finite_and_decreases_with_mass() {
        let survey = WiseSurvey::default();
        let d_heavy = max_detectable_distance(&survey, 15.0).expect("heavy detectable");
        let d_light = max_detectable_distance(&survey, 8.0).expect("light detectable");
        assert!(d_heavy.is_finite() && d_light.is_finite());
        // A more massive (larger, brighter) planet is detectable farther out.
        assert!(
            d_heavy > d_light,
            "heavier planet reaches farther: {d_heavy:.0} vs {d_light:.0} AU"
        );
        // The WISE W1 cold-body exclusion reaches the low hundreds of AU.
        assert!(
            (100.0..400.0).contains(&d_heavy),
            "d_max(15 M⊕) = {d_heavy:.0} AU"
        );
    }

    #[test]
    fn crossover_magnitude_equals_depth() {
        let survey = WiseSurvey::default();
        let d = max_detectable_distance(&survey, 12.0).expect("detectable");
        let m = w1_magnitude(12.0, d);
        assert!((m - survey.w1_depth).abs() < 1e-2, "W1({d:.1} AU) = {m:.3}");
    }

    #[test]
    fn footprint_gates_detection_probability() {
        let survey = WiseSurvey::default();
        let p9 = P9Thermal::new(15.0, 170.0);
        // Out of the galactic plane: high probability; in the plane: zero.
        assert!(detection_probability(&survey, &p9, 45.0) > 0.5);
        assert_eq!(detection_probability(&survey, &p9, 2.0), 0.0);
    }
}
