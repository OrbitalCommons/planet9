//! Proper motion and angular separation calculations for Planet Nine.
//!
//! For a circular orbit at heliocentric distance d (AU):
//!   µ_circular ≈ sqrt(GM_sun) / d^(3/2)  [rad/day]
//!
//! giving ~1.9'/yr at 500 AU and ~1.2'/yr at 700 AU.
//!
//! For eccentric orbits, the velocity at distance r with semi-major axis a
//! follows the vis-viva equation:
//!   v = sqrt(GM * (2/r - 1/a))
//!
//! This can be significantly higher than circular velocity when r is near
//! perihelion. Phan et al. (2025) use an angular separation window of
//! 42'–69.6' for the 500–700 AU distance range over 23 years, which
//! accounts for the full range of orbital eccentricities consistent with
//! P9 constraints.

use p9_core::constants::{GM_SUN, RAD2DEG, YEAR_DAYS};

/// Circular orbital velocity at heliocentric distance `d_au` in AU/day.
pub fn circular_velocity(d_au: f64) -> f64 {
    (GM_SUN / d_au).sqrt()
}

/// Velocity at heliocentric distance `r_au` for an orbit with semi-major
/// axis `a_au`, using the vis-viva equation: v² = GM(2/r - 1/a).
pub fn vis_viva_velocity(r_au: f64, a_au: f64) -> f64 {
    (GM_SUN * (2.0 / r_au - 1.0 / a_au)).sqrt()
}

/// Orbital period in years for an orbit with semi-major axis `a_au`.
pub fn orbital_period_years(a_au: f64) -> f64 {
    let period_days = 2.0 * std::f64::consts::PI * (a_au.powi(3) / GM_SUN).sqrt();
    period_days / YEAR_DAYS
}

/// Annual proper motion in arcminutes/year for an object at heliocentric
/// distance `r_au` with orbital velocity determined by the vis-viva equation
/// for semi-major axis `a_au`.
///
/// For a circular orbit, use `a_au = r_au`.
pub fn annual_proper_motion_arcmin(r_au: f64, a_au: f64) -> f64 {
    let v = vis_viva_velocity(r_au, a_au);
    let omega_rad_per_day = v / r_au;
    omega_rad_per_day * RAD2DEG * 60.0 * YEAR_DAYS
}

/// Annual proper motion for a circular orbit at distance `d_au`.
pub fn annual_proper_motion_circular(d_au: f64) -> f64 {
    annual_proper_motion_arcmin(d_au, d_au)
}

/// Angular separation in arcminutes between two survey epochs
/// separated by `baseline_years`, for an object at distance `r_au`
/// with semi-major axis `a_au`.
pub fn angular_separation_arcmin(r_au: f64, a_au: f64, baseline_years: f64) -> f64 {
    annual_proper_motion_arcmin(r_au, a_au) * baseline_years
}

/// Compute the distance range (AU) that corresponds to a given angular
/// separation range over a given baseline for CIRCULAR orbits, by
/// inverting the proper motion relation.
///
/// Returns (d_min, d_max) in AU.
pub fn distance_range_for_separation(
    sep_min_arcmin: f64,
    sep_max_arcmin: f64,
    baseline_years: f64,
) -> (f64, f64) {
    let k = GM_SUN.sqrt() * RAD2DEG * 60.0 * YEAR_DAYS * baseline_years;
    let d_min = (k / sep_max_arcmin).powf(2.0 / 3.0);
    let d_max = (k / sep_min_arcmin).powf(2.0 / 3.0);
    (d_min, d_max)
}

/// Proper motion data for a grid of distances.
///
/// Returns Vec of (distance_au, pm_circular, pm_eccentric, sep_circular, sep_eccentric)
/// where the eccentric case uses the Batygin & Brown (2016) nominal a=700 AU.
pub fn proper_motion_table(
    distances: &[f64],
    a_au: f64,
    baseline_years: f64,
) -> Vec<(f64, f64, f64, f64, f64)> {
    distances
        .iter()
        .map(|&r| {
            let mu_circ = annual_proper_motion_circular(r);
            let mu_ecc = annual_proper_motion_arcmin(r, a_au);
            (
                r,
                mu_circ,
                mu_ecc,
                mu_circ * baseline_years,
                mu_ecc * baseline_years,
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circular_proper_motion_at_500au() {
        let mu = annual_proper_motion_circular(500.0);
        // Circular orbit at 500 AU gives ~1.93'/yr
        assert!(
            (mu - 1.93).abs() < 0.1,
            "µ_circ at 500 AU = {mu:.2}'/yr, expected ~1.93"
        );
    }

    #[test]
    fn circular_proper_motion_at_700au() {
        let mu = annual_proper_motion_circular(700.0);
        // Circular orbit at 700 AU gives ~1.17'/yr
        assert!(
            (mu - 1.17).abs() < 0.1,
            "µ_circ at 700 AU = {mu:.2}'/yr, expected ~1.17"
        );
    }

    #[test]
    fn eccentric_orbit_higher_velocity_at_perihelion() {
        // An object at r=500 AU on an orbit with a=1500 AU, e=2/3
        // (perihelion = 500 AU) moves faster than circular
        let mu_circ = annual_proper_motion_circular(500.0);
        let mu_ecc = annual_proper_motion_arcmin(500.0, 1500.0);
        assert!(
            mu_ecc > mu_circ,
            "eccentric µ={mu_ecc:.2}' should exceed circular µ={mu_circ:.2}'"
        );
        // vis-viva ratio: sqrt(2/r - 1/a) / sqrt(1/r) = sqrt(2 - r/a)
        // = sqrt(2 - 500/1500) = sqrt(2 - 1/3) = sqrt(5/3) = 1.29
        let expected_ratio = (5.0_f64 / 3.0).sqrt();
        let actual_ratio = mu_ecc / mu_circ;
        assert!(
            (actual_ratio - expected_ratio).abs() < 0.01,
            "ratio = {actual_ratio:.3}, expected {expected_ratio:.3}"
        );
    }

    #[test]
    fn paper_separation_window_with_eccentric_orbits() {
        // Phan et al. use 42'–69.6' for 500–700 AU over 23 years.
        // The upper bound (69.6' at 500 AU → 3.03'/yr) requires
        // v/v_circ ≈ 1.57, or a ≈ d/(2 - 1.57²) ≈ negative → needs
        // very eccentric orbit or includes parallax contribution.
        //
        // With a=700 AU, e=0.6 (nominal P9): at r=500 AU,
        // v/v_circ = sqrt(2 - 500/700) = sqrt(1.286) = 1.13
        // → µ ≈ 1.93 * 1.13 = 2.19'/yr → sep ≈ 50' over 23yr
        let sep_ecc = angular_separation_arcmin(500.0, 700.0, 23.0);
        assert!(
            sep_ecc > 44.0 && sep_ecc < 55.0,
            "sep = {sep_ecc:.1}', expected ~50'"
        );
    }

    #[test]
    fn circular_separation_23yr_at_500au() {
        let sep = angular_separation_arcmin(500.0, 500.0, 23.0);
        assert!(
            (sep - 44.4).abs() < 2.0,
            "sep at 500 AU (circular) over 23yr = {sep:.1}'"
        );
    }

    #[test]
    fn circular_separation_23yr_at_700au() {
        let sep = angular_separation_arcmin(700.0, 700.0, 23.0);
        assert!(
            (sep - 26.8).abs() < 2.0,
            "sep at 700 AU (circular) over 23yr = {sep:.1}'"
        );
    }

    #[test]
    fn proper_motion_decreases_with_distance() {
        let mu_near = annual_proper_motion_circular(500.0);
        let mu_far = annual_proper_motion_circular(700.0);
        assert!(mu_near > mu_far);
    }

    #[test]
    fn distance_range_round_trip() {
        let baseline = 23.0;
        let sep_500 = angular_separation_arcmin(500.0, 500.0, baseline);
        let sep_700 = angular_separation_arcmin(700.0, 700.0, baseline);
        let (d_min, d_max) = distance_range_for_separation(sep_700, sep_500, baseline);
        assert!(
            (d_min - 500.0).abs() < 1.0,
            "d_min = {d_min:.1}, expected 500"
        );
        assert!(
            (d_max - 700.0).abs() < 1.0,
            "d_max = {d_max:.1}, expected 700"
        );
    }

    #[test]
    fn orbital_period_reasonable() {
        let p = orbital_period_years(600.0);
        assert!((p - 14_697.0).abs() < 100.0, "period at 600 AU = {p:.0} yr");
    }

    #[test]
    fn paper_search_window_distances() {
        // The paper's 42'–69.6' window for circular orbits maps to
        // ~371–525 AU. The paper's stated 500–700 AU range implies
        // their window accounts for eccentric orbits with higher velocities.
        let baseline = 23.0;
        let (d_min, d_max) = distance_range_for_separation(42.0, 69.6, baseline);
        assert!(d_min > 300.0 && d_min < 450.0, "d_min = {d_min:.0} AU");
        assert!(d_max > 450.0 && d_max < 600.0, "d_max = {d_max:.0} AU");
    }

    #[test]
    fn vis_viva_matches_circular_when_r_equals_a() {
        let v_circ = circular_velocity(600.0);
        let v_vv = vis_viva_velocity(600.0, 600.0);
        assert!(
            (v_circ - v_vv).abs() < 1e-12,
            "circular and vis-viva should match: {v_circ} vs {v_vv}"
        );
    }
}
