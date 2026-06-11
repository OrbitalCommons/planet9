//! Apparent sky motion of Planet Nine between the IRAS (1983) and AKARI
//! (2006) epochs: heliocentric orbital motion plus geocentric parallax.
//!
//! Heliocentric transverse (sky-plane) angular rate of a Kepler orbit:
//!
//!   µ = h / r² = √(GM·a(1−e²)) / r²   [rad/day]
//!
//! (the angular-momentum law — NOT the total vis-viva speed divided by r,
//! which wrongly counts the radial velocity as sky motion). For a circular
//! orbit this reduces to √(GM/d³): ~1.9'/yr at 500 AU, ~1.2'/yr at 700 AU.
//!
//! The apparent positions catalogued by the two surveys are *geocentric*:
//! Earth's orbital position contributes a parallax of ±1/d rad (±6.9' at
//! 500 AU) at each epoch, up to ~13.8' between epochs. Heliocentric motion
//! alone — even parabolic, √2 × circular = 62.8' over 23 yr at 500 AU — can
//! never reach the 69.6' upper bound of Phan et al. (2025); the parallax
//! term is what closes the gap.
//!
//! Model (documented simplifications): the object moves in the ecliptic
//! plane at its heliocentric transverse rate (radial drift over 23 yr is
//! ≲1% of r for a ~15,000 yr orbit and is neglected); Earth is on a
//! circular 1 AU orbit; the survey epochs sample arbitrary Earth orbital
//! phases (both surveys scanned for ~1 yr). Eccentricities are scanned up
//! to `E_MAX` = 0.7, the approximate upper bound of published Planet Nine
//! orbit fits (Batygin et al. 2019; Brown & Batygin 2021).

use p9_core::constants::{GM_SUN, RAD2DEG, YEAR_DAYS};

/// Upper bound on Planet Nine eccentricity scanned by the separation-window
/// model (Batygin et al. 2019 / Brown & Batygin 2021 posteriors give
/// e ≲ 0.6–0.7).
pub const E_MAX: f64 = 0.7;

/// Arcminutes per radian.
const ARCMIN_PER_RAD: f64 = RAD2DEG * 60.0;

/// Circular orbital velocity at heliocentric distance `d_au` in AU/day.
pub fn circular_velocity(d_au: f64) -> f64 {
    (GM_SUN / d_au).sqrt()
}

/// Total speed at heliocentric distance `r_au` for an orbit with semi-major
/// axis `a_au`, from the vis-viva equation: v² = GM(2/r − 1/a). This is the
/// *speed*, not the sky-plane rate — use `transverse_rate_arcmin_yr` for
/// proper motion.
pub fn vis_viva_velocity(r_au: f64, a_au: f64) -> f64 {
    (GM_SUN * (2.0 / r_au - 1.0 / a_au)).sqrt()
}

/// Orbital period in years for an orbit with semi-major axis `a_au`.
pub fn orbital_period_years(a_au: f64) -> f64 {
    let period_days = 2.0 * std::f64::consts::PI * (a_au.powi(3) / GM_SUN).sqrt();
    period_days / YEAR_DAYS
}

/// Heliocentric transverse angular rate µ = h/r² = √(GM·a(1−e²))/r² in
/// rad/day, for an object at heliocentric distance `r_au` on an orbit
/// (`a_au`, `e`).
pub fn transverse_rate_rad_day(r_au: f64, a_au: f64, e: f64) -> f64 {
    (GM_SUN * a_au * (1.0 - e * e)).sqrt() / (r_au * r_au)
}

/// Heliocentric transverse angular rate in arcmin/yr.
pub fn transverse_rate_arcmin_yr(r_au: f64, a_au: f64, e: f64) -> f64 {
    transverse_rate_rad_day(r_au, a_au, e) * ARCMIN_PER_RAD * YEAR_DAYS
}

/// Annual proper motion for a circular orbit at distance `d_au` (arcmin/yr).
pub fn annual_proper_motion_circular(d_au: f64) -> f64 {
    transverse_rate_arcmin_yr(d_au, d_au, 0.0)
}

/// Heliocentric-only angular separation accumulated over `baseline_years`
/// (arcmin) — no parallax. Lower-level building block kept for comparison
/// against the geocentric model.
pub fn heliocentric_separation_arcmin(r_au: f64, a_au: f64, e: f64, baseline_years: f64) -> f64 {
    transverse_rate_arcmin_yr(r_au, a_au, e) * baseline_years
}

/// Apparent *geocentric* angular separation (arcmin) between the two survey
/// epochs for an object at heliocentric distance `d_au` on an orbit
/// (`a_au`, `e`), with Earth at orbital phases `earth_phase1`/`earth_phase2`
/// (radians) at the two epochs.
///
/// Geometry (coplanar, documented in the module header): the object starts
/// at (d, 0), sweeps the heliocentric angle µ·Δt, and Earth sits on the
/// unit circle at each epoch. The separation is the angle between the two
/// geocentric direction vectors.
pub fn apparent_separation_arcmin(
    d_au: f64,
    a_au: f64,
    e: f64,
    baseline_years: f64,
    earth_phase1: f64,
    earth_phase2: f64,
) -> f64 {
    let theta = transverse_rate_rad_day(d_au, a_au, e) * baseline_years * YEAR_DAYS;

    // Object heliocentric positions at the two epochs.
    let (p1x, p1y) = (d_au, 0.0);
    let (p2x, p2y) = (d_au * theta.cos(), d_au * theta.sin());

    // Earth heliocentric positions (circular 1 AU orbit).
    let (e1x, e1y) = (earth_phase1.cos(), earth_phase1.sin());
    let (e2x, e2y) = (earth_phase2.cos(), earth_phase2.sin());

    // Geocentric direction vectors.
    let (g1x, g1y) = (p1x - e1x, p1y - e1y);
    let (g2x, g2y) = (p2x - e2x, p2y - e2y);

    let dot = g1x * g2x + g1y * g2y;
    let n1 = (g1x * g1x + g1y * g1y).sqrt();
    let n2 = (g2x * g2x + g2y * g2y).sqrt();
    let angle = (dot / (n1 * n2)).clamp(-1.0, 1.0).acos();
    angle * ARCMIN_PER_RAD
}

/// Range of apparent two-epoch separations (arcmin) achievable at
/// heliocentric distance `d_au`, scanning eccentricity branches
/// (circular; perihelion and aphelion of an e = `e_max` orbit) and all
/// Earth phase pairs on a grid.
///
/// Returns `(sep_min, sep_max)`.
pub fn separation_window(d_au: f64, e_max: f64, baseline_years: f64) -> (f64, f64) {
    const N_PHASE: usize = 48;
    let mut sep_min = f64::INFINITY;
    let mut sep_max = f64::NEG_INFINITY;

    // (a, e) branches: circular, observed at perihelion of an eccentric
    // orbit (fastest sky motion), observed at aphelion (slowest).
    let branches = [
        (d_au, 0.0),
        (d_au / (1.0 - e_max), e_max),
        (d_au / (1.0 + e_max), e_max),
    ];

    for &(a, e) in &branches {
        for i in 0..N_PHASE {
            let phi1 = i as f64 * std::f64::consts::TAU / N_PHASE as f64;
            for j in 0..N_PHASE {
                let phi2 = j as f64 * std::f64::consts::TAU / N_PHASE as f64;
                let s = apparent_separation_arcmin(d_au, a, e, baseline_years, phi1, phi2);
                sep_min = sep_min.min(s);
                sep_max = sep_max.max(s);
            }
        }
    }
    (sep_min, sep_max)
}

/// Implied heliocentric distance for an observed two-epoch separation:
/// the *largest* distance at which the separation is achievable within the
/// model (perihelion of an e ≤ `E_MAX` orbit, most favorable parallax) —
/// the convention of Phan et al. (2025), whose 42'–69.6' window maps to
/// 500–700 AU. Solved by bisection on the monotone-decreasing maximum of
/// `separation_window`.
pub fn implied_distance(sep_arcmin: f64, baseline_years: f64) -> f64 {
    let max_sep = |d: f64| separation_window(d, E_MAX, baseline_years).1;
    let (mut lo, mut hi) = (50.0_f64, 5000.0_f64);
    if max_sep(hi) >= sep_arcmin {
        return hi;
    }
    if max_sep(lo) <= sep_arcmin {
        return lo;
    }
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if max_sep(mid) >= sep_arcmin {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Implied heliocentric distance from a circular-orbit, parallax-free
/// proper-motion inversion: d = (√GM / µ)^{2/3} with µ in rad/day. This is
/// the *lower-bound* convention (heliocentric circular motion only); the
/// geocentric model in `implied_distance` is the paper's convention.
///
/// Single shared implementation (previously triplicated across this crate).
pub fn implied_distance_circular(proper_motion_arcmin_yr: f64) -> f64 {
    let mu_rad_day = proper_motion_arcmin_yr / (ARCMIN_PER_RAD * YEAR_DAYS);
    (GM_SUN.sqrt() / mu_rad_day).powf(2.0 / 3.0)
}

/// Distance range (AU) implied by a separation window over a baseline,
/// using the geocentric maximum-motion inversion of `implied_distance`:
/// larger separations map to smaller distances.
///
/// Returns `(d_min, d_max)` for `(sep_max, sep_min)`. With the paper's
/// 42'–69.6' window and 23 yr this recovers ≈ (510, 740) AU, consistent
/// with the paper's 500–700 AU search range.
pub fn distance_range_for_separation(
    sep_min_arcmin: f64,
    sep_max_arcmin: f64,
    baseline_years: f64,
) -> (f64, f64) {
    (
        implied_distance(sep_max_arcmin, baseline_years),
        implied_distance(sep_min_arcmin, baseline_years),
    )
}

/// Proper motion data for a grid of distances.
///
/// Returns Vec of (distance_au, pm_circular, pm_eccentric, sep_circular,
/// sep_eccentric) where the eccentric case puts the object at the
/// perihelion of an orbit with semi-major axis `a_au` (e = 1 − r/a),
/// using the heliocentric transverse rate.
pub fn proper_motion_table(
    distances: &[f64],
    a_au: f64,
    baseline_years: f64,
) -> Vec<(f64, f64, f64, f64, f64)> {
    distances
        .iter()
        .map(|&r| {
            let mu_circ = annual_proper_motion_circular(r);
            let e = (1.0 - r / a_au).max(0.0);
            let mu_ecc = transverse_rate_arcmin_yr(r, a_au, e);
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
        assert!(
            (mu - 1.93).abs() < 0.1,
            "µ_circ at 500 AU = {mu:.2}'/yr, expected ~1.93"
        );
    }

    #[test]
    fn circular_proper_motion_at_700au() {
        let mu = annual_proper_motion_circular(700.0);
        assert!(
            (mu - 1.17).abs() < 0.1,
            "µ_circ at 700 AU = {mu:.2}'/yr, expected ~1.17"
        );
    }

    #[test]
    fn transverse_rate_at_perihelion_scales_sqrt_one_plus_e() {
        // At perihelion r = a(1−e): µ/µ_circ(r) = √(a(1−e²)/r) = √(1+e).
        let (r, e) = (500.0, 0.5);
        let a = r / (1.0 - e);
        let ratio = transverse_rate_arcmin_yr(r, a, e) / annual_proper_motion_circular(r);
        assert!(
            (ratio - 1.5_f64.sqrt()).abs() < 1e-9,
            "ratio = {ratio:.4}, expected √1.5"
        );
    }

    #[test]
    fn transverse_rate_below_total_speed_rate_off_apsis() {
        // The old (wrong) model used v_total/r; off the apsides the radial
        // component makes that strictly larger than the true sky rate h/r².
        let (r, a, e) = (500.0, 600.0, 0.5);
        let mu_transverse = transverse_rate_rad_day(r, a, e);
        let mu_total_over_r = vis_viva_velocity(r, a) / r;
        assert!(
            mu_transverse < mu_total_over_r,
            "transverse {mu_transverse:.3e} should be < v/r {mu_total_over_r:.3e}"
        );
    }

    #[test]
    fn heliocentric_motion_alone_cannot_reach_69_6_arcmin() {
        // Even a parabolic flyby (e → 1, perihelion branch) gives at most
        // √2 × circular = 62.8' at 500 AU over 23 yr: the 69.6' bound of the
        // paper requires the parallax term.
        let e = 0.999;
        let sep = heliocentric_separation_arcmin(500.0, 500.0 / (1.0 - e), e, 23.0);
        assert!(sep < 69.6, "heliocentric-only sep = {sep:.1}'");
        assert!((sep - 62.8).abs() < 1.0, "parabolic limit = {sep:.1}'");
    }

    #[test]
    fn parallax_widens_the_window_at_500au() {
        // For circular orbits at 500 AU the geocentric window spans
        // ~±13.8' (two ±6.9' epochs) around the 44.4' heliocentric drift.
        let (lo, hi) = separation_window(500.0, 0.0, 23.0);
        assert!(hi - lo > 20.0, "window [{lo:.1}, {hi:.1}] too narrow");
        let sep_helio = heliocentric_separation_arcmin(500.0, 500.0, 0.0, 23.0);
        assert!(lo < sep_helio && hi > sep_helio);
    }

    #[test]
    fn separation_window_accommodates_paper_range() {
        // Phan et al. window: 42'–69.6' for 500–700 AU over 23 yr.
        let (_, max_500) = separation_window(500.0, E_MAX, 23.0);
        assert!(max_500 >= 69.6, "max sep at 500 AU = {max_500:.1}'");
        let (min_700, max_700) = separation_window(700.0, E_MAX, 23.0);
        assert!(
            min_700 <= 42.0 && max_700 >= 42.0,
            "42' not within [{min_700:.1}, {max_700:.1}] at 700 AU"
        );
    }

    #[test]
    fn implied_distance_recovers_paper_range() {
        // Inverting the paper's window ends recovers ~500–700 AU.
        let d_at_696 = implied_distance(69.6, 23.0);
        assert!(
            d_at_696 > 480.0 && d_at_696 < 560.0,
            "d(69.6') = {d_at_696:.0} AU, expected ~510"
        );
        let d_at_42 = implied_distance(42.0, 23.0);
        assert!(
            d_at_42 > 680.0 && d_at_42 < 790.0,
            "d(42') = {d_at_42:.0} AU, expected ~740"
        );
    }

    #[test]
    fn distance_range_for_paper_window() {
        let (d_min, d_max) = distance_range_for_separation(42.0, 69.6, 23.0);
        assert!(d_min < d_max);
        assert!(d_min > 480.0 && d_min < 560.0, "d_min = {d_min:.0}");
        assert!(d_max > 680.0 && d_max < 790.0, "d_max = {d_max:.0}");
    }

    #[test]
    fn circular_separation_23yr_at_500au() {
        let sep = heliocentric_separation_arcmin(500.0, 500.0, 0.0, 23.0);
        assert!(
            (sep - 44.4).abs() < 2.0,
            "sep at 500 AU (circular) over 23yr = {sep:.1}'"
        );
    }

    #[test]
    fn circular_separation_23yr_at_700au() {
        let sep = heliocentric_separation_arcmin(700.0, 700.0, 0.0, 23.0);
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
    fn implied_distance_circular_round_trip() {
        let d = 500.0;
        let mu = annual_proper_motion_circular(d);
        let d_back = implied_distance_circular(mu);
        assert!((d_back - d).abs() < 1.0, "d_back = {d_back:.1}");
    }

    #[test]
    fn apparent_separation_equals_heliocentric_when_parallax_cancels() {
        // With Earth at the same phase at both epochs the parallax nearly
        // cancels (small secondary term from the object's sky motion).
        let sep_geo = apparent_separation_arcmin(500.0, 500.0, 0.0, 23.0, 0.3, 0.3);
        let sep_helio = heliocentric_separation_arcmin(500.0, 500.0, 0.0, 23.0);
        assert!(
            (sep_geo - sep_helio).abs() < 1.0,
            "geo {sep_geo:.2}' vs helio {sep_helio:.2}'"
        );
    }

    #[test]
    fn orbital_period_reasonable() {
        let p = orbital_period_years(600.0);
        assert!((p - 14_697.0).abs() < 100.0, "period at 600 AU = {p:.0} yr");
    }

    #[test]
    fn vis_viva_matches_circular_when_r_equals_a() {
        let v_circ = circular_velocity(600.0);
        let v_vv = vis_viva_velocity(600.0, 600.0);
        assert!((v_circ - v_vv).abs() < 1e-12);
    }
}
