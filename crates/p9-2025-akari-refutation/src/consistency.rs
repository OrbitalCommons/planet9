//! The refutation: is the Phan et al. IRAS–AKARI candidate *pair* consistent
//! with a single bound Planet Nine?
//!
//! The Phan scheme links one IRAS source (mid-epoch 1983-06-24) to one AKARI
//! source (mid-epoch 2006-12-31), separated by `47.46′` on the sky over a
//! ~23.5 yr baseline (proper motion `≈ 2.06′/yr`), reused from
//! `p9_2025_iras_akari::orbital_constraints`.
//!
//! Chen et al. (§3.2) establish that the angular motion of a bound Planet
//! Nine decomposes into **parallax** and **proper motion**, and that the
//! parallax term sets the scale on which a genuine object is detectable.
//! Two independent budgets bound what a single bound object can do between
//! the two *real* epochs:
//!
//! 1. **Parallax budget.** With Earth's geometry fixed at the two real
//!    mid-epochs (a near-half-integer-year baseline, so Earth sits on nearly
//!    opposite sides of its orbit), the largest two-epoch shift a stationary
//!    distant object can show from parallax alone is the full parallactic
//!    diameter `2/d` rad — `≈ 13.7′` at 500 AU, `≈ 9.8′` at 700 AU. The
//!    observed 47.46′ is several times larger: parallax cannot account for
//!    the pair.
//!
//! 2. **Proper-motion budget.** The observed `2.06′/yr` exceeds the circular
//!    heliocentric transverse rate everywhere in the 500–700 AU band
//!    (`1.17–1.93′/yr`); matching it requires an eccentric orbit observed
//!    near perihelion *and* a constructively-aligned parallax. The implied
//!    circular-orbit distance falls *below* the 500 AU floor of the
//!    predicted P9 region.
//!
//! The bound-object two-epoch separation envelope at the real epochs (with
//! parallax orientation fixed by Earth and only the object's own orbit
//! scanned) is computed in [`bound_object_envelope_arcmin`]; the candidate
//! sits at its extreme and is reproducible only with a fine-tuned geometry.
//! The single headline inconsistency metric is the ratio of the observed
//! separation to the maximum bound-object parallax —
//! [`parallax_inconsistency_ratio`].

use nalgebra::Vector3;
use p9_2025_iras_akari::orbital_constraints::{
    derive_constraints, CandidatePair, TwoEpochConstraints,
};
use p9_2025_iras_akari::proper_motion::{
    annual_proper_motion_circular, implied_distance_circular, transverse_rate_rad_day,
};
use p9_2025_iras_akari::survey_model::{AkariFisSurvey, IrasSurvey};
use p9_core::constants::RAD2DEG;
use p9_core::coords::observer::{EarthProvider, EarthState, Time, Timescale};

/// Arcminutes per radian.
const ARCMIN_PER_RAD: f64 = RAD2DEG * 60.0;

/// The Phan et al. candidate's two-epoch constraints, reused verbatim from
/// the `p9-2025-iras-akari` reproduction (encoded Table 2 positions and
/// epochs). Separation `≈ 47.46′`, proper motion `≈ 2.06′/yr`, moving
/// south-west.
pub fn candidate_constraints() -> TwoEpochConstraints {
    derive_constraints(&CandidatePair::paper_candidate())
}

/// Earth heliocentric ecliptic positions (AU) at the two real survey
/// mid-epochs from `earth`. The near-half-integer-year baseline places Earth
/// on nearly opposite sides of its orbit at the two epochs.
pub fn earth_states_at_mid_epochs(earth: &mut impl EarthProvider) -> (EarthState, EarthState) {
    let ts = Timescale::default();
    let t1 = IrasSurvey::default().mid_epoch(&ts);
    let t2 = AkariFisSurvey::default().mid_epoch(&ts);
    (earth.earth_state(&t1), earth.earth_state(&t2))
}

/// The two real survey mid-epoch [`Time`]s (IRAS, AKARI).
pub fn mid_epochs() -> (Time, Time) {
    let ts = Timescale::default();
    (
        IrasSurvey::default().mid_epoch(&ts),
        AkariFisSurvey::default().mid_epoch(&ts),
    )
}

/// Maximum two-epoch angular separation (arcmin) a *stationary* distant
/// object at heliocentric distance `d_au` can show from **parallax alone**,
/// given Earth's positions `e1`, `e2` at the two real epochs. Scans the
/// object's ecliptic longitude on a fine grid; the maximum is the full
/// parallactic diameter `≈ 2/d` rad.
///
/// This is the parallax budget of Chen et al. §3.2 evaluated at the real
/// IRAS/AKARI geometry rather than for a generic six-month interval.
pub fn max_parallax_separation_arcmin(d_au: f64, e1: &Vector3<f64>, e2: &Vector3<f64>) -> f64 {
    const N: usize = 720;
    let mut max_sep = 0.0_f64;
    for k in 0..N {
        let lam = k as f64 * std::f64::consts::TAU / N as f64;
        // Object fixed in heliocentric ecliptic space (no orbital motion):
        // isolates the parallax contribution.
        let p = Vector3::new(d_au * lam.cos(), d_au * lam.sin(), 0.0);
        let g1 = p - e1;
        let g2 = p - e2;
        let sep = g1.angle(&g2) * ARCMIN_PER_RAD;
        max_sep = max_sep.max(sep);
    }
    max_sep
}

/// Two-epoch geocentric separation (arcmin) for a bound object at distance
/// `d_au` on orbit (`a_au`, `e`), starting at ecliptic longitude `lambda1`,
/// observed from Earth states `e1`/`e2` separated by `baseline_days`. Combines
/// parallax (Earth geometry) and proper motion (heliocentric transverse rate
/// reused from `p9-2025-iras-akari`) in the ecliptic plane.
fn bound_object_separation_arcmin(
    d_au: f64,
    a_au: f64,
    e: f64,
    lambda1: f64,
    e1: &Vector3<f64>,
    e2: &Vector3<f64>,
    baseline_days: f64,
) -> f64 {
    let mu = transverse_rate_rad_day(d_au, a_au, e); // rad/day
    let theta = mu * baseline_days;
    let p1 = Vector3::new(d_au * lambda1.cos(), d_au * lambda1.sin(), 0.0);
    let lam2 = lambda1 + theta;
    let p2 = Vector3::new(d_au * lam2.cos(), d_au * lam2.sin(), 0.0);
    let g1 = p1 - e1;
    let g2 = p2 - e2;
    g1.angle(&g2) * ARCMIN_PER_RAD
}

/// The bound-object two-epoch separation envelope (arcmin) at the real
/// IRAS/AKARI epochs: parallax orientation is *fixed* by Earth, and only the
/// object's own orbit is scanned — ecliptic longitude on a grid and the
/// circular / perihelion / aphelion branches of an `e = e_max` orbit. This
/// is the genuine "consistent with a bound object" window (unlike the
/// free-Earth-phase scan of `p9-2025-iras-akari`, which over-counts
/// parallax it cannot actually realise at fixed epochs).
///
/// Returns `(sep_min, sep_max)`.
pub fn bound_object_envelope_arcmin(
    d_au: f64,
    e_max: f64,
    e1: &Vector3<f64>,
    e2: &Vector3<f64>,
    baseline_days: f64,
) -> (f64, f64) {
    const N_LAMBDA: usize = 720;
    let branches = [
        (d_au, 0.0),
        (d_au / (1.0 - e_max), e_max),
        (d_au / (1.0 + e_max), e_max),
    ];
    let mut sep_min = f64::INFINITY;
    let mut sep_max = f64::NEG_INFINITY;
    for &(a, e) in &branches {
        for k in 0..N_LAMBDA {
            let lambda1 = k as f64 * std::f64::consts::TAU / N_LAMBDA as f64;
            let s = bound_object_separation_arcmin(d_au, a, e, lambda1, e1, e2, baseline_days);
            sep_min = sep_min.min(s);
            sep_max = sep_max.max(s);
        }
    }
    (sep_min, sep_max)
}

/// The headline inconsistency metric: the ratio of the observed candidate
/// separation to the maximum two-epoch parallax a bound object at `d_au`
/// could show at the real epochs. A ratio `≫ 1` means parallax — the
/// dominant term for a genuine bound P9 per Chen et al. — cannot account for
/// the candidate's motion.
pub fn parallax_inconsistency_ratio(d_au: f64, e1: &Vector3<f64>, e2: &Vector3<f64>) -> f64 {
    let observed = candidate_constraints().separation_arcmin;
    observed / max_parallax_separation_arcmin(d_au, e1, e2)
}

/// Baseline in days between the two real survey mid-epochs.
pub fn baseline_days() -> f64 {
    let (t1, t2) = mid_epochs();
    t2.tdb() - t1.tdb()
}

/// Result of the bound-object consistency check at one assumed distance.
#[derive(Debug, Clone)]
pub struct ConsistencyVerdict {
    /// Assumed heliocentric distance (AU).
    pub d_au: f64,
    /// Observed candidate separation (arcmin).
    pub observed_arcmin: f64,
    /// Maximum bound-object parallax at the real epochs (arcmin).
    pub max_parallax_arcmin: f64,
    /// Maximum bound-object total (parallax + proper motion) separation
    /// at the real epochs over the scanned orbits (arcmin).
    pub max_bound_total_arcmin: f64,
    /// Observed / max-parallax ratio (the headline inconsistency metric).
    pub parallax_excess_ratio: f64,
    /// Whether the observed separation lies inside the bound-object
    /// envelope at all (i.e. reproducible only at the envelope extreme).
    pub within_bound_envelope: bool,
}

/// Run the full consistency check at distance `d_au` using Earth states from
/// `earth` at the real survey mid-epochs.
pub fn check(earth: &mut impl EarthProvider, d_au: f64, e_max: f64) -> ConsistencyVerdict {
    let (s1, s2) = earth_states_at_mid_epochs(earth);
    let (e1, e2) = (s1.0, s2.0);
    let bdays = baseline_days();
    let observed = candidate_constraints().separation_arcmin;
    let max_par = max_parallax_separation_arcmin(d_au, &e1, &e2);
    let (env_min, env_max) = bound_object_envelope_arcmin(d_au, e_max, &e1, &e2, bdays);
    ConsistencyVerdict {
        d_au,
        observed_arcmin: observed,
        max_parallax_arcmin: max_par,
        max_bound_total_arcmin: env_max,
        parallax_excess_ratio: observed / max_par,
        within_bound_envelope: observed >= env_min && observed <= env_max,
    }
}

/// Implied heliocentric distance (AU) from the candidate's proper motion
/// under the circular-orbit, parallax-free inversion — the *lower bound* on
/// distance. Reuses the `p9-2025-iras-akari` shared inversion.
pub fn implied_circular_distance() -> f64 {
    implied_distance_circular(candidate_constraints().proper_motion_arcmin_yr)
}

/// Maximum circular-orbit proper motion (arcmin/yr) anywhere in the
/// `[d_min, d_max]` band: the slowest-allowed bound-object rate floor is at
/// `d_max`, the fastest circular rate at `d_min`.
pub fn max_circular_proper_motion(d_min: f64, _d_max: f64) -> f64 {
    annual_proper_motion_circular(d_min)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::published;
    use p9_core::constants::YEAR_DAYS;
    use p9_core::coords::observer::CircularEarth;

    /// Earth states from the ephemeris when available, else the analytic
    /// circular-Earth fallback (the parallax budget is set by Earth's
    /// orbital geometry, which the circular model reproduces to ~0.03 AU —
    /// well inside the arcmin tolerances pinned here).
    fn earth_states() -> (Vector3<f64>, Vector3<f64>) {
        use p9_core::coords::observer::EphemerisEarth;
        if let Ok(mut eph) = EphemerisEarth::try_new() {
            let (s1, s2) = earth_states_at_mid_epochs(&mut eph);
            (s1.0, s2.0)
        } else {
            let mut circ = CircularEarth;
            let (s1, s2) = earth_states_at_mid_epochs(&mut circ);
            (s1.0, s2.0)
        }
    }

    #[test]
    fn candidate_separation_and_pm_match_iras_akari_crate() {
        let c = candidate_constraints();
        assert!(
            (c.separation_arcmin - 47.46).abs() < 0.5,
            "separation = {:.2}'",
            c.separation_arcmin
        );
        assert!(
            (c.proper_motion_arcmin_yr - 2.06).abs() < 0.1,
            "pm = {:.2}'/yr",
            c.proper_motion_arcmin_yr
        );
    }

    #[test]
    fn real_baseline_is_about_23_5_years() {
        let yr = baseline_days() / YEAR_DAYS;
        assert!((yr - 23.52).abs() < 0.05, "baseline = {yr:.2} yr");
    }

    #[test]
    fn max_parallax_is_the_diameter_two_over_d() {
        // The parallax-only maximum at the real epochs should be close to
        // the full parallactic diameter 2/d rad (Earth on nearly opposite
        // sides over the ~23.5 yr near-half-integer baseline).
        let (e1, e2) = earth_states();
        for &d in &[500.0, 600.0, 700.0] {
            let max_par = max_parallax_separation_arcmin(d, &e1, &e2);
            let diameter = 2.0 / d * ARCMIN_PER_RAD;
            assert!(
                (max_par - diameter).abs() < 0.6,
                "d={d}: max parallax {max_par:.2}' vs diameter {diameter:.2}'"
            );
        }
    }

    #[test]
    fn observed_separation_far_exceeds_parallax_budget() {
        // Headline refutation: 47.46' is several times larger than the
        // largest parallax a bound object could show at any distance in the
        // search band. Parallax — the dominant term for a real bound P9 —
        // cannot explain the pair.
        let (e1, e2) = earth_states();
        for &d in &[
            published::DISTANCE_MIN_AU,
            600.0,
            published::DISTANCE_MAX_AU,
        ] {
            let ratio = parallax_inconsistency_ratio(d, &e1, &e2);
            assert!(
                ratio > 3.0,
                "d={d}: observed/parallax = {ratio:.2}, expected >3 (inconsistent)"
            );
        }
        // At 500 AU (most generous, largest parallax) the excess is ~3.5x.
        let ratio_500 = parallax_inconsistency_ratio(published::DISTANCE_MIN_AU, &e1, &e2);
        assert!(
            (ratio_500 - 3.45).abs() < 0.3,
            "500 AU parallax excess = {ratio_500:.2}, expected ~3.45"
        );
    }

    #[test]
    fn candidate_proper_motion_exceeds_circular_bound_in_band() {
        // The observed 2.06'/yr is faster than a circular orbit anywhere in
        // the 500-700 AU band: matching it needs an eccentric orbit observed
        // near perihelion (a special configuration), not a generic bound P9.
        let pm = candidate_constraints().proper_motion_arcmin_yr;
        let max_circ =
            max_circular_proper_motion(published::DISTANCE_MIN_AU, published::DISTANCE_MAX_AU);
        assert!(
            pm > max_circ,
            "candidate pm {pm:.2}'/yr should exceed circular max {max_circ:.2}'/yr"
        );
    }

    #[test]
    fn implied_circular_distance_is_below_search_floor() {
        // Interpreting the candidate as a circular-orbit bound object puts it
        // at ~478 AU — below the 500 AU floor of the predicted P9 region.
        let d = implied_circular_distance();
        assert!(
            d < published::DISTANCE_MIN_AU,
            "implied circular distance = {d:.0} AU, expected < 500"
        );
        assert!((d - 478.0).abs() < 20.0, "implied distance = {d:.0} AU");
    }

    #[test]
    fn magnitude_alone_is_not_a_clean_discriminator_at_500au() {
        // Honest caveat: when the object's own ecliptic longitude is free,
        // the *magnitude* of the bound-object envelope at 500 AU spans the
        // observed 47.46' (its circular upper edge ~59'). This is exactly why
        // the magnitude looked "consistent" to Phan et al. — and why the
        // refutation rests on the parallax budget and direction, not the
        // separation magnitude alone.
        let (e1, e2) = earth_states();
        let bdays = baseline_days();
        let observed = candidate_constraints().separation_arcmin;
        let (env_min, env_max) = bound_object_envelope_arcmin(500.0, 0.7, &e1, &e2, bdays);
        assert!(
            observed > env_min && observed < env_max,
            "observed {observed:.1}' lies inside 500 AU magnitude envelope [{env_min:.1}, {env_max:.1}]"
        );
    }

    #[test]
    fn at_700au_bound_envelope_cannot_reach_observed() {
        // At the far edge of the search band a bound object simply cannot
        // produce 47.46' at the real epochs, even at e_max perihelion.
        let (e1, e2) = earth_states();
        let bdays = baseline_days();
        let observed = candidate_constraints().separation_arcmin;
        let (_, max_700) = bound_object_envelope_arcmin(700.0, 0.7, &e1, &e2, bdays);
        assert!(
            max_700 < observed,
            "700 AU eccentric max {max_700:.1}' should be < observed {observed:.1}'"
        );
    }

    #[test]
    fn check_verdict_flags_inconsistency() {
        let mut circ = CircularEarth;
        // The parallax budget is exceeded at every distance in the band.
        let v600 = check(&mut circ, 600.0, 0.7);
        assert!(v600.parallax_excess_ratio > 3.0);
        assert!(v600.observed_arcmin > v600.max_parallax_arcmin);

        // At 700 AU the bound-object total (parallax + proper motion) cannot
        // reach the observed separation even at e_max perihelion, so the
        // candidate falls outside the envelope entirely.
        let v700 = check(&mut circ, 700.0, 0.7);
        assert!(
            v700.max_bound_total_arcmin < v700.observed_arcmin,
            "700 AU bound total {:.1}' vs observed {:.1}'",
            v700.max_bound_total_arcmin,
            v700.observed_arcmin
        );
        assert!(!v700.within_bound_envelope);
    }
}
