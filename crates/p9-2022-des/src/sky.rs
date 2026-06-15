//! Orbit -> apparent sky-position mapping for the DES footprint model.
//!
//! The footprint of `survey_model` is defined in *equatorial* coordinates
//! (RA/dec), so orbital sky positions must be rotated from the heliocentric
//! ecliptic frame to the equatorial frame before any footprint test (the
//! analogous fix to p9-2021-ztf's `sky.rs`; substituting ecliptic latitude
//! for declination misclassifies exactly the survey-boundary regions that
//! drive unique-coverage estimates). The rotation itself comes from
//! p9-core's `coords::sky` (starfield's ECLIPJ2000 frame matrix).
//!
//! Apparent (geocentric) positions are computed per observing night: the
//! geocentric vector is the heliocentric vector minus Earth's orbital
//! position, which captures both the annual parallactic oscillation
//! (~1/r radians, i.e. +/-0.14 deg at 400 AU) and the slow orbital drift
//! of the object across the 6-year survey baseline.
//!
//! Two Earth models are available:
//! - `apparent_position_with_earth_deg` (primary): a real Earth state from
//!   `p9_core::coords::observer` (DE421 via `EphemerisEarth`), with
//!   light-time retardation and stellar aberration applied through the
//!   shared apparent chain. Per-night Earth states are precomputed once
//!   per survey grid (`night_earth_states`) so the kernel is not hit in
//!   the per-orbit loops.
//! - `apparent_position_deg` (analytic fallback): a circular, coplanar
//!   1 AU Earth with phase zero at the survey start, accurate to ~0.03 deg
//!   in parallax — kept for the speed-critical inner MC loops.

use nalgebra::Vector3;
use p9_core::constants::{GM_SUN, YEAR_DAYS};
use p9_core::coords::observer::{
    apparent_direction_from, EarthProvider, EarthState, Time, Timescale,
};
use p9_core::coords::sky::ecliptic_vec_to_equatorial_deg;
use p9_core::types::OrbitalElements;

/// Start of the DES wide survey (Y1 began 2013-08-31; the six observing
/// seasons run through early 2019). Anchor for the night grid's day
/// offsets.
pub fn survey_start(ts: &Timescale) -> Time {
    ts.utc((2013, 8, 31))
}

/// Apparent geocentric equatorial (RA, dec) in degrees of an orbit
/// `t_days` after its stored epoch, plus the heliocentric distance (AU).
///
/// The object's mean anomaly is advanced by its mean motion; Earth is
/// modeled on a circular coplanar 1-AU orbit (phase zero at t = 0), which is
/// accurate to ~0.03 deg in parallax — ample for footprint membership.
pub fn apparent_position_deg(elem: &OrbitalElements, t_days: f64) -> (f64, f64, f64) {
    let pos = heliocentric_position(elem, t_days);
    let r_helio = pos.norm();

    let theta = std::f64::consts::TAU * t_days / YEAR_DAYS;
    let (sin_t, cos_t) = theta.sin_cos();
    let geo = pos - Vector3::new(cos_t, sin_t, 0.0);

    let (ra, dec) = ecliptic_vec_to_equatorial_deg(&geo);
    (ra, dec, r_helio)
}

/// Heliocentric ecliptic position (AU) of an orbit `t_days` after its
/// stored epoch (mean anomaly advanced at the mean motion).
fn heliocentric_position(elem: &OrbitalElements, t_days: f64) -> Vector3<f64> {
    let n = (GM_SUN / (elem.a * elem.a * elem.a)).sqrt(); // rad/day
    let mut advanced = *elem;
    advanced.mean_anomaly = (elem.mean_anomaly + n * t_days).rem_euclid(std::f64::consts::TAU);
    advanced.to_state_vector(GM_SUN).pos
}

/// Earth states (heliocentric ecliptic J2000) at `offsets_days` after the
/// anchor epoch `start`. Precompute these once per night grid; the grid is
/// shared by every orbit in a recovery run.
pub fn earth_states_at(
    earth: &mut impl EarthProvider,
    ts: &Timescale,
    start: &Time,
    offsets_days: &[f64],
) -> Vec<EarthState> {
    offsets_days
        .iter()
        .map(|&dt| earth.earth_state(&ts.tt_jd(start.tt() + dt, None)))
        .collect()
}

/// Ephemeris-grade apparent geocentric equatorial (RA, dec) in degrees of
/// an orbit `t_days` after its stored epoch, observed from the given Earth
/// state (from `earth_states_at`), plus the heliocentric distance (AU).
///
/// Unlike the analytic `apparent_position_deg`, this uses a real Earth
/// position and the full apparent chain: light-time retardation (~2-4 days
/// at P9 distances) and stellar aberration.
pub fn apparent_position_with_earth_deg(
    elem: &OrbitalElements,
    earth_state: &EarthState,
    t_days: f64,
) -> (f64, f64, f64) {
    let r_helio = heliocentric_position(elem, t_days).norm();
    let geo = apparent_direction_from(earth_state, |dt| heliocentric_position(elem, t_days + dt));
    let (ra, dec) = ecliptic_vec_to_equatorial_deg(&geo);
    (ra, dec, r_helio)
}

/// Sun–object–Earth phase angle (radians) of an orbit `t_days` after its
/// stored epoch, observed from the given Earth state (instantaneous
/// geometry — the few-day light-time displacement changes α by far less
/// than the H-G law resolves at these distances). Bounded by
/// asin(1 AU / r): below 0.35° everywhere in the BB21 reference
/// population (minimum r ≈ 176 AU), ~0.11° at the fiducial 500 AU.
pub fn phase_angle_with_earth(
    elem: &OrbitalElements,
    earth_state: &EarthState,
    t_days: f64,
) -> f64 {
    let helio = heliocentric_position(elem, t_days);
    let delta = (helio - earth_state.0).norm();
    p9_core::analysis::photometry::phase_angle(helio.norm(), delta, earth_state.0.norm())
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    #[test]
    fn test_ecliptic_pole_maps_to_obliquity_complement() {
        // The north ecliptic pole sits at dec = 90 - 23.4393 = 66.56 deg.
        let (_, dec) = ecliptic_vec_to_equatorial_deg(&Vector3::new(0.0, 0.0, 1.0));
        assert!((dec - (90.0 - 23.439_291)).abs() < 1e-6, "dec = {dec}");
    }

    #[test]
    fn test_equinox_direction_unchanged() {
        // The vernal-equinox direction is shared by both frames.
        let (ra, dec) = ecliptic_vec_to_equatorial_deg(&Vector3::new(1.0, 0.0, 0.0));
        assert!(ra.abs() < 1e-9 && dec.abs() < 1e-9);
    }

    #[test]
    fn test_ecliptic_plane_declination_bounded_by_obliquity() {
        // Positions in the ecliptic plane reach at most |dec| = obliquity.
        for k in 0..36 {
            let lambda = k as f64 * 10.0 * DEG2RAD;
            let (_, dec) =
                ecliptic_vec_to_equatorial_deg(&Vector3::new(lambda.cos(), lambda.sin(), 0.0));
            assert!(dec.abs() <= 23.44 + 1e-9, "dec = {dec}");
        }
    }

    #[test]
    fn ephemeris_apparent_position_close_to_analytic() {
        use p9_core::coords::observer::EphemerisEarth;
        let mut earth = match EphemerisEarth::try_new() {
            Ok(e) => e,
            Err(_) => {
                eprintln!("skipping: no ephemeris kernel in starfield cache");
                return;
            }
        };
        let ts = Timescale::default();
        let start = survey_start(&ts);
        let elem = OrbitalElements {
            a: 400.0,
            e: 0.0,
            i: 0.3,
            omega: 0.0,
            omega_big: 1.0,
            mean_anomaly: 2.0,
        };
        // The analytic model's Earth phase is arbitrary (zero at t = 0),
        // so the two paths can differ by up to the full parallax
        // amplitude, 2·asin(1/400) ≈ 0.29 deg, plus aberration (~0.006
        // deg) — but no more.
        let offsets = [0.0, 200.0, 800.0, 1500.0, 2100.0];
        let states = earth_states_at(&mut earth, &ts, &start, &offsets);
        for (&t_days, state) in offsets.iter().zip(&states) {
            let (ra_a, dec_a, r_a) = apparent_position_deg(&elem, t_days);
            let (ra_e, dec_e, r_e) = apparent_position_with_earth_deg(&elem, state, t_days);
            assert!((r_a - r_e).abs() < 1e-9);
            let sep = p9_core::coords::sky::angular_distance(
                &Vector3::new(
                    dec_a.to_radians().cos() * ra_a.to_radians().cos(),
                    dec_a.to_radians().cos() * ra_a.to_radians().sin(),
                    dec_a.to_radians().sin(),
                ),
                &Vector3::new(
                    dec_e.to_radians().cos() * ra_e.to_radians().cos(),
                    dec_e.to_radians().cos() * ra_e.to_radians().sin(),
                    dec_e.to_radians().sin(),
                ),
            )
            .to_degrees();
            assert!(
                sep < 0.30,
                "ephemeris vs analytic apparent position differs by {sep:.4} deg at t = {t_days}"
            );
        }
    }

    /// Issue #41: measure the phase-induced Δmag over the reference
    /// population at the real DES observing geometry. The phase angle is
    /// bounded by asin(1 AU/r), tiny at P9 distances, but the H-G law's
    /// opposition surge (dΔm/dα → ∞ at α = 0) makes the darkening tens of
    /// mmag, not the <1 mmag a linear phase law would give: measured over
    /// 500 BB21-posterior orbits × the 40-night DES grid (DE421 Earth
    /// states), max Δm = 0.0763 mag at max α = 0.329° (the population's
    /// minimum heliocentric distance is 176 AU; at the fiducial 500 AU the
    /// surge is 0.039 mag). The survey-recovery pins absorb this: a
    /// ≤0.08 mag dimming on the ephemeris path leaves the computed DES
    /// recovery and the combined exclusion fractions inside their existing
    /// regression tolerances (verified — no pins re-pinned). Kernel-gated.
    #[test]
    fn phase_induced_dmag_small_over_reference_population() {
        use p9_core::analysis::photometry::{hg_phase_factor, DEFAULT_SLOPE_G};
        use p9_core::coords::observer::EphemerisEarth;
        use p9_core::data::reference_population::generate_reference_population;
        use rand::SeedableRng;

        let mut earth = match EphemerisEarth::try_new() {
            Ok(e) => e,
            Err(_) => {
                eprintln!("skipping: no ephemeris kernel in starfield cache");
                return;
            }
        };
        let nights =
            crate::survey_model::DesSurvey::default().observing_epochs_with_earth(&mut earth);
        let mut rng = rand::rngs::StdRng::seed_from_u64(41);
        let population = generate_reference_population(500, &mut rng);

        let mut max_dm = 0.0_f64;
        let mut max_alpha = 0.0_f64;
        let mut min_r = f64::INFINITY;
        for obj in &population {
            let elem = OrbitalElements {
                a: obj.a,
                e: obj.e,
                i: obj.i,
                omega: obj.omega,
                omega_big: obj.omega_big,
                mean_anomaly: obj.mean_anomaly,
            };
            for (_, t_days, state) in &nights {
                let alpha = phase_angle_with_earth(&elem, state, *t_days);
                let dm = -2.5 * hg_phase_factor(DEFAULT_SLOPE_G, alpha).log10();
                max_alpha = max_alpha.max(alpha);
                max_dm = max_dm.max(dm);
                min_r = min_r.min(heliocentric_position(&elem, *t_days).norm());
            }
        }
        eprintln!(
            "max Δm = {max_dm:.4} mag, max α = {:.4}°, min r = {min_r:.1} AU",
            max_alpha.to_degrees()
        );

        // Geometry: phase angle bounded by the annual maximum
        // asin(1 AU / min r) = asin(1/176) = 0.326° (plus Earth's
        // eccentricity); measured 0.329°.
        assert!(
            max_alpha.to_degrees() < 0.35,
            "max phase angle = {:.4}°",
            max_alpha.to_degrees()
        );
        assert!(min_r > 150.0, "min r = {min_r:.1} AU");
        // Photometry: the opposition surge contributes a few tens of mmag
        // (measured 0.0763 mag); it must never grow to a level that would
        // move survey completeness materially (>0.1 mag).
        assert!(
            (0.03..0.10).contains(&max_dm),
            "max phase-induced Δm = {max_dm:.4} mag"
        );
    }

    #[test]
    fn test_parallax_amplitude_scales_inversely_with_distance() {
        // Over one year the apparent RA of a distant object oscillates with
        // amplitude ~ (1 AU / r) radians.
        let elem = OrbitalElements {
            a: 400.0,
            e: 0.0,
            i: 0.0,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 90.0 * DEG2RAD,
        };
        let mut ra_min = f64::INFINITY;
        let mut ra_max = f64::NEG_INFINITY;
        for k in 0..=24 {
            let (ra, _, _) = apparent_position_deg(&elem, k as f64 * YEAR_DAYS / 24.0);
            ra_min = ra_min.min(ra);
            ra_max = ra_max.max(ra);
        }
        let half_amplitude = (ra_max - ra_min) / 2.0;
        let expected = (1.0_f64 / 400.0).asin() / DEG2RAD; // ~0.143 deg
        assert!(
            (half_amplitude - expected).abs() < 0.05,
            "parallax half-amplitude = {half_amplitude:.3} deg, expected ~{expected:.3}"
        );
    }
}
