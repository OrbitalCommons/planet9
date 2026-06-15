//! DES observing-night grid anchor and per-night Earth states.
//!
//! The generic orbit -> apparent sky-position helpers
//! (`apparent_position_deg`, `apparent_position_with_earth_deg`,
//! `phase_angle_with_earth`) live in `p9_core::coords::sky`. This module keeps
//! only the DES-specific night-grid pieces: the survey-start anchor and the
//! per-night Earth-state precompute.
//!
//! The footprint of `survey_model` is defined in *equatorial* coordinates
//! (RA/dec), so orbital sky positions are rotated from the heliocentric
//! ecliptic frame to the equatorial frame (in core's `coords::sky`) before any
//! footprint test.
//!
//! Two Earth models are available in core:
//! - `apparent_position_with_earth_deg` (primary): a real Earth state from
//!   `p9_core::coords::observer` (DE421 via `EphemerisEarth`), with
//!   light-time retardation and stellar aberration. Per-night Earth states are
//!   precomputed once per survey grid (`earth_states_at`) so the kernel is not
//!   hit in the per-orbit loops.
//! - `apparent_position_deg` (analytic fallback): a circular, coplanar 1 AU
//!   Earth with phase zero at the survey start, accurate to ~0.03 deg in
//!   parallax — kept for the speed-critical inner MC loops.

use p9_core::coords::observer::{EarthProvider, EarthState, Time, Timescale};

/// Start of the DES wide survey (Y1 began 2013-08-31; the six observing
/// seasons run through early 2019). Anchor for the night grid's day
/// offsets.
pub fn survey_start(ts: &Timescale) -> Time {
    ts.utc((2013, 8, 31))
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

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::coords::sky::{
        angular_distance, apparent_position_deg, apparent_position_with_earth_deg,
        phase_angle_with_earth,
    };
    use p9_core::types::OrbitalElements;

    #[test]
    fn ephemeris_apparent_position_close_to_analytic() {
        use nalgebra::Vector3;
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
            let sep = angular_distance(
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
                let (_, _, r_helio) = apparent_position_with_earth_deg(&elem, state, *t_days);
                min_r = min_r.min(r_helio);
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
}
