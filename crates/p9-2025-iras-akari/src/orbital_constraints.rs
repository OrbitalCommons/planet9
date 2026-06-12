//! Orbital constraints derived from the candidate pair.
//!
//! With only two sky positions separated by 23 years, a full Keplerian
//! orbit (6 elements) cannot be determined — the paper notes that 9
//! detections are required. However, the angular separation and direction
//! of motion constrain a combination of distance, velocity, and sky
//! trajectory.
//!
//! From the candidate pair we can derive:
//! - Direction of motion on the sky (position angle)
//! - Annual proper motion (angular rate)
//! - Implied heliocentric distance (given assumed orbit shape)
//! - Constraints on (a, e) pairs consistent with the observed motion

use crate::proper_motion;
use crate::survey_model::angular_separation_arcmin;

/// The candidate pair from Phan et al. (2025) Table 2.
pub struct CandidatePair {
    /// IRAS source RA (degrees)
    pub iras_ra_deg: f64,
    /// IRAS source Dec (degrees)
    pub iras_dec_deg: f64,
    /// IRAS epoch (fractional year)
    pub iras_epoch: f64,
    /// AKARI source RA (degrees)
    pub akari_ra_deg: f64,
    /// AKARI source Dec (degrees)
    pub akari_dec_deg: f64,
    /// AKARI epoch (fractional year)
    pub akari_epoch: f64,
    /// IRAS 60 µm flux (Jy)
    pub iras_flux_60: f64,
    /// IRAS 100 µm flux (Jy)
    pub iras_flux_100: f64,
    /// AKARI 65 µm flux (Jy)
    pub akari_flux_65: f64,
    /// AKARI 90 µm flux (Jy)
    pub akari_flux_90: f64,
}

impl CandidatePair {
    /// The single good candidate from Table 2.
    pub fn paper_candidate() -> Self {
        Self {
            iras_ra_deg: 35.74075,
            iras_dec_deg: -48.5125,
            iras_epoch: 1983.5,
            akari_ra_deg: 35.18379,
            akari_dec_deg: -49.2135,
            akari_epoch: 2006.5,
            iras_flux_60: 0.24,
            iras_flux_100: 0.52,
            akari_flux_65: 0.09,
            akari_flux_90: 0.27,
        }
    }
}

/// Constraints derivable from a two-epoch observation.
#[derive(Debug, Clone)]
pub struct TwoEpochConstraints {
    /// Angular separation in arcminutes
    pub separation_arcmin: f64,
    /// Annual proper motion in arcmin/yr
    pub proper_motion_arcmin_yr: f64,
    /// Position angle of motion (degrees, N through E)
    pub position_angle_deg: f64,
    /// Epoch baseline in years
    pub baseline_years: f64,
    /// Direction of RA change (negative = westward)
    pub delta_ra_arcmin: f64,
    /// Direction of Dec change (negative = southward)
    pub delta_dec_arcmin: f64,
}

/// Distance-orbit constraints for a grid of assumed semi-major axes.
#[derive(Debug, Clone)]
pub struct DistanceOrbitConstraint {
    /// Assumed semi-major axis (AU)
    pub a_au: f64,
    /// Heliocentric distance at observation (AU), perihelion branch
    pub r_au: f64,
    /// Eccentricity required to reproduce the observed sky rate at
    /// perihelion of an orbit with this semi-major axis
    pub implied_eccentricity: f64,
    /// Whether this (a, e) combination is physically plausible
    pub plausible: bool,
}

/// Compute two-epoch constraints from the candidate pair.
pub fn derive_constraints(pair: &CandidatePair) -> TwoEpochConstraints {
    let baseline = pair.akari_epoch - pair.iras_epoch;

    // Angular separation: single shared Vincenty implementation.
    let sep_arcmin = angular_separation_arcmin(
        pair.iras_ra_deg,
        pair.iras_dec_deg,
        pair.akari_ra_deg,
        pair.akari_dec_deg,
    );

    let ra1 = pair.iras_ra_deg.to_radians();
    let dec1 = pair.iras_dec_deg.to_radians();
    let ra2 = pair.akari_ra_deg.to_radians();
    let dec2 = pair.akari_dec_deg.to_radians();
    let d_ra = ra2 - ra1;

    // Position angle (bearing from IRAS position to AKARI position)
    let pa_y = dec2.cos() * d_ra.sin();
    let pa_x = dec1.cos() * dec2.sin() - dec1.sin() * dec2.cos() * d_ra.cos();
    let pa_rad = pa_y.atan2(pa_x);
    let pa_deg = pa_rad.to_degrees().rem_euclid(360.0);

    // RA and Dec components of motion
    let delta_ra_deg =
        (pair.akari_ra_deg - pair.iras_ra_deg) * pair.iras_dec_deg.to_radians().cos();
    let delta_dec_deg = pair.akari_dec_deg - pair.iras_dec_deg;

    TwoEpochConstraints {
        separation_arcmin: sep_arcmin,
        proper_motion_arcmin_yr: sep_arcmin / baseline,
        position_angle_deg: pa_deg,
        baseline_years: baseline,
        delta_ra_arcmin: delta_ra_deg * 60.0,
        delta_dec_arcmin: delta_dec_deg * 60.0,
    }
}

/// Invert the circular-orbit, parallax-free proper motion to get the
/// implied (lower-bound) distance. Delegates to the single shared
/// implementation in `proper_motion`.
pub fn implied_distance_circular(proper_motion_arcmin_yr: f64) -> f64 {
    proper_motion::implied_distance_circular(proper_motion_arcmin_yr)
}

/// Explore the (a, e) parameter space consistent with the observed proper
/// motion, assuming the object is observed near perihelion (the favorable
/// configuration: at perihelion the velocity is purely transverse and the
/// sky rate is maximal for a given orbit).
///
/// For each assumed semi-major axis `a`, solve for the eccentricity at
/// which the perihelion transverse rate
///
///   µ(e) = √(GM·a(1−e²)) / (a(1−e))²
///
/// (monotonically increasing in e) matches the observed rate, by
/// bisection. No solution with e ∈ [0, 0.95] → implausible.
pub fn explore_orbit_space(
    constraints: &TwoEpochConstraints,
    a_values: &[f64],
) -> Vec<DistanceOrbitConstraint> {
    let mu_obs = constraints.proper_motion_arcmin_yr;
    a_values
        .iter()
        .map(|&a| {
            let rate_at = |e: f64| {
                let r = a * (1.0 - e);
                proper_motion::transverse_rate_arcmin_yr(r, a, e)
            };
            const E_LIM: f64 = 0.95;
            if rate_at(0.0) > mu_obs || rate_at(E_LIM) < mu_obs {
                // Even circular at r = a is too fast, or e = 0.95 too slow.
                return DistanceOrbitConstraint {
                    a_au: a,
                    r_au: f64::NAN,
                    implied_eccentricity: f64::NAN,
                    plausible: false,
                };
            }
            let (mut lo, mut hi) = (0.0_f64, E_LIM);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if rate_at(mid) < mu_obs {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let e = 0.5 * (lo + hi);
            DistanceOrbitConstraint {
                a_au: a,
                r_au: a * (1.0 - e),
                implied_eccentricity: e,
                plausible: true,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::coords::sky::equatorial_to_ecliptic_deg;

    #[test]
    fn candidate_separation_matches_paper() {
        let pair = CandidatePair::paper_candidate();
        let c = derive_constraints(&pair);
        assert!(
            (c.separation_arcmin - 47.46).abs() < 0.5,
            "separation = {:.2}', expected ~47.46'",
            c.separation_arcmin
        );
    }

    #[test]
    fn candidate_proper_motion() {
        let pair = CandidatePair::paper_candidate();
        let c = derive_constraints(&pair);
        let expected_mu = 47.46 / 23.0; // ~2.06'/yr
        assert!(
            (c.proper_motion_arcmin_yr - expected_mu).abs() < 0.1,
            "µ = {:.2}'/yr, expected ~{:.2}",
            c.proper_motion_arcmin_yr,
            expected_mu
        );
    }

    #[test]
    fn candidate_moves_south_and_west() {
        let pair = CandidatePair::paper_candidate();
        let c = derive_constraints(&pair);
        // AKARI RA < IRAS RA → moved west (negative RA direction)
        assert!(c.delta_ra_arcmin < 0.0, "should move westward");
        // AKARI Dec < IRAS Dec → moved south
        assert!(c.delta_dec_arcmin < 0.0, "should move southward");
    }

    #[test]
    fn candidate_position_angle() {
        let pair = CandidatePair::paper_candidate();
        let c = derive_constraints(&pair);
        // Moving south-west → PA should be roughly 190°–220°
        assert!(
            c.position_angle_deg > 180.0 && c.position_angle_deg < 250.0,
            "PA = {:.1}°, expected ~190-220° (south-west)",
            c.position_angle_deg
        );
    }

    #[test]
    fn implied_circular_distance_is_lower_bound() {
        let pair = CandidatePair::paper_candidate();
        let c = derive_constraints(&pair);
        let d_circ = implied_distance_circular(c.proper_motion_arcmin_yr);
        // Circular/parallax-free inversion: ~420 AU (the lower bound).
        assert!(
            d_circ > 350.0 && d_circ < 500.0,
            "implied circular distance = {d_circ:.0} AU"
        );
        // The geocentric maximum-distance inversion lands in the paper's
        // 500–700 AU range.
        let d_geo = proper_motion::implied_distance(c.separation_arcmin, c.baseline_years);
        assert!(
            d_geo > 500.0 && d_geo < 720.0,
            "implied geocentric distance = {d_geo:.0} AU"
        );
        assert!(d_geo > d_circ);
    }

    #[test]
    fn implied_distance_ephemeris_in_paper_range() {
        use p9_core::coords::observer::EphemerisEarth;
        let mut eph = match EphemerisEarth::try_new() {
            Ok(e) => e,
            Err(_) => {
                eprintln!("skipping: no ephemeris kernel in starfield cache");
                return;
            }
        };
        // Real-Earth inversion of the candidate's 47.46' separation:
        // ~694 AU, +3.4% over the circular-Earth fallback's ~675 AU
        // (real 23.5 yr mid-to-mid baseline + the wider AKARI window),
        // still inside the paper's 500-700(+) AU search interpretation.
        let pair = CandidatePair::paper_candidate();
        let c = derive_constraints(&pair);
        let d_eph = proper_motion::implied_distance_ephemeris(&mut eph, c.separation_arcmin);
        let d_fallback = proper_motion::implied_distance(c.separation_arcmin, c.baseline_years);
        assert!(
            d_eph > 500.0 && d_eph < 720.0,
            "ephemeris implied distance = {d_eph:.0} AU"
        );
        assert!(
            (d_eph - d_fallback).abs() / d_fallback < 0.05,
            "ephemeris {d_eph:.0} vs fallback {d_fallback:.0} AU"
        );
    }

    #[test]
    fn orbit_space_for_nominal_p9() {
        let pair = CandidatePair::paper_candidate();
        let c = derive_constraints(&pair);
        let a_values: Vec<f64> = (300..=1200).step_by(100).map(|a| a as f64).collect();
        let orbits = explore_orbit_space(&c, &a_values);

        // At least some should be plausible
        let n_plausible = orbits.iter().filter(|o| o.plausible).count();
        assert!(
            n_plausible > 0,
            "no plausible orbits found in 300-1200 AU range"
        );

        // For a = 700 AU (nominal P9): perihelion-branch solution at
        // e ≈ 0.26, r ≈ 520 AU.
        let a700 = orbits
            .iter()
            .find(|o| (o.a_au - 700.0).abs() < 1.0)
            .unwrap();
        assert!(a700.plausible, "a=700 AU should be plausible");
        assert!(
            a700.implied_eccentricity > 0.1 && a700.implied_eccentricity < 0.4,
            "e = {:.2} at a=700, expected ~0.26",
            a700.implied_eccentricity
        );
        assert!(
            (a700.r_au - 520.0).abs() < 40.0,
            "r = {:.0} at a=700, expected ~520",
            a700.r_au
        );
        // Round trip: the solved (r, a, e) reproduces the observed rate.
        let mu = proper_motion::transverse_rate_arcmin_yr(
            a700.r_au,
            a700.a_au,
            a700.implied_eccentricity,
        );
        assert!((mu - c.proper_motion_arcmin_yr).abs() < 1e-6);
    }

    #[test]
    fn ecliptic_coordinates_of_candidate() {
        let pair = CandidatePair::paper_candidate();
        // Midpoint position
        let mid_ra = (pair.iras_ra_deg + pair.akari_ra_deg) / 2.0;
        let mid_dec = (pair.iras_dec_deg + pair.akari_dec_deg) / 2.0;
        let (ecl_lon, ecl_lat) = equatorial_to_ecliptic_deg(mid_ra, mid_dec);
        // At RA~35.5°, Dec~-48.9°, ecliptic latitude should be far from
        // the ecliptic plane (|β| >> 0), consistent with P9's predicted
        // high-inclination orbit
        assert!(
            ecl_lat.abs() > 20.0,
            "ecliptic latitude = {ecl_lat:.1}°, expected well off the ecliptic"
        );
        // Longitude should be in the ~0-60° range
        assert!(
            ecl_lon > 0.0 && ecl_lon < 90.0,
            "ecliptic longitude = {ecl_lon:.1}°"
        );
    }

    #[test]
    fn flux_ratio_is_consistent_with_cold_source() {
        let pair = CandidatePair::paper_candidate();
        // A cold blackbody (30-50 K) should have S_100/S_60 > 1
        // (more flux at longer wavelength, near or below the Rayleigh-Jeans regime)
        let iras_ratio = pair.iras_flux_100 / pair.iras_flux_60;
        assert!(
            iras_ratio > 1.0,
            "IRAS S_100/S_60 = {iras_ratio:.2}, expected > 1 for cold source"
        );

        // AKARI S_90/S_65 should similarly be > 1
        let akari_ratio = pair.akari_flux_90 / pair.akari_flux_65;
        assert!(
            akari_ratio > 1.0,
            "AKARI S_90/S_65 = {akari_ratio:.2}, expected > 1 for cold source"
        );
    }
}
