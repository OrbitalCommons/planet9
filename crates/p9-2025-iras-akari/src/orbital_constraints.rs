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

use p9_core::constants::{GM_SUN, RAD2DEG, YEAR_DAYS};

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
    /// Implied heliocentric distance (AU) from circular-orbit proper motion inversion
    pub r_circular_au: f64,
    /// Implied eccentricity if distance = r_circular and semi-major axis = a
    pub implied_eccentricity: f64,
    /// Whether this (a, e) combination is physically plausible
    pub plausible: bool,
}

/// Compute two-epoch constraints from the candidate pair.
pub fn derive_constraints(pair: &CandidatePair) -> TwoEpochConstraints {
    let baseline = pair.akari_epoch - pair.iras_epoch;

    let ra1 = pair.iras_ra_deg.to_radians();
    let dec1 = pair.iras_dec_deg.to_radians();
    let ra2 = pair.akari_ra_deg.to_radians();
    let dec2 = pair.akari_dec_deg.to_radians();

    // Angular separation (Vincenty)
    let d_ra = ra2 - ra1;
    let num = ((dec2.cos() * d_ra.sin()).powi(2)
        + (dec1.cos() * dec2.sin() - dec1.sin() * dec2.cos() * d_ra.cos()).powi(2))
    .sqrt();
    let den = dec1.sin() * dec2.sin() + dec1.cos() * dec2.cos() * d_ra.cos();
    let sep_rad = num.atan2(den);
    let sep_arcmin = sep_rad.to_degrees() * 60.0;

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

/// Invert the circular-orbit proper motion to get implied distance.
///
/// µ = sqrt(GM) / r^(3/2) * conversion
/// r = (sqrt(GM) * conversion / µ_per_day)^(2/3)
pub fn implied_distance_circular(proper_motion_arcmin_yr: f64) -> f64 {
    let mu_per_day = proper_motion_arcmin_yr / (RAD2DEG * 60.0 * YEAR_DAYS);
    // mu_per_day = sqrt(GM) / r^(3/2)
    // r^(3/2) = sqrt(GM) / mu_per_day
    (GM_SUN.sqrt() / mu_per_day).powf(2.0 / 3.0)
}

/// Explore the (a, e) parameter space consistent with the observed proper motion.
///
/// For each assumed semi-major axis `a`, compute the heliocentric distance `r`
/// that would produce the observed proper motion for a circular orbit, then
/// derive the eccentricity needed: e = |1 - r/a| (approximately).
pub fn explore_orbit_space(
    constraints: &TwoEpochConstraints,
    a_values: &[f64],
) -> Vec<DistanceOrbitConstraint> {
    a_values
        .iter()
        .map(|&a| {
            // For the observed µ, the vis-viva equation gives:
            // µ = v/r = sqrt(GM*(2/r - 1/a)) / r
            // We need to solve for r given µ and a.
            // Use numerical inversion.
            let r = solve_distance_for_motion(constraints.proper_motion_arcmin_yr, a);
            let e = if r <= a {
                // r = a(1-e*cos(E)), but simplest constraint: perihelion <= r <= aphelion
                // e_min such that a(1-e) <= r: e >= 1 - r/a
                // e_max such that r <= a(1+e): e >= r/a - 1
                (1.0 - r / a).abs()
            } else {
                // r > a: only possible if r is near aphelion
                r / a - 1.0
            };
            let plausible = e < 1.0 && r > 0.0 && a > 0.0;
            DistanceOrbitConstraint {
                a_au: a,
                r_circular_au: r,
                implied_eccentricity: e,
                plausible,
            }
        })
        .collect()
}

/// Solve for heliocentric distance r given observed proper motion and semi-major axis a.
///
/// µ = sqrt(GM * (2/r - 1/a)) / r  →  µ²r² = GM(2/r - 1/a)
/// µ²r³ + GM/a = 2*GM*... This is cubic in r; we solve numerically.
fn solve_distance_for_motion(mu_arcmin_yr: f64, a_au: f64) -> f64 {
    let mu_rad_per_day = mu_arcmin_yr / (RAD2DEG * 60.0 * YEAR_DAYS);
    let mu2 = mu_rad_per_day * mu_rad_per_day;

    // f(r) = mu² * r² - GM*(2/r - 1/a)
    // f(r) = mu² * r² - 2*GM/r + GM/a
    // We want f(r) = 0, with r > 0.
    // Multiply by r: mu² * r³ + GM/a * r - 2*GM = 0

    let gm = GM_SUN;
    // Coefficients of cubic: mu²*r³ + 0*r² + (GM/a)*r - 2*GM = 0
    // Use Newton's method starting from the circular estimate.
    let r_circ = (gm.sqrt() / mu_rad_per_day).powf(2.0 / 3.0);
    let mut r = r_circ;

    for _ in 0..100 {
        let f = mu2 * r.powi(3) + (gm / a_au) * r - 2.0 * gm;
        let fp = 3.0 * mu2 * r.powi(2) + gm / a_au;
        if fp.abs() < 1e-30 {
            break;
        }
        let dr = f / fp;
        r -= dr;
        if r < 1.0 {
            r = 1.0;
        }
        if dr.abs() < 1e-10 {
            break;
        }
    }
    r
}

/// Ecliptic coordinates (approximate) from RA/Dec.
///
/// Uses the obliquity of the ecliptic ε = 23.4393° to convert
/// equatorial coordinates to ecliptic coordinates.
pub fn equatorial_to_ecliptic(ra_deg: f64, dec_deg: f64) -> (f64, f64) {
    let eps = 23.4393_f64.to_radians();
    let ra = ra_deg.to_radians();
    let dec = dec_deg.to_radians();

    let sin_lambda = ra.sin() * eps.cos() + dec.tan() * eps.sin();
    let cos_lambda = ra.cos();
    let lambda = sin_lambda.atan2(cos_lambda);

    let sin_beta = dec.sin() * eps.cos() - dec.cos() * eps.sin() * ra.sin();
    let beta = sin_beta.asin();

    (lambda.to_degrees().rem_euclid(360.0), beta.to_degrees())
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn implied_circular_distance() {
        let pair = CandidatePair::paper_candidate();
        let c = derive_constraints(&pair);
        let d = implied_distance_circular(c.proper_motion_arcmin_yr);
        // For µ ≈ 2.06'/yr circular orbit, d ≈ 420 AU
        assert!(
            d > 350.0 && d < 500.0,
            "implied circular distance = {d:.0} AU"
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

        // For a = 700 AU (nominal P9), check eccentricity is moderate
        let a700 = orbits
            .iter()
            .find(|o| (o.a_au - 700.0).abs() < 1.0)
            .unwrap();
        assert!(a700.plausible, "a=700 AU should be plausible");
        assert!(
            a700.implied_eccentricity < 0.8,
            "e = {:.2} at a=700, expected moderate",
            a700.implied_eccentricity
        );
    }

    #[test]
    fn ecliptic_coordinates_of_candidate() {
        let pair = CandidatePair::paper_candidate();
        // Midpoint position
        let mid_ra = (pair.iras_ra_deg + pair.akari_ra_deg) / 2.0;
        let mid_dec = (pair.iras_dec_deg + pair.akari_dec_deg) / 2.0;
        let (ecl_lon, ecl_lat) = equatorial_to_ecliptic(mid_ra, mid_dec);
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

    #[test]
    fn solve_distance_recovers_circular_orbit() {
        // Compute µ for a circular orbit at d=500 AU, then solve with a=500.
        // Should recover r=500.
        let d = 500.0;
        let mu = crate::proper_motion::annual_proper_motion_circular(d);
        let r_solved = solve_distance_for_motion(mu, d);
        assert!(
            (r_solved - d).abs() < 1.0,
            "r_solved={r_solved:.1}, expected {d}"
        );
    }
}
