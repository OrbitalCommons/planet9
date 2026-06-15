//! Ephemeris constraint on Planet Nine's *current orbital phase*.
//!
//! Cassini Earth–Saturn range residuals (Fienga et al. 2016) and Saturn-
//! precession bounds (Iorio 2026) favor P9 being near true anomaly
//! ν ≈ 108–129° on the Brown & Batygin orbit right now (and disfavor large
//! complementary arcs). That is a statement about *where on its orbit* P9 sits
//! today — i.e. a prior on the current true anomaly — which directly cuts the
//! sky RA range, not just the declination.
//!
//! The favored interval and tolerance are sourced from
//! `p9_core::data::ephemeris_constraint`; here we turn them into (1) a smooth
//! ν-weight to reweight the otherwise-uniform orbital phase, and (2) the
//! favored-ν sky arc for plotting.

use nalgebra::Vector3;
use p9_core::constants::DEG2RAD;
use p9_core::coords::sky::ecliptic_vec_to_equatorial_deg;
use p9_core::data::ephemeris_constraint::{FAVORED_INTERVAL_DEG, TRUE_ANOMALY_TOLERANCE_DEG};
use p9_core::types::{position_at_true_anomaly, P9Params};

/// Favored true-anomaly interval (deg), from Fienga et al. 2016 via
/// `p9_core::data::ephemeris_constraint`.
pub fn favored_interval_deg() -> (f64, f64) {
    FAVORED_INTERVAL_DEG
}

/// 1σ tolerance (deg) on the favored interval edges.
pub fn tolerance_deg() -> f64 {
    TRUE_ANOMALY_TOLERANCE_DEG / 2.0
}

/// Smallest absolute angular separation (deg) of `nu_deg` from the interval
/// `[lo, hi]` on the circle.
fn distance_to_interval(nu_deg: f64, lo: f64, hi: f64) -> f64 {
    let n = nu_deg.rem_euclid(360.0);
    if n >= lo && n <= hi {
        return 0.0;
    }
    let d1 = ((n - lo + 540.0).rem_euclid(360.0) - 180.0).abs();
    let d2 = ((n - hi + 540.0).rem_euclid(360.0) - 180.0).abs();
    d1.min(d2)
}

/// Ephemeris weight in [0, 1] for a sample at true anomaly `nu_deg`: flat at 1
/// inside the favored interval, Gaussian wings outside with σ = [`tolerance_deg`].
pub fn nu_weight(nu_deg: f64) -> f64 {
    let (lo, hi) = favored_interval_deg();
    let d = distance_to_interval(nu_deg, lo, hi);
    if d == 0.0 {
        1.0
    } else {
        let sigma = tolerance_deg().max(1.0);
        (-0.5 * (d / sigma).powi(2)).exp()
    }
}

/// Sky arc (RA, Dec in deg) traced by the favored-ν zone (interval ± tolerance)
/// on a given orbit — for plotting the "ephemeris says here" band.
pub fn favored_arc(orbit: &P9Params, n: usize) -> Vec<[f64; 2]> {
    let (lo, hi) = favored_interval_deg();
    let tol = tolerance_deg();
    let (a, b) = (lo - tol, hi + tol);
    (0..n)
        .map(|k| {
            let nu = (a + (b - a) * k as f64 / (n.max(2) - 1) as f64) * DEG2RAD;
            let pos: Vector3<f64> = position_at_true_anomaly(orbit, nu).position;
            let (ra, dec) = ecliptic_vec_to_equatorial_deg(&pos);
            [ra, dec]
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weight_peaks_in_favored_interval() {
        let (lo, hi) = favored_interval_deg();
        assert_eq!(nu_weight((lo + hi) / 2.0), 1.0);
        assert_eq!(nu_weight(lo), 1.0);
        // Aphelion (180°) is outside the favored zone → suppressed.
        assert!(nu_weight(180.0) < 0.5);
        // The opposite side of the orbit is strongly suppressed.
        assert!(nu_weight(300.0) < 0.05);
    }

    #[test]
    fn distance_wraps_on_circle() {
        // 350° is 11° from 1°-edge interval... sanity on wrap.
        assert!((distance_to_interval(350.0, 0.0, 10.0) - 10.0).abs() < 1e-9);
        assert_eq!(distance_to_interval(5.0, 0.0, 10.0), 0.0);
    }

    #[test]
    fn favored_arc_is_on_sky() {
        let arc = favored_arc(&P9Params::mcmc_2021(), 16);
        assert_eq!(arc.len(), 16);
        for p in &arc {
            assert!((0.0..=360.0).contains(&p[0]) && (-90.0..=90.0).contains(&p[1]));
        }
    }
}
