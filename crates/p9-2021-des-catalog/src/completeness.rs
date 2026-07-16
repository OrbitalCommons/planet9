//! Distant-object detection completeness over the DES footprint, and a
//! Planet-Nine footprint-coverage fraction.
//!
//! The six-year DES catalog is complete to r ≈ 23.8 (the shared DES depth in
//! `p9_core::analysis::surveys`). For a distant object the apparent magnitude
//! is dominated by the inverse-fourth-power dimming with heliocentric distance,
//! so the completeness curve translates directly into a *maximum detectable
//! distance* at a given absolute magnitude H. Combined with the fraction of
//! the sky the ~5000 deg² footprint covers (modulo a galactic-plane mask), this
//! gives the catalog's distant-object completeness — the quantity the paper's
//! detection-efficiency characterization feeds into the clustering / exclusion
//! analyses.
//!
//! We compute:
//!   * `apparent_mag` — H, heliocentric and geocentric distance to r;
//!   * `magnitude_completeness` — logistic detection efficiency at r = 23.8
//!     (the Belyakov/Bernardinelli logit form reused from `p9-2022-des`);
//!   * `sky_coverage_fraction` — the 5000 deg² footprint as a fraction of the
//!     full sky after a galactic-plane cut;
//!   * `footprint_completeness` — a Monte-Carlo distant-object completeness
//!     over an isotropic shell: P(in footprint, off the galactic plane) ×
//!     P(bright enough at r = 23.8). A number in (0, 1).

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use p9_core::analysis::surveys::des_footprint_contains;
use p9_core::analysis::surveys::limiting_magnitude;
use p9_core::constants::DEG2RAD;

/// Shared DES r-band 50%-completeness depth (23.8) from the p9-core survey
/// table — the single source of truth.
pub fn des_depth_r() -> f64 {
    limiting_magnitude("DES").expect("DES depth in p9-core survey table")
}

/// Apparent r magnitude of a distant solar-system object:
///   r = H + 5 log10(d_helio · d_geo)        [d in AU]
/// using the standard outer-system approximation d_geo ≈ d_helio − 1.
pub fn apparent_mag(h_mag: f64, d_helio_au: f64) -> f64 {
    let d_geo = (d_helio_au - 1.0).max(0.1);
    h_mag + 5.0 * (d_helio_au * d_geo).log10()
}

/// Logistic single-epoch detection completeness at apparent magnitude `m`,
/// anchored on the DES r depth m50 = 23.8 (Belyakov/Bernardinelli logit form,
/// asymptotic completeness c = 0.95, slope k = 3.0, ~0.7 mag transition).
pub fn magnitude_completeness(apparent_mag: f64) -> f64 {
    let c = 0.95;
    let k = 3.0;
    let m50 = des_depth_r();
    c / (1.0 + (k * (apparent_mag - m50)).exp())
}

/// Maximum heliocentric distance (AU) at which an object of absolute magnitude
/// `h_mag` reaches a target completeness `target` (default-use ~0.5 at m50).
/// Solved by bisection on the monotone apparent-magnitude relation.
pub fn max_detectable_distance(h_mag: f64, target: f64) -> f64 {
    let mut lo = 2.0_f64;
    let mut hi = 5000.0_f64;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if magnitude_completeness(apparent_mag(h_mag, mid)) > target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Equatorial (RA, Dec) radians to galactic latitude b (radians). Uses the
/// J2000 north galactic pole (RA 192.859°, Dec 27.128°).
pub fn galactic_latitude(ra: f64, dec: f64) -> f64 {
    let ra_ngp = 192.859_48 * DEG2RAD;
    let dec_ngp = 27.128_25 * DEG2RAD;
    (dec.sin() * dec_ngp.sin() + dec.cos() * dec_ngp.cos() * (ra - ra_ngp).cos()).asin()
}

/// Fraction of the full sphere covered by the DES footprint after masking the
/// galactic plane (|b| < `b_cut_deg`), by uniform-sphere Monte Carlo.
pub fn sky_coverage_fraction(b_cut_deg: f64, seed: u64, n: usize) -> f64 {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut hits = 0usize;
    for _ in 0..n {
        let ra_deg = rng.gen::<f64>() * 360.0;
        let z: f64 = rng.gen_range(-1.0..1.0);
        let dec = z.asin();
        let dec_deg = dec / DEG2RAD;
        let b = galactic_latitude(ra_deg * DEG2RAD, dec) / DEG2RAD;
        if b.abs() >= b_cut_deg && des_footprint_contains(ra_deg, dec_deg) {
            hits += 1;
        }
    }
    hits as f64 / n as f64
}

/// Distant-object completeness over the footprint for a population at fixed
/// absolute magnitude `h_mag` distributed isotropically on a shell at
/// heliocentric distance `d_helio_au`: the product of footprint coverage
/// (off the galactic plane) and the per-object magnitude completeness at that
/// distance. Returns a fraction in (0, 1).
pub fn footprint_completeness(
    h_mag: f64,
    d_helio_au: f64,
    b_cut_deg: f64,
    seed: u64,
    n: usize,
) -> f64 {
    let mag = apparent_mag(h_mag, d_helio_au);
    let bright = magnitude_completeness(mag);
    let coverage = sky_coverage_fraction(b_cut_deg, seed, n);
    (coverage * bright).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn des_depth_is_23_8() {
        assert_relative_eq!(des_depth_r(), 23.8, epsilon = 1e-10);
    }

    #[test]
    fn completeness_is_half_at_m50() {
        // At apparent r = 23.8 the logit gives c/2 = 0.475.
        assert_relative_eq!(magnitude_completeness(23.8), 0.475, epsilon = 1e-9);
    }

    #[test]
    fn bright_objects_fully_complete_faint_objects_lost() {
        assert!(magnitude_completeness(20.0) > 0.94, "bright should be ~c");
        assert!(magnitude_completeness(26.0) < 0.02, "faint should be lost");
    }

    #[test]
    fn apparent_mag_dims_with_distance() {
        // A Sedna-like H = 1.6 at 80 AU vs 500 AU.
        let near = apparent_mag(1.6, 80.0);
        let far = apparent_mag(1.6, 500.0);
        assert!(far > near);
        // ~inverse-fourth-power: ~10 mag fainter per decade in distance.
        assert!(
            (far - near) > 7.0 && (far - near) < 9.0,
            "delta = {}",
            far - near
        );
    }

    #[test]
    fn p9_reach_falls_short_of_hypothesized_distance() {
        // A Planet Nine at H ~ 4 is detectable at the 50% level only out to
        // ~100 AU in DES — well inside its hypothesized ~400-800 AU semi-major
        // axis. This is the honest physics behind "no Planet Nine in DES":
        // the r = 23.8 depth simply does not reach a P9 at its predicted
        // distance unless it is unusually bright/large.
        let d = max_detectable_distance(4.0, 0.475);
        assert!(d > 50.0 && d < 200.0, "max detectable distance = {d:.0} AU");
        assert!(
            d < 400.0,
            "DES reach {d:.0} AU should fall short of P9's ~400+ AU"
        );
    }

    #[test]
    fn footprint_coverage_fraction_is_a_probability() {
        // ~5000 deg^2 is ~12% of the sky; a galactic cut trims it modestly.
        let f = sky_coverage_fraction(10.0, 2021, 400_000);
        assert!(f > 0.0 && f < 1.0, "coverage = {f}");
        // 5000/41253 = 0.121 full-sky; with the galactic mask it lands below.
        assert!(
            (f - 0.121).abs() < 0.04,
            "coverage = {f:.3} vs 0.121 full-sky"
        );
    }

    #[test]
    fn distant_object_completeness_in_open_unit_interval() {
        // A distant, faint extreme TNO (H = 7) at 300 AU within DES.
        let c = footprint_completeness(7.0, 300.0, 10.0, 99, 200_000);
        assert!(c > 0.0 && c < 1.0, "completeness = {c}");
    }

    #[test]
    fn footprint_completeness_matches_factorization() {
        // Cross-check: footprint completeness = coverage * magnitude completeness.
        let h = 6.0;
        let d = 250.0;
        let cov = sky_coverage_fraction(10.0, 55, 300_000);
        let comp = footprint_completeness(h, d, 10.0, 55, 300_000);
        let expected = cov * magnitude_completeness(apparent_mag(h, d));
        assert_relative_eq!(comp, expected, epsilon = 1e-12);
    }
}
