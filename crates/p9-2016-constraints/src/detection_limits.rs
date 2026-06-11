//! Detection limit assessment from Brown & Batygin (2016) Section 4.
//!
//! Predicts Planet Nine's apparent brightness and maps which regions of
//! its orbit have been surveyed to sufficient depth.
//!
//! TODO: Full sky coverage analysis requires survey footprint data from
//! starfield-datasources (WISE, CRTS, Pan-STARRS, DES).

use std::f64::consts::PI;

use p9_core::analysis::photometry;
use p9_core::constants::EARTH_RADIUS_KM;
use p9_core::types::P9Params;

/// Predict V-band apparent magnitude of Planet Nine at a given true anomaly,
/// using the shared p9-core photometry (H from diameter and albedo, then the
/// reflected-sunlight m = H + 5 log10(r·Δ) law with Δ at opposition; the
/// phase function is ≈ 1 near opposition at these distances).
pub fn predicted_v_magnitude(p9: &P9Params, true_anomaly: f64, albedo: f64, radius_km: f64) -> f64 {
    let r = p9.a * (1.0 - p9.e * p9.e) / (1.0 + p9.e * true_anomaly.cos());
    let h = photometry::absolute_magnitude(radius_km, albedo);
    photometry::apparent_magnitude(h, r, photometry::opposition_delta(r))
}

/// Planet Nine's physical radius from mass, via the shared Neptune-anchored
/// mass-radius relation (the previous local 3.0 R⊕ · M^0.27 fit made a
/// 10 M⊕ planet larger than Neptune, ~40% too big, skewing detectability
/// optimistic).
pub fn estimate_radius_km(mass_earth: f64) -> f64 {
    photometry::mass_radius_neptunian(mass_earth) * EARTH_RADIUS_KM
}

/// Compute the sky position (ecliptic longitude, latitude) of Planet Nine
/// at a given true anomaly.
///
/// Returns (lambda, beta) in radians.
pub fn sky_position(p9: &P9Params, true_anomaly: f64) -> (f64, f64) {
    let r = p9.a * (1.0 - p9.e * p9.e) / (1.0 + p9.e * true_anomaly.cos());

    // Position in orbital plane
    let x_orb = r * true_anomaly.cos();
    let y_orb = r * true_anomaly.sin();

    // Rotate to ecliptic
    let cos_w = p9.omega.cos();
    let sin_w = p9.omega.sin();
    let cos_i = p9.i.cos();
    let sin_i = p9.i.sin();
    let cos_o = p9.omega_big.cos();
    let sin_o = p9.omega_big.sin();

    let x = (cos_o * cos_w - sin_o * sin_w * cos_i) * x_orb
        + (-cos_o * sin_w - sin_o * cos_w * cos_i) * y_orb;
    let y = (sin_o * cos_w + cos_o * sin_w * cos_i) * x_orb
        + (-sin_o * sin_w + cos_o * cos_w * cos_i) * y_orb;
    let z = sin_w * sin_i * x_orb + cos_w * sin_i * y_orb;

    let dist = (x * x + y * y + z * z).sqrt();
    let lambda = y.atan2(x);
    let beta = (z / dist).asin();

    (lambda, beta)
}

/// Brightness curve: V magnitude vs true anomaly for the full orbit.
pub fn brightness_curve(p9: &P9Params, albedo: f64, n_points: usize) -> Vec<(f64, f64, f64)> {
    let radius_km = estimate_radius_km(p9.mass_earth);

    (0..n_points)
        .map(|i| {
            let nu = (i as f64 / n_points as f64) * 2.0 * PI - PI;
            let v_mag = predicted_v_magnitude(p9, nu, albedo, radius_km);
            let r = p9.a * (1.0 - p9.e * p9.e) / (1.0 + p9.e * nu.cos());
            (nu, r, v_mag)
        })
        .collect()
}

/// Survey depth limits (approximate, from the paper's analysis).
/// Returns the limiting V magnitude for each survey.
pub fn survey_depths() -> Vec<(&'static str, f64)> {
    vec![
        ("WISE W1", 16.5),
        ("WISE W2", 15.5),
        ("CRTS", 19.5),
        ("Pan-STARRS 3π", 21.5),
        ("DES", 23.5),
    ]
}

/// Fraction of the orbit brighter than a given V magnitude limit.
pub fn detectable_fraction(p9: &P9Params, v_limit: f64, albedo: f64) -> f64 {
    let n = 1000;
    let radius_km = estimate_radius_km(p9.mass_earth);
    let n_detectable = (0..n)
        .filter(|&i| {
            let nu = (i as f64 / n as f64) * 2.0 * PI - PI;
            predicted_v_magnitude(p9, nu, albedo, radius_km) < v_limit
        })
        .count();

    n_detectable as f64 / n as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_brightness_at_perihelion_vs_aphelion() {
        let p9 = P9Params::nominal_2016();

        let radius_km = estimate_radius_km(p9.mass_earth);
        let v_peri = predicted_v_magnitude(&p9, 0.0, 0.5, radius_km);
        let v_apo = predicted_v_magnitude(&p9, PI, 0.5, radius_km);

        // Should be much brighter at perihelion
        assert!(
            v_peri < v_apo,
            "Perihelion V={:.1} should be < aphelion V={:.1}",
            v_peri,
            v_apo
        );

        // Paper predicts V ~ 22-25 for nominal parameters
        assert!(
            v_peri > 15.0 && v_peri < 30.0,
            "Perihelion V={:.1} out of expected range",
            v_peri
        );
    }

    #[test]
    fn test_radius_estimate() {
        let r_10 = estimate_radius_km(10.0);
        let r_1 = estimate_radius_km(1.0);

        // 10 Earth mass should have larger radius than 1 Earth mass
        assert!(r_10 > r_1);
        // Should be in Neptune-size range (20,000-30,000 km) for 10 M_Earth
        assert!(
            r_10 > 15_000.0 && r_10 < 40_000.0,
            "R(10 M_Earth) = {:.0} km",
            r_10
        );
    }

    #[test]
    fn test_sky_position_coplanar() {
        let mut p9 = P9Params::nominal_2016();
        p9.i = 0.0;
        p9.omega = 0.0;
        p9.omega_big = 0.0;

        let (_lambda, beta) = sky_position(&p9, 0.0);
        // Coplanar at perihelion with ω=Ω=0 should be on ecliptic
        assert!(
            beta.abs() < 0.01,
            "β should be ~0 for coplanar: {:.4}",
            beta
        );
    }

    #[test]
    fn test_brightness_curve_monotonic_near_perihelion() {
        let p9 = P9Params::nominal_2016();
        let curve = brightness_curve(&p9, 0.5, 100);

        assert_eq!(curve.len(), 100);

        // Find brightest point (lowest V)
        let brightest = curve
            .iter()
            .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap())
            .unwrap();
        // Should be near perihelion (ν ≈ 0)
        assert!(
            brightest.0.abs() < 0.5,
            "Brightest at ν={:.2} rad, expected near 0",
            brightest.0
        );
    }
}
