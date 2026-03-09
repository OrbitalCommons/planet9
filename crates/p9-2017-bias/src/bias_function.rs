//! Observational bias function computation.
//!
//! The bias function B(ϖ, Ω, i) gives the relative probability that an object
//! with orbital elements (a, e, H) would be discovered at a given (ϖ, Ω, i).
//!
//! Key insight from the paper: distant eccentric objects are only observable
//! near perihelion, so the bias is dominated by sky coverage near the
//! perihelion direction.
//!
//! TODO: Full implementation requires MPC discovery catalog with survey
//! pointings and depths. Currently uses a simplified perihelion-distance
//! bias model.

use std::f64::consts::PI;

use p9_core::constants::DEG2RAD;

/// Simplified bias function based on perihelion visibility.
///
/// For highly eccentric orbits, objects spend most of their time near aphelion
/// and are only bright enough to detect near perihelion. The bias is primarily
/// a function of where perihelion falls on the sky.
///
/// This simplified model assumes:
/// 1. Object is only detectable within ±30° of perihelion along the orbit
/// 2. Detection probability scales as r^{-4} (apparent brightness ∝ r^{-4})
/// 3. Sky coverage is uniform (simplified — real surveys are patchy)
///
/// Returns a bias weight in [0, 1].
pub fn perihelion_bias(
    a: f64,
    e: f64,
    _varpi: f64,
    omega: f64,
    i: f64,
    galactic_latitude_cutoff: f64,
) -> f64 {
    let q = a * (1.0 - e);

    // Maximum heliocentric distance for detection
    // Paper uses ~90 AU as practical limit
    let r_max_detect = 90.0;

    if q > r_max_detect {
        return 0.0;
    }

    // Brightness factor: brighter at perihelion
    let brightness_factor = (r_max_detect / q).powi(4).min(100.0) / 100.0;

    // Galactic plane avoidance: surveys avoid |b| < cutoff
    // Perihelion ecliptic latitude depends on i and ω
    let beta_peri = omega.sin() * i.sin();
    let galactic_penalty = if beta_peri.abs() < galactic_latitude_cutoff.sin() {
        0.3 // Reduced but not zero (some surveys cover the plane)
    } else {
        1.0
    };

    brightness_factor * galactic_penalty
}

/// Compute bias weights for a grid of (ϖ, Ω) values.
///
/// Returns a 2D array of bias values for each (ϖ, Ω) combination,
/// with inclination marginalized over a sin(i) distribution.
pub fn bias_grid(a: f64, e: f64, n_varpi: usize, n_omega_big: usize) -> Vec<Vec<f64>> {
    let d_varpi = 2.0 * PI / n_varpi as f64;
    let d_omega = 2.0 * PI / n_omega_big as f64;
    let galactic_cutoff = 10.0 * DEG2RAD;

    let n_i_samples = 10;

    let mut grid = vec![vec![0.0; n_omega_big]; n_varpi];

    for j in 0..n_varpi {
        let varpi = (j as f64 + 0.5) * d_varpi;
        for k in 0..n_omega_big {
            let omega_big = (k as f64 + 0.5) * d_omega;
            let omega = varpi - omega_big;

            // Marginalize over inclination with sin(i) prior
            let mut total_weight = 0.0;
            let mut total_sin = 0.0;
            for l in 0..n_i_samples {
                let i = (l as f64 + 0.5) / n_i_samples as f64 * PI;
                let sin_i = i.sin();
                let w = perihelion_bias(a, e, varpi, omega, i, galactic_cutoff);
                total_weight += w * sin_i;
                total_sin += sin_i;
            }

            grid[j][k] = total_weight / total_sin;
        }
    }

    grid
}

/// Compute the effective bias weight for a specific KBO.
///
/// This is the bias at the object's actual position relative to the
/// average bias across all positions.
pub fn bias_weight(a: f64, e: f64, varpi: f64, omega: f64, i: f64) -> f64 {
    let galactic_cutoff = 10.0 * DEG2RAD;
    let actual = perihelion_bias(a, e, varpi, omega, i, galactic_cutoff);

    // Average bias over all ϖ values
    let n_samples = 36;
    let avg: f64 = (0..n_samples)
        .map(|j| {
            let v = 2.0 * PI * j as f64 / n_samples as f64;
            let w = v - (varpi - omega); // Keep ω fixed, vary ϖ
            perihelion_bias(a, e, v, w, i, galactic_cutoff)
        })
        .sum::<f64>()
        / n_samples as f64;

    if avg > 0.0 {
        actual / avg
    } else {
        1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perihelion_bias_detectable() {
        // Typical Sedna-like object: a=500, e=0.85, q=75 AU
        let b = perihelion_bias(500.0, 0.85, 0.0, 0.0, 0.1, 10.0 * DEG2RAD);
        assert!(b > 0.0, "Sedna-like object should be detectable");
        assert!(b <= 1.0);
    }

    #[test]
    fn test_perihelion_bias_too_distant() {
        // Object with q > 90 AU is undetectable
        let b = perihelion_bias(1000.0, 0.85, 0.0, 0.0, 0.1, 10.0 * DEG2RAD);
        // q = 1000 * 0.15 = 150 AU > 90 AU
        assert!(
            b == 0.0,
            "Object with q=150 AU should be undetectable, got b={:.4}",
            b
        );
    }

    #[test]
    fn test_bias_grid_shape() {
        let grid = bias_grid(300.0, 0.85, 36, 36);
        assert_eq!(grid.len(), 36);
        assert_eq!(grid[0].len(), 36);

        // All values should be non-negative
        for row in &grid {
            for &val in row {
                assert!(val >= 0.0);
            }
        }
    }

    #[test]
    fn test_bias_weight_around_unity() {
        // For a typical object, bias weight should be around 1.0
        let w = bias_weight(300.0, 0.85, 1.0, 5.5, 0.3);
        assert!(w > 0.0);
        assert!(w < 10.0, "Weight {:.2} seems too high", w);
    }
}
