//! Synthetic reference population generation for Planet Nine.
//!
//! Generates a large population of synthetic Planet Nine realizations
//! drawn from the posterior emulator, computing observable properties
//! (apparent magnitude, heliocentric distance) for each.
//!
//! Photometry comes from `p9_core::analysis::photometry` (Neptune-anchored
//! mass-radius relation, reflected-sunlight distance law, Neptune's geometric
//! albedo 0.41) — the previous inline model invented an albedo ~ U(0.3, 0.7)
//! and a super-Earth mass-radius exponent inconsistent with the rest of the
//! workspace.

use rand::Rng;
use serde::{Deserialize, Serialize};

use p9_core::analysis::photometry::{planet_apparent_magnitude, ALBEDO_NEPTUNE};
use p9_core::types::{solve_kepler, P9Params};

use crate::posterior::{mcmc_2021_posterior, sample_from_posterior, P9Posterior};

/// A single synthetic Planet Nine realization with observable properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceP9 {
    /// Mass in Earth masses
    pub mass: f64,
    /// Semi-major axis in AU
    pub a: f64,
    /// Eccentricity
    pub e: f64,
    /// Inclination in radians
    pub i: f64,
    /// Argument of perihelion in radians
    pub omega: f64,
    /// Longitude of ascending node in radians
    pub omega_big: f64,
    /// Mean anomaly in radians
    pub mean_anomaly: f64,
    /// Geometric albedo (Neptune's V-band value, 0.41)
    pub albedo: f64,
    /// Apparent V-band magnitude at current position
    pub v_magnitude: f64,
}

/// Apparent V magnitude of a planet of `mass_earth` at heliocentric distance
/// `helio_distance_au`, observed at opposition.
///
/// Delegates to `p9_core::analysis::photometry::planet_apparent_magnitude`.
pub fn brightness_at_position(mass_earth: f64, helio_distance_au: f64, albedo: f64) -> f64 {
    planet_apparent_magnitude(mass_earth, albedo, helio_distance_au)
}

/// Generate a reference population of synthetic Planet Nine objects.
pub fn generate_reference_population<R: Rng>(n: usize, rng: &mut R) -> Vec<ReferenceP9> {
    let posterior = mcmc_2021_posterior();
    generate_reference_population_with_posterior(n, &posterior, rng)
}

/// Generate a reference population using a specific posterior.
pub fn generate_reference_population_with_posterior<R: Rng>(
    n: usize,
    posterior: &P9Posterior,
    rng: &mut R,
) -> Vec<ReferenceP9> {
    let mut population = Vec::with_capacity(n);

    for _ in 0..n {
        let params = sample_from_posterior(posterior, rng);
        let albedo = ALBEDO_NEPTUNE;

        let helio_dist = heliocentric_distance(&params);
        let v_mag = brightness_at_position(params.mass_earth, helio_dist, albedo);

        population.push(ReferenceP9 {
            mass: params.mass_earth,
            a: params.a,
            e: params.e,
            i: params.i,
            omega: params.omega,
            omega_big: params.omega_big,
            mean_anomaly: params.mean_anomaly,
            albedo,
            v_magnitude: v_mag,
        });
    }

    population
}

/// Compute heliocentric distance from orbital elements at the given mean
/// anomaly, via the safeguarded p9-core Kepler solver (the previous local
/// Newton loop started at E₀ = M, which stalls for e → 1 near perihelion).
pub fn heliocentric_distance(params: &P9Params) -> f64 {
    let e = params.e;
    let ea = solve_kepler(e, params.mean_anomaly);
    params.a * (1.0 - e * ea.cos())
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;
    use rand::SeedableRng;

    #[test]
    fn test_brightness_perihelion() {
        let v = brightness_at_position(6.2, 300.0, ALBEDO_NEPTUNE);
        assert!(v > 15.0 && v < 30.0, "V magnitude at perihelion: {:.1}", v);
    }

    #[test]
    fn test_brightness_reflected_light_scaling() {
        // Reflected sunlight: doubling the distance dims by ~10 log10(2) ≈ 3 mag.
        let v_300 = brightness_at_position(6.2, 300.0, ALBEDO_NEPTUNE);
        let v_600 = brightness_at_position(6.2, 600.0, ALBEDO_NEPTUNE);
        let dm = v_600 - v_300;
        assert!((dm - 3.01).abs() < 0.05, "Δm for 2x distance = {dm:.2}");
    }

    #[test]
    fn test_brightness_mass_dependence() {
        let v_light = brightness_at_position(3.0, 300.0, ALBEDO_NEPTUNE);
        let v_heavy = brightness_at_position(10.0, 300.0, ALBEDO_NEPTUNE);
        assert!(
            v_light > v_heavy,
            "Heavier planet should be brighter: {:.1} vs {:.1}",
            v_light,
            v_heavy
        );
    }

    #[test]
    fn test_generate_population() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let pop = generate_reference_population(100, &mut rng);

        assert_eq!(pop.len(), 100);
        for obj in &pop {
            assert!(obj.mass > 0.0);
            assert!(obj.a > 0.0);
            assert!(obj.e > 0.0 && obj.e < 1.0);
            assert!((obj.albedo - ALBEDO_NEPTUNE).abs() < 1e-12);
            assert!(obj.v_magnitude.is_finite());
        }
    }

    #[test]
    fn test_heliocentric_distance_perihelion() {
        let params = P9Params {
            mass_earth: 6.2,
            a: 380.0,
            e: 0.21,
            i: 16.0 * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        };
        let r = heliocentric_distance(&params);
        let expected_q = 380.0 * (1.0 - 0.21);
        assert!(
            (r - expected_q).abs() < 0.01,
            "At M=0, should be at perihelion: r={:.2}, q={:.2}",
            r,
            expected_q
        );
    }

    #[test]
    fn test_heliocentric_distance_high_e_aphelion() {
        // The safeguarded solver must handle e = 0.95 at aphelion.
        let params = P9Params {
            mass_earth: 6.2,
            a: 380.0,
            e: 0.95,
            i: 0.0,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: std::f64::consts::PI,
        };
        let r = heliocentric_distance(&params);
        let expected = 380.0 * 1.95;
        assert!((r - expected).abs() < 0.01, "r = {r:.2}");
    }
}
