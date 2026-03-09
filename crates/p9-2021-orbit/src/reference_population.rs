//! Synthetic reference population generation for Planet Nine.
//!
//! Generates a large population of synthetic Planet Nine realizations
//! drawn from the MCMC posterior, computing observable properties
//! (apparent magnitude, sky position) for each.

use rand::Rng;
use serde::{Deserialize, Serialize};

use p9_core::types::P9Params;

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
    /// Geometric albedo (assumed)
    pub albedo: f64,
    /// Apparent V-band magnitude at current position
    pub v_magnitude: f64,
}

/// Estimate V-band apparent magnitude of Planet Nine.
///
/// Uses the IAU standard absolute magnitude H = -2.5*log10(p) - 5*log10(D/1329)
/// where p is geometric albedo and D is diameter in km.
///
/// For distant objects at opposition (delta ~ r):
///   V = H + 5 * log10(r^2) = H + 10 * log10(r)
///
/// Planet radius estimated from mass using R ~ (M/M_earth)^0.55 * R_earth
/// (appropriate for super-Earth / mini-Neptune range).
pub fn brightness_at_position(mass_earth: f64, helio_distance_au: f64, albedo: f64) -> f64 {
    let r_earth_km = 6371.0;
    let radius_km = r_earth_km * mass_earth.powf(0.55);
    let diameter_km = 2.0 * radius_km;

    // IAU standard absolute magnitude: H = -2.5 * log10(p) - 5 * log10(D_km / 1329)
    let h_abs = -2.5 * albedo.log10() - 5.0 * (diameter_km / 1329.0).log10();

    // Apparent magnitude: V = H + 5 * log10(r * delta)
    // For distant objects at opposition, delta ~ r
    h_abs + 5.0 * (helio_distance_au * helio_distance_au).log10()
}

/// Generate a reference population of synthetic Planet Nine objects.
///
/// Each realization is drawn from the MCMC posterior. The heliocentric
/// distance is computed from the orbital elements and a random true anomaly
/// (uniformly distributed in mean anomaly), and the apparent magnitude
/// is estimated assuming a geometric albedo of 0.5.
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
        let albedo = 0.3 + rng.gen::<f64>() * 0.4;

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

/// Compute heliocentric distance from orbital elements at the given mean anomaly.
///
/// Solves Kepler's equation to get the true anomaly, then computes
/// r = a(1-e^2) / (1 + e*cos(nu)).
fn heliocentric_distance(params: &P9Params) -> f64 {
    let e = params.e;
    let m = params.mean_anomaly;

    let ea = solve_kepler(e, m);

    let nu = 2.0 * ((1.0 + e).sqrt() * (ea / 2.0).sin()).atan2((1.0 - e).sqrt() * (ea / 2.0).cos());

    let p = params.a * (1.0 - e * e);
    p / (1.0 + e * nu.cos())
}

/// Solve Kepler's equation M = E - e*sin(E) via Newton-Raphson.
fn solve_kepler(e: f64, m: f64) -> f64 {
    let mut ea = m;
    for _ in 0..50 {
        let delta = (ea - e * ea.sin() - m) / (1.0 - e * ea.cos());
        ea -= delta;
        if delta.abs() < 1e-15 {
            break;
        }
    }
    ea
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;
    use rand::SeedableRng;

    #[test]
    fn test_brightness_perihelion() {
        let v = brightness_at_position(6.2, 300.0, 0.5);
        assert!(v > 15.0 && v < 30.0, "V magnitude at perihelion: {:.1}", v);
    }

    #[test]
    fn test_brightness_aphelion() {
        let v_peri = brightness_at_position(6.2, 300.0, 0.5);
        let v_apo = brightness_at_position(6.2, 460.0, 0.5);
        assert!(
            v_apo > v_peri,
            "Should be fainter at aphelion: {:.1} vs {:.1}",
            v_apo,
            v_peri
        );
    }

    #[test]
    fn test_brightness_mass_dependence() {
        let v_light = brightness_at_position(3.0, 300.0, 0.5);
        let v_heavy = brightness_at_position(10.0, 300.0, 0.5);
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
            assert!(obj.albedo > 0.0 && obj.albedo < 1.0);
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
}
