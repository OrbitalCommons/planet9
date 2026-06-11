//! Inner Oort Cloud model.
//!
//! Generates a synthetic inner Oort cloud (IOC) population based on the
//! birth cluster environment described in Batygin & Brown (2021). The Sun's
//! birth cluster is modeled as a Plummer sphere with M ~ 1200 Msun and
//! characteristic radius ~ 0.35 pc. Stellar encounters within this cluster
//! lift planetesimals from the protoplanetary disk into IOC orbits with
//! semi-major axes in the range 2000-20000 AU.

use rand::Rng;
use rand_distr::{Distribution, Normal, Uniform};
use serde::{Deserialize, Serialize};

use p9_core::constants::*;
use p9_core::types::OrbitalElements;

/// Birth cluster parameters for IOC generation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OortCloudConfig {
    /// Total cluster mass in solar masses (nominal: 1200 Msun)
    pub cluster_mass_msun: f64,
    /// Plummer sphere characteristic radius in parsecs (nominal: 0.35 pc)
    pub cluster_radius_pc: f64,
    /// Minimum semi-major axis for IOC objects (AU)
    pub a_min: f64,
    /// Maximum semi-major axis for IOC objects (AU)
    pub a_max: f64,
    /// Minimum perihelion distance (AU) — objects must be detached
    pub q_min: f64,
    /// Maximum perihelion distance (AU)
    pub q_max: f64,
    /// Inclination dispersion (radians) — isotropic for IOC
    pub sigma_i: f64,
    /// Number of particles to generate
    pub n_particles: usize,
}

impl OortCloudConfig {
    /// Nominal configuration from Batygin & Brown (2021).
    pub fn nominal() -> Self {
        Self {
            cluster_mass_msun: 1200.0,
            cluster_radius_pc: 0.35,
            a_min: 2000.0,
            a_max: 20000.0,
            q_min: 40.0,
            q_max: 500.0,
            sigma_i: 40.0 * DEG2RAD,
            n_particles: 500,
        }
    }

    /// Compact cluster variant (denser environment, more IOC production).
    pub fn compact_cluster() -> Self {
        Self {
            cluster_mass_msun: 2000.0,
            cluster_radius_pc: 0.2,
            a_min: 2000.0,
            a_max: 20000.0,
            q_min: 40.0,
            q_max: 500.0,
            sigma_i: 45.0 * DEG2RAD,
            n_particles: 500,
        }
    }

    /// Cluster velocity dispersion in AU/day (estimated from virial theorem).
    /// sigma_v ~ sqrt(G * M / r)
    pub fn velocity_dispersion_au_day(&self) -> f64 {
        let r_au = self.cluster_radius_pc * PC_AU;
        (GM_SUN * self.cluster_mass_msun / r_au).sqrt()
    }

    /// Typical encounter timescale in days.
    /// t_enc ~ r / sigma_v
    pub fn encounter_timescale_days(&self) -> f64 {
        let r_au = self.cluster_radius_pc * PC_AU;
        let sigma_v = self.velocity_dispersion_au_day();
        r_au / sigma_v
    }
}

/// Generate a synthetic inner Oort cloud population.
///
/// Creates IOC objects with semi-major axes uniformly distributed in
/// [a_min, a_max], perihelion distances in [q_min, q_max], and
/// isotropic angular distributions appropriate for a population
/// emplaced by stellar encounters in the birth cluster.
pub fn generate_ioc_population<R: Rng>(
    config: &OortCloudConfig,
    rng: &mut R,
) -> Vec<OrbitalElements> {
    let a_dist = Uniform::new(config.a_min, config.a_max);
    let q_dist = Uniform::new(config.q_min, config.q_max);
    let angle_dist = Uniform::new(0.0, TWO_PI);
    let incl_normal = Normal::new(0.0, config.sigma_i).unwrap();

    let mut population = Vec::with_capacity(config.n_particles);

    while population.len() < config.n_particles {
        let a = a_dist.sample(rng);
        let q = q_dist.sample(rng);
        let e = 1.0 - q / a;

        if e < 0.0 || e >= 1.0 {
            continue;
        }

        // IOC inclinations are broadly distributed but not fully isotropic;
        // use a wide Gaussian centered on the ecliptic
        let i = incl_normal.sample(rng).abs().min(std::f64::consts::PI);

        let omega = angle_dist.sample(rng);
        let omega_big = angle_dist.sample(rng);
        let mean_anomaly = angle_dist.sample(rng);

        population.push(OrbitalElements {
            a,
            e,
            i,
            omega_big,
            omega,
            mean_anomaly,
        });
    }

    population
}

/// Fraction of scattered planetesimals trapped into the inner Oort cloud
/// during the Sun's birth-cluster phase.
///
/// Documented literature assumption, not a computed quantity: cluster N-body
/// simulations (Brasser, Duncan & Levison 2006; cited by Batygin & Brown
/// 2021 for their IOC source population) find a trapping efficiency of a few
/// per cent for ~10³ M☉ embedded clusters. Computing it would require the
/// full cluster encounter simulation, which is outside this crate's scope.
/// (Replaces an invented 0.05·√(ρ/ρ₀) scaling whose only test was circular.)
pub const IOC_TRAPPING_FRACTION_ASSUMED: f64 = 0.05;

/// Generate a distant scattered-disk comparison population: a ∈ (250, 550)
/// AU (the paper's distant-belt membership range), q ∈ (30, 50) AU, moderate
/// inclinations, uniform angles. Used as the reduced-scale simulated
/// population for the f_ϖ and semi-major-axis comparisons (previously the
/// comparison histogram was hard-coded counts).
pub fn generate_scattered_disk<R: Rng>(n: usize, rng: &mut R) -> Vec<OrbitalElements> {
    let a_dist = Uniform::new(250.0, 550.0);
    let q_dist = Uniform::new(30.0, 50.0);
    let angle_dist = Uniform::new(0.0, TWO_PI);
    let incl_normal = Normal::new(0.0, 15.0 * DEG2RAD).unwrap();

    let mut population = Vec::with_capacity(n);
    while population.len() < n {
        let a = a_dist.sample(rng);
        let q = q_dist.sample(rng);
        let e = 1.0 - q / a;
        if !(0.0..1.0).contains(&e) {
            continue;
        }
        population.push(OrbitalElements {
            a,
            e,
            i: incl_normal.sample(rng).abs().min(std::f64::consts::PI),
            omega: angle_dist.sample(rng),
            omega_big: angle_dist.sample(rng),
            mean_anomaly: angle_dist.sample(rng),
        });
    }
    population
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_nominal_config() {
        let config = OortCloudConfig::nominal();
        assert!((config.cluster_mass_msun - 1200.0).abs() < 1e-10);
        assert!((config.cluster_radius_pc - 0.35).abs() < 1e-10);
        assert!(config.a_min < config.a_max);
    }

    #[test]
    fn test_velocity_dispersion() {
        let config = OortCloudConfig::nominal();
        let sigma = config.velocity_dispersion_au_day();
        assert!(sigma > 0.0);
        // Convert to km/s: should be ~1-5 km/s for a modest cluster
        let sigma_kms = sigma / KMS_TO_AUDAY;
        assert!(sigma_kms > 0.1, "sigma_v = {sigma_kms} km/s too low");
        assert!(sigma_kms < 20.0, "sigma_v = {sigma_kms} km/s too high");
    }

    #[test]
    fn test_generate_ioc_population_count() {
        let mut rng = StdRng::seed_from_u64(42);
        let config = OortCloudConfig {
            n_particles: 100,
            ..OortCloudConfig::nominal()
        };
        let pop = generate_ioc_population(&config, &mut rng);
        assert_eq!(pop.len(), 100);
    }

    #[test]
    fn test_ioc_orbital_ranges() {
        let mut rng = StdRng::seed_from_u64(123);
        let config = OortCloudConfig::nominal();
        let pop = generate_ioc_population(&config, &mut rng);

        for elem in &pop {
            assert!(
                elem.a >= config.a_min && elem.a <= config.a_max,
                "a = {} out of range",
                elem.a
            );
            let q = elem.perihelion();
            assert!(
                q >= config.q_min * 0.99 && q <= config.q_max * 1.01,
                "q = {} out of range [{}, {}]",
                q,
                config.q_min,
                config.q_max
            );
            assert!(elem.e >= 0.0 && elem.e < 1.0, "e = {} invalid", elem.e);
            assert!(elem.i >= 0.0, "i must be non-negative");
        }
    }

    #[test]
    fn test_scattered_disk_population_ranges() {
        let mut rng = StdRng::seed_from_u64(3);
        let pop = generate_scattered_disk(200, &mut rng);
        assert_eq!(pop.len(), 200);
        for e in &pop {
            assert!(e.a >= 250.0 && e.a <= 550.0);
            let q = e.perihelion();
            assert!(q >= 29.9 && q <= 50.1, "q = {q}");
            assert!(e.e > 0.0 && e.e < 1.0);
        }
    }
}
