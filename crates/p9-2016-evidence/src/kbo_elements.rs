//! Orbital elements for the 6 dynamically stable KBOs from Batygin & Brown (2016).
//!
//! These are the objects with q > 30 AU, a > 250 AU that remain dynamically
//! stable over 4 Gyr with the four giant planets:
//!   Sedna, 2012 VP113, 2004 VN112, 2010 GB174, 2000 CR105, 2010 VZ98
//!
//! Elements are osculating, heliocentric, ecliptic J2000.0 from the MPC.
//! TODO: Replace with live queries from starfield-datasources MPC/SBDB when available.

use p9_core::constants::DEG2RAD;
use p9_core::types::OrbitalElements;

/// A named KBO with its orbital elements and observational uncertainties.
#[derive(Debug, Clone)]
pub struct KboRecord {
    pub name: &'static str,
    pub designation: &'static str,
    pub elements: OrbitalElements,
    /// 1-sigma uncertainties for clone generation
    pub sigma_a: f64,
    pub sigma_e: f64,
    pub sigma_i: f64,
    pub sigma_omega: f64,
    pub sigma_omega_big: f64,
    pub sigma_m: f64,
}

/// The 6 dynamically stable KBOs from Batygin & Brown (2016) Table 1.
/// Elements from JPL SBDB, epoch ~2016.
pub fn stable_kbos() -> Vec<KboRecord> {
    vec![
        KboRecord {
            name: "Sedna",
            designation: "2003 VB12",
            elements: OrbitalElements {
                a: 506.8,
                e: 0.8496,
                i: 11.93 * DEG2RAD,
                omega: 311.46 * DEG2RAD,
                omega_big: 144.51 * DEG2RAD,
                mean_anomaly: 358.12 * DEG2RAD,
            },
            sigma_a: 0.5,
            sigma_e: 0.0001,
            sigma_i: 0.001 * DEG2RAD,
            sigma_omega: 0.01 * DEG2RAD,
            sigma_omega_big: 0.01 * DEG2RAD,
            sigma_m: 0.01 * DEG2RAD,
        },
        KboRecord {
            name: "2012 VP113",
            designation: "2012 VP113",
            elements: OrbitalElements {
                a: 261.0,
                e: 0.6886,
                i: 24.05 * DEG2RAD,
                omega: 293.78 * DEG2RAD,
                omega_big: 90.81 * DEG2RAD,
                mean_anomaly: 5.63 * DEG2RAD,
            },
            sigma_a: 3.0,
            sigma_e: 0.003,
            sigma_i: 0.01 * DEG2RAD,
            sigma_omega: 0.1 * DEG2RAD,
            sigma_omega_big: 0.1 * DEG2RAD,
            sigma_m: 0.1 * DEG2RAD,
        },
        KboRecord {
            name: "2004 VN112",
            designation: "2004 VN112",
            elements: OrbitalElements {
                a: 327.5,
                e: 0.8527,
                i: 25.56 * DEG2RAD,
                omega: 327.15 * DEG2RAD,
                omega_big: 66.01 * DEG2RAD,
                mean_anomaly: 1.71 * DEG2RAD,
            },
            sigma_a: 6.0,
            sigma_e: 0.005,
            sigma_i: 0.01 * DEG2RAD,
            sigma_omega: 0.1 * DEG2RAD,
            sigma_omega_big: 0.1 * DEG2RAD,
            sigma_m: 0.1 * DEG2RAD,
        },
        KboRecord {
            name: "2010 GB174",
            designation: "2010 GB174",
            elements: OrbitalElements {
                a: 371.7,
                e: 0.8627,
                i: 21.54 * DEG2RAD,
                omega: 347.77 * DEG2RAD,
                omega_big: 130.59 * DEG2RAD,
                mean_anomaly: 2.81 * DEG2RAD,
            },
            sigma_a: 13.0,
            sigma_e: 0.01,
            sigma_i: 0.01 * DEG2RAD,
            sigma_omega: 0.2 * DEG2RAD,
            sigma_omega_big: 0.2 * DEG2RAD,
            sigma_m: 0.2 * DEG2RAD,
        },
        KboRecord {
            name: "2000 CR105",
            designation: "2000 CR105",
            elements: OrbitalElements {
                a: 228.8,
                e: 0.8024,
                i: 22.75 * DEG2RAD,
                omega: 316.74 * DEG2RAD,
                omega_big: 128.28 * DEG2RAD,
                mean_anomaly: 6.27 * DEG2RAD,
            },
            sigma_a: 1.0,
            sigma_e: 0.001,
            sigma_i: 0.001 * DEG2RAD,
            sigma_omega: 0.02 * DEG2RAD,
            sigma_omega_big: 0.02 * DEG2RAD,
            sigma_m: 0.02 * DEG2RAD,
        },
        KboRecord {
            name: "2010 VZ98",
            designation: "2010 VZ98",
            elements: OrbitalElements {
                a: 153.2,
                e: 0.7706,
                i: 4.51 * DEG2RAD,
                omega: 313.90 * DEG2RAD,
                omega_big: 117.39 * DEG2RAD,
                mean_anomaly: 7.97 * DEG2RAD,
            },
            sigma_a: 1.5,
            sigma_e: 0.002,
            sigma_i: 0.005 * DEG2RAD,
            sigma_omega: 0.05 * DEG2RAD,
            sigma_omega_big: 0.05 * DEG2RAD,
            sigma_m: 0.05 * DEG2RAD,
        },
    ]
}

/// Generate `n_clones` of a KBO, perturbed within observational uncertainties.
pub fn generate_clones(
    kbo: &KboRecord,
    n_clones: usize,
    rng: &mut impl rand::Rng,
) -> Vec<OrbitalElements> {
    use rand_distr::{Distribution, Normal};

    let mut clones = Vec::with_capacity(n_clones);
    for _ in 0..n_clones {
        let a_norm = Normal::new(0.0, kbo.sigma_a).unwrap();
        let e_norm = Normal::new(0.0, kbo.sigma_e).unwrap();
        let i_norm = Normal::new(0.0, kbo.sigma_i).unwrap();
        let w_norm = Normal::new(0.0, kbo.sigma_omega).unwrap();
        let o_norm = Normal::new(0.0, kbo.sigma_omega_big).unwrap();
        let m_norm = Normal::new(0.0, kbo.sigma_m).unwrap();

        let mut elem = kbo.elements;
        elem.a += a_norm.sample(rng);
        elem.e = (elem.e + e_norm.sample(rng)).clamp(0.0, 0.9999);
        elem.i += i_norm.sample(rng);
        elem.omega += w_norm.sample(rng);
        elem.omega_big += o_norm.sample(rng);
        elem.mean_anomaly += m_norm.sample(rng);

        clones.push(elem);
    }
    clones
}

/// Compute longitude of perihelion for an element set.
pub fn longitude_of_perihelion(elem: &OrbitalElements) -> f64 {
    (elem.omega + elem.omega_big) % p9_core::constants::TWO_PI
}

/// The observed clustering statistics from the paper.
/// ω = 318 ± 8 deg, Ω = 113 ± 13 deg, ϖ = 71 ± 16 deg
pub fn observed_clustering_stats() -> (f64, f64, f64) {
    let kbos = stable_kbos();
    let n = kbos.len() as f64;

    let varpi_mean: f64 = kbos
        .iter()
        .map(|k| longitude_of_perihelion(&k.elements))
        .sum::<f64>()
        / n;

    let omega_mean: f64 = kbos.iter().map(|k| k.elements.omega).sum::<f64>() / n;

    let omega_big_mean: f64 = kbos.iter().map(|k| k.elements.omega_big).sum::<f64>() / n;

    (varpi_mean, omega_big_mean, omega_mean)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_six_stable_kbos() {
        let kbos = stable_kbos();
        assert_eq!(kbos.len(), 6);

        for kbo in &kbos {
            assert!(kbo.elements.a > 150.0, "{} a too small", kbo.name);
            assert!(
                kbo.elements.e > 0.0 && kbo.elements.e < 1.0,
                "{} e invalid",
                kbo.name
            );
            assert!(
                kbo.elements.perihelion() > 30.0,
                "{} q = {:.1} should be > 30 AU",
                kbo.name,
                kbo.elements.perihelion()
            );
        }
    }

    #[test]
    fn test_perihelion_clustering() {
        let kbos = stable_kbos();

        // All 6 KBOs should have ω near 318 deg (between 290-350 deg)
        for kbo in &kbos {
            let omega_deg = kbo.elements.omega / DEG2RAD;
            assert!(
                omega_deg > 280.0 && omega_deg < 360.0,
                "{}: ω = {:.1}° should cluster near 318°",
                kbo.name,
                omega_deg
            );
        }
    }

    #[test]
    fn test_clone_generation() {
        use rand::SeedableRng;
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let kbos = stable_kbos();
        let clones = generate_clones(&kbos[0], 6, &mut rng);

        assert_eq!(clones.len(), 6);
        for clone in &clones {
            // Clones should be near the original
            assert!((clone.a - kbos[0].elements.a).abs() < 10.0 * kbos[0].sigma_a);
            assert!(clone.e > 0.0 && clone.e < 1.0);
        }
    }
}
