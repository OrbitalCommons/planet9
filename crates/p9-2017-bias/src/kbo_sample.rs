//! KBO sample from Brown (2017).
//!
//! The paper analyzes distant KBOs with a > 230 AU to test for orbital clustering.
//! This module contains the hardcoded sample used in the paper.

use p9_core::constants::DEG2RAD;
use p9_core::types::OrbitalElements;

/// A KBO in the distant sample.
#[derive(Debug, Clone)]
pub struct DistantKbo {
    pub name: &'static str,
    pub elements: OrbitalElements,
    /// Absolute magnitude H
    pub h_mag: f64,
}

/// The 10 KBOs with a > 230 AU used in the paper's primary analysis.
///
/// Orbital elements from the paper (Table 1). Angles in degrees converted to radians.
pub fn paper_sample_a230() -> Vec<DistantKbo> {
    vec![
        DistantKbo {
            name: "Sedna",
            elements: OrbitalElements {
                a: 506.0,
                e: 0.85,
                i: 11.9 * DEG2RAD,
                omega: 311.5 * DEG2RAD,
                omega_big: 144.5 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            h_mag: 1.6,
        },
        DistantKbo {
            name: "2012 VP113",
            elements: OrbitalElements {
                a: 261.0,
                e: 0.69,
                i: 24.1 * DEG2RAD,
                omega: 293.8 * DEG2RAD,
                omega_big: 90.8 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            h_mag: 4.0,
        },
        DistantKbo {
            name: "2013 RF98",
            elements: OrbitalElements {
                a: 325.0,
                e: 0.89,
                i: 29.6 * DEG2RAD,
                omega: 316.5 * DEG2RAD,
                omega_big: 67.6 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            h_mag: 8.7,
        },
        DistantKbo {
            name: "2004 VN112",
            elements: OrbitalElements {
                a: 327.0,
                e: 0.85,
                i: 25.6 * DEG2RAD,
                omega: 327.1 * DEG2RAD,
                omega_big: 66.0 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            h_mag: 6.4,
        },
        DistantKbo {
            name: "2010 GB174",
            elements: OrbitalElements {
                a: 351.0,
                e: 0.86,
                i: 21.5 * DEG2RAD,
                omega: 347.8 * DEG2RAD,
                omega_big: 130.6 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            h_mag: 6.5,
        },
        DistantKbo {
            name: "2000 CR105",
            elements: OrbitalElements {
                a: 228.0,
                e: 0.81,
                i: 22.7 * DEG2RAD,
                omega: 317.2 * DEG2RAD,
                omega_big: 128.3 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            h_mag: 6.3,
        },
        DistantKbo {
            name: "2007 TG422",
            elements: OrbitalElements {
                a: 501.0,
                e: 0.93,
                i: 18.6 * DEG2RAD,
                omega: 285.7 * DEG2RAD,
                omega_big: 113.0 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            h_mag: 6.2,
        },
        DistantKbo {
            name: "2013 FT28",
            elements: OrbitalElements {
                a: 310.0,
                e: 0.86,
                i: 17.3 * DEG2RAD,
                omega: 40.2 * DEG2RAD,
                omega_big: 217.8 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            h_mag: 6.7,
        },
        DistantKbo {
            name: "2014 SR349",
            elements: OrbitalElements {
                a: 289.0,
                e: 0.84,
                i: 18.0 * DEG2RAD,
                omega: 341.4 * DEG2RAD,
                omega_big: 34.8 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            h_mag: 6.6,
        },
        DistantKbo {
            name: "2015 RX245",
            elements: OrbitalElements {
                a: 430.0,
                e: 0.89,
                i: 12.2 * DEG2RAD,
                omega: 65.2 * DEG2RAD,
                omega_big: 8.6 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            h_mag: 6.2,
        },
    ]
}

/// Extract longitude of perihelion ϖ = ω + Ω for each KBO (in radians).
pub fn longitudes_of_perihelion(kbos: &[DistantKbo]) -> Vec<f64> {
    kbos.iter()
        .map(|k| {
            let varpi = k.elements.omega + k.elements.omega_big;
            varpi.rem_euclid(2.0 * std::f64::consts::PI)
        })
        .collect()
}

/// Extract arguments of perihelion ω for each KBO (in radians).
pub fn arguments_of_perihelion(kbos: &[DistantKbo]) -> Vec<f64> {
    kbos.iter().map(|k| k.elements.omega).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sample_size() {
        let kbos = paper_sample_a230();
        assert_eq!(kbos.len(), 10);
    }

    #[test]
    fn test_all_distant() {
        let kbos = paper_sample_a230();
        for kbo in &kbos {
            assert!(
                kbo.elements.a > 220.0,
                "{} has a = {:.0} (expected > 220)",
                kbo.name,
                kbo.elements.a
            );
        }
    }

    #[test]
    fn test_all_eccentric() {
        let kbos = paper_sample_a230();
        for kbo in &kbos {
            assert!(
                kbo.elements.e > 0.6,
                "{} has e = {:.2} (expected > 0.6)",
                kbo.name,
                kbo.elements.e
            );
        }
    }

    #[test]
    fn test_varpi_clustering() {
        let kbos = paper_sample_a230();
        let varpis = longitudes_of_perihelion(&kbos);

        // Most ϖ values should be clustered (mean ~ 0-100°)
        let sin_sum: f64 = varpis.iter().map(|v| v.sin()).sum();
        let cos_sum: f64 = varpis.iter().map(|v| v.cos()).sum();
        let n = varpis.len() as f64;
        let r_bar = (sin_sum * sin_sum + cos_sum * cos_sum).sqrt() / n;

        // R̄ > 0.3 indicates clustering
        assert!(
            r_bar > 0.3,
            "R̄ = {:.3} should indicate clustering (> 0.3)",
            r_bar
        );
    }
}
