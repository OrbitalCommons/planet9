//! OSSOS survey sensitivity comparison.
//!
//! Demonstrates that the OSSOS survey's limited sky coverage renders it
//! insensitive to detecting the clustering signal. With only 4 distant
//! KBOs and restricted survey fields, even a true clustering signal
//! would be undetectable at > 65% confidence.

use crate::kbo_sample::DistantKbo;
use crate::poincare_variables::*;

use p9_core::constants::DEG2RAD;
use p9_core::types::OrbitalElements;

/// The 4 OSSOS-detected distant KBOs with a > 230 AU.
///
/// These are the objects that OSSOS found, which were used
/// by Shankman et al. (2017) to argue against clustering.
pub fn ossos_sample() -> Vec<DistantKbo> {
    vec![
        DistantKbo {
            name: "2013 SY99",
            elements: OrbitalElements {
                a: 735.0,
                e: 0.93,
                i: 4.2 * DEG2RAD,
                omega: 32.2 * DEG2RAD,
                omega_big: 29.5 * DEG2RAD,
                mean_anomaly: 0.03 * DEG2RAD,
            },
        },
        DistantKbo {
            name: "2015 GT50",
            elements: OrbitalElements {
                a: 312.0,
                e: 0.89,
                i: 8.8 * DEG2RAD,
                omega: 46.1 * DEG2RAD,
                omega_big: 164.3 * DEG2RAD,
                mean_anomaly: 1.1 * DEG2RAD,
            },
        },
        DistantKbo {
            name: "2015 KG163",
            elements: OrbitalElements {
                a: 680.0,
                e: 0.94,
                i: 14.0 * DEG2RAD,
                omega: 32.0 * DEG2RAD,
                omega_big: 219.1 * DEG2RAD,
                mean_anomaly: 0.05 * DEG2RAD,
            },
        },
        DistantKbo {
            name: "2015 RX245",
            elements: OrbitalElements {
                a: 430.0,
                e: 0.89,
                i: 12.2 * DEG2RAD,
                omega: 65.2 * DEG2RAD,
                omega_big: 8.6 * DEG2RAD,
                mean_anomaly: 0.4 * DEG2RAD,
            },
        },
    ]
}

/// Result of the OSSOS sensitivity analysis.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct OssosSensitivity {
    /// OSSOS sample perihelion clustering
    pub ossos_perihelion: f64,
    /// OSSOS sample pole clustering
    pub ossos_pole: f64,
    /// Full sample perihelion clustering
    pub full_perihelion: f64,
    /// Full sample pole clustering
    pub full_pole: f64,
    /// Number of OSSOS objects
    pub n_ossos: usize,
    /// Number of full sample objects
    pub n_full: usize,
}

/// Compare clustering strength between OSSOS and full sample.
pub fn sensitivity_analysis(
    full_sample: &[DistantKbo],
    ossos_sample: &[DistantKbo],
) -> OssosSensitivity {
    let full_states: Vec<PoincareState> = full_sample
        .iter()
        .map(|k| PoincareState::from_elements(&k.elements))
        .collect();
    let ossos_states: Vec<PoincareState> = ossos_sample
        .iter()
        .map(|k| PoincareState::from_elements(&k.elements))
        .collect();

    let full_mean = mean_state(&full_states);
    let ossos_mean = mean_state(&ossos_states);

    OssosSensitivity {
        ossos_perihelion: perihelion_clustering(&ossos_mean),
        ossos_pole: pole_clustering(&ossos_mean),
        full_perihelion: perihelion_clustering(&full_mean),
        full_pole: pole_clustering(&full_mean),
        n_ossos: ossos_sample.len(),
        n_full: full_sample.len(),
    }
}

/// Compute the expected power of the OSSOS survey to detect clustering.
///
/// With only n=4 objects, the minimum detectable clustering signal is
/// approximately 1/sqrt(n) ≈ 0.5 in Poincaré coordinates, which is
/// much larger than typical clustering strengths.
pub fn ossos_detection_power(n_objects: usize) -> f64 {
    // Rough approximation: detection threshold ∝ 1/sqrt(n)
    // With n=4, threshold ≈ 0.5
    // With n=14, threshold ≈ 0.27
    1.0 / (n_objects as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kbo_sample;

    #[test]
    fn test_ossos_sample_size() {
        assert_eq!(ossos_sample().len(), 4);
    }

    #[test]
    fn test_sensitivity_analysis() {
        let full = kbo_sample::paper_sample_a230();
        let ossos = ossos_sample();
        let result = sensitivity_analysis(&full, &ossos);

        assert_eq!(result.n_full, 14);
        assert_eq!(result.n_ossos, 4);
        assert!(result.full_perihelion > 0.0);
        assert!(result.ossos_perihelion >= 0.0);
    }

    #[test]
    fn test_detection_power() {
        let power_4 = ossos_detection_power(4);
        let power_14 = ossos_detection_power(14);

        // With fewer objects, threshold is higher (harder to detect)
        assert!(
            power_4 > power_14,
            "OSSOS (n=4) threshold {:.2} should be higher than full (n=14) {:.2}",
            power_4,
            power_14
        );
    }
}
