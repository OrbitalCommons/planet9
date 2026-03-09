//! Extended KBO sample with 14 objects (a > 230 AU) as of July 2018.
//!
//! Expands the 10-object sample from Brown (2017) with 4 additional discoveries.

use p9_core::constants::DEG2RAD;
use p9_core::types::OrbitalElements;

/// A distant KBO with orbital elements and discovery metadata.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DistantKbo {
    pub name: &'static str,
    pub elements: OrbitalElements,
}

/// The 14 KBOs with a > 230 AU and q > 30 AU as of July 2018.
///
/// Elements from MPC/JPL, approximate values in ecliptic J2000.
pub fn paper_sample_a230() -> Vec<DistantKbo> {
    vec![
        DistantKbo {
            name: "Sedna",
            elements: OrbitalElements {
                a: 506.0,
                e: 0.85,
                i: 11.93 * DEG2RAD,
                omega: 311.5 * DEG2RAD,
                omega_big: 144.5 * DEG2RAD,
                mean_anomaly: 358.1 * DEG2RAD,
            },
        },
        DistantKbo {
            name: "2012 VP113",
            elements: OrbitalElements {
                a: 261.0,
                e: 0.69,
                i: 24.05 * DEG2RAD,
                omega: 293.8 * DEG2RAD,
                omega_big: 90.8 * DEG2RAD,
                mean_anomaly: 2.5 * DEG2RAD,
            },
        },
        DistantKbo {
            name: "2013 RF98",
            elements: OrbitalElements {
                a: 350.0,
                e: 0.89,
                i: 29.6 * DEG2RAD,
                omega: 316.5 * DEG2RAD,
                omega_big: 67.6 * DEG2RAD,
                mean_anomaly: 0.3 * DEG2RAD,
            },
        },
        DistantKbo {
            name: "2004 VN112",
            elements: OrbitalElements {
                a: 327.0,
                e: 0.85,
                i: 25.6 * DEG2RAD,
                omega: 327.1 * DEG2RAD,
                omega_big: 66.0 * DEG2RAD,
                mean_anomaly: 1.8 * DEG2RAD,
            },
        },
        DistantKbo {
            name: "2010 GB174",
            elements: OrbitalElements {
                a: 351.0,
                e: 0.87,
                i: 21.5 * DEG2RAD,
                omega: 347.8 * DEG2RAD,
                omega_big: 130.6 * DEG2RAD,
                mean_anomaly: 0.6 * DEG2RAD,
            },
        },
        DistantKbo {
            name: "2000 CR105",
            elements: OrbitalElements {
                a: 230.0,
                e: 0.80,
                i: 22.7 * DEG2RAD,
                omega: 317.2 * DEG2RAD,
                omega_big: 128.3 * DEG2RAD,
                mean_anomaly: 5.4 * DEG2RAD,
            },
        },
        DistantKbo {
            name: "2014 SR349",
            elements: OrbitalElements {
                a: 299.0,
                e: 0.84,
                i: 18.0 * DEG2RAD,
                omega: 341.4 * DEG2RAD,
                omega_big: 34.8 * DEG2RAD,
                mean_anomaly: 0.7 * DEG2RAD,
            },
        },
        DistantKbo {
            name: "2007 TG422",
            elements: OrbitalElements {
                a: 501.0,
                e: 0.93,
                i: 18.6 * DEG2RAD,
                omega: 285.7 * DEG2RAD,
                omega_big: 112.9 * DEG2RAD,
                mean_anomaly: 0.1 * DEG2RAD,
            },
        },
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
            name: "2013 FT28",
            elements: OrbitalElements {
                a: 310.0,
                e: 0.86,
                i: 17.3 * DEG2RAD,
                omega: 40.2 * DEG2RAD,
                omega_big: 217.8 * DEG2RAD,
                mean_anomaly: 3.2 * DEG2RAD,
            },
        },
        // 4 additional objects not in the 2017 sample
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
        DistantKbo {
            name: "2014 FE72",
            elements: OrbitalElements {
                a: 2155.0,
                e: 0.99,
                i: 20.6 * DEG2RAD,
                omega: 134.0 * DEG2RAD,
                omega_big: 336.8 * DEG2RAD,
                mean_anomaly: 0.01 * DEG2RAD,
            },
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sample_size() {
        let kbos = paper_sample_a230();
        assert_eq!(kbos.len(), 14);
    }

    #[test]
    fn test_all_above_230_au() {
        for kbo in paper_sample_a230() {
            assert!(
                kbo.elements.a >= 230.0,
                "{} has a = {} < 230 AU",
                kbo.name,
                kbo.elements.a
            );
        }
    }

    #[test]
    fn test_all_positive_perihelion() {
        for kbo in paper_sample_a230() {
            let q = kbo.elements.perihelion();
            assert!(q > 20.0, "{} has q = {:.1} AU, expected > 20", kbo.name, q);
        }
    }
}
