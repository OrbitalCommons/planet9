//! Known highly inclined TNOs from the paper.
//!
//! These are the observed objects that the simulation aims to reproduce:
//! - Drac (2008 KV42): retrograde centaur, i ≈ 103°
//! - Niku (2011 KT19): high-i TNO, i ≈ 110°
//! - 2016 NM56: retrograde, i ≈ 144°

use p9_core::constants::DEG2RAD;
use p9_core::types::OrbitalElements;

/// A known high-inclination TNO for comparison.
#[derive(Debug, Clone)]
pub struct KnownTno {
    pub name: &'static str,
    pub designation: &'static str,
    pub elements: OrbitalElements,
}

/// Return the three key high-inclination TNOs discussed in the paper.
pub fn paper_tnos() -> Vec<KnownTno> {
    vec![
        KnownTno {
            name: "Drac",
            designation: "2008 KV42",
            elements: OrbitalElements {
                a: 41.4,
                e: 0.49,
                i: 103.4 * DEG2RAD,
                omega: 262.8 * DEG2RAD,
                omega_big: 260.9 * DEG2RAD,
                mean_anomaly: 0.0,
            },
        },
        KnownTno {
            name: "Niku",
            designation: "2011 KT19",
            elements: OrbitalElements {
                a: 35.6,
                e: 0.33,
                i: 110.1 * DEG2RAD,
                omega: 83.4 * DEG2RAD,
                omega_big: 243.8 * DEG2RAD,
                mean_anomaly: 0.0,
            },
        },
        KnownTno {
            name: "2016 NM56",
            designation: "2016 NM56",
            elements: OrbitalElements {
                a: 74.0,
                e: 0.88,
                i: 144.0 * DEG2RAD,
                omega: 130.0 * DEG2RAD,
                omega_big: 250.0 * DEG2RAD,
                mean_anomaly: 0.0,
            },
        },
    ]
}

/// Additional high-inclination centaurs and TNOs for extended comparison.
pub fn extended_high_i_objects() -> Vec<KnownTno> {
    let mut objects = paper_tnos();
    objects.push(KnownTno {
        name: "2010 WG9",
        designation: "2010 WG9",
        elements: OrbitalElements {
            a: 53.7,
            e: 0.65,
            i: 70.3 * DEG2RAD,
            omega: 17.9 * DEG2RAD,
            omega_big: 92.1 * DEG2RAD,
            mean_anomaly: 0.0,
        },
    });
    objects.push(KnownTno {
        name: "2013 LA2",
        designation: "2013 LA2",
        elements: OrbitalElements {
            a: 43.4,
            e: 0.77,
            i: 152.0 * DEG2RAD,
            omega: 23.0 * DEG2RAD,
            omega_big: 66.0 * DEG2RAD,
            mean_anomaly: 0.0,
        },
    });
    objects
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_paper_tnos() {
        let tnos = paper_tnos();
        assert_eq!(tnos.len(), 3);

        for tno in &tnos {
            // All should have i > 90° (retrograde or near-retrograde)
            assert!(
                tno.elements.i > 90.0 * DEG2RAD,
                "{} has i = {:.1}° (expected > 90°)",
                tno.name,
                tno.elements.i / DEG2RAD
            );
            // All should have a < 100 AU (decoupled from P9)
            assert!(
                tno.elements.a < 100.0,
                "{} has a = {:.1} AU (expected < 100)",
                tno.name,
                tno.elements.a
            );
        }
    }

    #[test]
    fn test_drac_is_retrograde() {
        let tnos = paper_tnos();
        let drac = &tnos[0];
        assert!(drac.elements.i > 90.0 * DEG2RAD);
        assert!((drac.elements.a - 41.4).abs() < 1.0);
    }
}
