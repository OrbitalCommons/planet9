//! Orbital elements for the 6 dynamically stable KBOs from Batygin & Brown
//! (2016).
//!
//! These are the objects with q > 30 AU, a > 150 AU the paper argues remain
//! dynamically stable over 4 Gyr with the four giant planets:
//!   Sedna, 2012 VP113, 2004 VN112, 2010 GB174, 2000 CR105, 2010 VZ98
//!
//! Elements are osculating, heliocentric, ecliptic J2000.0 from JPL SBDB
//! (epoch ~2016). The struct is data-only: published orbit-fit covariances
//! are not transcribed here.
//!
//! Overlap with [`crate::data::etno`]: this is a *distinct* table, not a
//! subset of `BROWN_2017_SAMPLE`. It pins the 6-object stability sample from
//! the 2016 evidence paper at higher angular/eccentricity precision, carries
//! per-object mean anomalies and SBDB designations, and includes 2010 VZ98
//! (a ≈ 153 AU), which is below the a ≳ 230 AU cut of the Brown (2017)
//! clustering sample. The shared objects also differ numerically (e.g. Sedna
//! a = 506.8 here vs 506.0 in the ETNO table). The two tables are kept
//! separate on purpose — do not merge them without re-pinning every
//! downstream regression.

use crate::constants::*;
use crate::types::OrbitalElements;

/// A named KBO with its orbital elements.
#[derive(Debug, Clone)]
pub struct KboRecord {
    pub name: &'static str,
    pub designation: &'static str,
    pub elements: OrbitalElements,
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
        },
    ]
}

/// Compute longitude of perihelion for an element set.
pub fn longitude_of_perihelion(elem: &OrbitalElements) -> f64 {
    (elem.omega + elem.omega_big).rem_euclid(TWO_PI)
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
}
