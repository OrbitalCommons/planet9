//! Extreme trans-Neptunian object (ETNO) observational table.
//!
//! Single vetted source for the workspace. Four crates previously carried
//! their own copies of this sample; the p9-2018-resonance copy had scrambled
//! semi-major axes (e.g. 2013 RF98 at 780 AU instead of ~325 AU) that
//! displaced every implied-a₉ resonance peak by tens of AU.
//!
//! Provenance: Brown (2017) Table 1 (the a > 230 AU clustering sample, with
//! 2000 CR105 at a ≈ 228 AU included as in the paper), cross-checked against
//! the JPL Small-Body Database (osculating elements; angles rounded to
//! 0.1 deg, a to the AU). Mean anomalies are not part of the clustering
//! analyses and are set to 0.

use crate::constants::{DEG2RAD, TWO_PI};
use crate::types::OrbitalElements;

/// An ETNO in the vetted sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Etno {
    pub name: &'static str,
    /// Semi-major axis (AU)
    pub a: f64,
    /// Eccentricity
    pub e: f64,
    /// Inclination (degrees)
    pub i_deg: f64,
    /// Argument of perihelion (degrees)
    pub omega_deg: f64,
    /// Longitude of ascending node (degrees)
    pub omega_big_deg: f64,
    /// Absolute magnitude H
    pub h_mag: f64,
}

impl Etno {
    /// Orbital elements (radians/AU), mean anomaly set to 0.
    pub fn elements(&self) -> OrbitalElements {
        OrbitalElements {
            a: self.a,
            e: self.e,
            i: self.i_deg * DEG2RAD,
            omega: self.omega_deg * DEG2RAD,
            omega_big: self.omega_big_deg * DEG2RAD,
            mean_anomaly: 0.0,
        }
    }

    /// Perihelion distance q = a(1 − e) in AU.
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
    }

    /// Longitude of perihelion ϖ = ω + Ω in radians, wrapped to [0, 2π).
    pub fn longitude_of_perihelion(&self) -> f64 {
        ((self.omega_deg + self.omega_big_deg) * DEG2RAD).rem_euclid(TWO_PI)
    }
}

/// The 10 ETNOs with a ≳ 230 AU from Brown (2017) Table 1.
pub const BROWN_2017_SAMPLE: [Etno; 10] = [
    Etno {
        name: "Sedna",
        a: 506.0,
        e: 0.85,
        i_deg: 11.9,
        omega_deg: 311.5,
        omega_big_deg: 144.5,
        h_mag: 1.6,
    },
    Etno {
        name: "2012 VP113",
        a: 261.0,
        e: 0.69,
        i_deg: 24.1,
        omega_deg: 293.8,
        omega_big_deg: 90.8,
        h_mag: 4.0,
    },
    Etno {
        name: "2013 RF98",
        a: 325.0,
        e: 0.89,
        i_deg: 29.6,
        omega_deg: 316.5,
        omega_big_deg: 67.6,
        h_mag: 8.7,
    },
    Etno {
        name: "2004 VN112",
        a: 327.0,
        e: 0.85,
        i_deg: 25.6,
        omega_deg: 327.1,
        omega_big_deg: 66.0,
        h_mag: 6.4,
    },
    Etno {
        name: "2010 GB174",
        a: 351.0,
        e: 0.86,
        i_deg: 21.5,
        omega_deg: 347.8,
        omega_big_deg: 130.6,
        h_mag: 6.5,
    },
    Etno {
        name: "2000 CR105",
        a: 228.0,
        e: 0.81,
        i_deg: 22.7,
        omega_deg: 317.2,
        omega_big_deg: 128.3,
        h_mag: 6.3,
    },
    Etno {
        name: "2007 TG422",
        a: 501.0,
        e: 0.93,
        i_deg: 18.6,
        omega_deg: 285.7,
        omega_big_deg: 113.0,
        h_mag: 6.2,
    },
    Etno {
        name: "2013 FT28",
        a: 310.0,
        e: 0.86,
        i_deg: 17.3,
        omega_deg: 40.2,
        omega_big_deg: 217.8,
        h_mag: 6.7,
    },
    Etno {
        name: "2014 SR349",
        a: 289.0,
        e: 0.84,
        i_deg: 18.0,
        omega_deg: 341.4,
        omega_big_deg: 34.8,
        h_mag: 6.6,
    },
    Etno {
        name: "2015 RX245",
        a: 430.0,
        e: 0.89,
        i_deg: 12.2,
        omega_deg: 65.2,
        omega_big_deg: 8.6,
        h_mag: 6.2,
    },
];

/// Longitudes of perihelion ϖ (radians) for the whole sample.
pub fn longitudes_of_perihelion() -> Vec<f64> {
    BROWN_2017_SAMPLE
        .iter()
        .map(|k| k.longitude_of_perihelion())
        .collect()
}

/// Semi-major axes (AU) for the whole sample — the single source for
/// implied-resonance analyses.
pub fn semi_major_axes() -> Vec<f64> {
    BROWN_2017_SAMPLE.iter().map(|k| k.a).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::circular::{mean_resultant_length, rayleigh_p_value};

    #[test]
    fn test_sample_size_and_membership() {
        assert_eq!(BROWN_2017_SAMPLE.len(), 10);
        let names: Vec<&str> = BROWN_2017_SAMPLE.iter().map(|k| k.name).collect();
        assert!(names.contains(&"Sedna"));
        assert!(names.contains(&"2012 VP113"));
    }

    #[test]
    fn test_elements_sanity() {
        for kbo in &BROWN_2017_SAMPLE {
            assert!(kbo.a > 220.0, "{}: a = {}", kbo.name, kbo.a);
            assert!(kbo.e > 0.6 && kbo.e < 1.0, "{}: e = {}", kbo.name, kbo.e);
            // All members are detached/scattered objects with q > 30 AU
            // (no Neptune crossers, unlike the corrupted tables this replaces)
            assert!(
                kbo.perihelion() > 30.0,
                "{}: q = {:.1}",
                kbo.name,
                kbo.perihelion()
            );
            assert!(kbo.i_deg >= 0.0 && kbo.i_deg < 40.0);
        }
    }

    #[test]
    fn test_known_values_pinned() {
        let sedna = &BROWN_2017_SAMPLE[0];
        assert_eq!(sedna.a, 506.0);
        assert!((sedna.perihelion() - 75.9).abs() < 0.1);
        // 2013 RF98 at ~325 AU (the corrupted table said 780)
        let rf98 = BROWN_2017_SAMPLE
            .iter()
            .find(|k| k.name == "2013 RF98")
            .unwrap();
        assert_eq!(rf98.a, 325.0);
    }

    #[test]
    fn test_varpi_clustering_present() {
        // The sample's defining feature: clustered longitudes of perihelion.
        let varpis = longitudes_of_perihelion();
        let r_bar = mean_resultant_length(&varpis);
        assert!(r_bar > 0.3, "R̄ = {r_bar:.3}");
        // Rayleigh test (vs isotropy, before bias correction) is significant.
        assert!(rayleigh_p_value(&varpis) < 0.05);
    }
}
