//! The DES six-year extreme-TNO (a > 150 AU, q > 30 AU) subset.
//!
//! Provenance: Bernardinelli et al. (2021), "A search for trans-Neptunian
//! objects with the six-year Dark Energy Survey" (arXiv:2109.03758). The full
//! catalog has 815 TNOs; the paper isolates 16 *extreme* TNOs (a > 150 AU,
//! q > 30 AU) for the clustering / Planet-Nine discussion.
//!
//! We encode the extreme-TNO orbital elements that are also published in the
//! DES isotropy table (Bernardinelli et al. 2020, reproduced in
//! `p9-2020-des-isotropy::des_sample`): seven Y4-discovery extreme TNOs plus
//! the DES recovery 2013 RF98. These are the well-vetted DES eTNOs with
//! published elements; they form a documented subset of the paper's 16.
//! Absolute magnitudes H are carried for the completeness analysis. Element
//! 1σ fit uncertainties are dropped (central values pinned).
//!
//! Residual: this table has 8 objects vs the paper's full 16 extreme TNOs (the
//! remaining 8 are fainter / lower-significance detections without a single
//! widely cross-referenced element set). The clustering conclusion —
//! consistency with isotropy under the DES selection — is robust across the
//! vetted subset; see `super::reference::N_EXTREME_TNOS` for the labelled 16.

use p9_core::constants::{DEG2RAD, TWO_PI};

/// A DES extreme TNO with the elements needed for clustering and completeness.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DesTno {
    pub name: &'static str,
    /// Semi-major axis (AU).
    pub a: f64,
    /// Eccentricity.
    pub e: f64,
    /// Inclination (degrees).
    pub i_deg: f64,
    /// Longitude of ascending node Ω (degrees).
    pub node_deg: f64,
    /// Argument of perihelion ω (degrees).
    pub argp_deg: f64,
    /// Absolute magnitude H (r-band).
    pub h_mag: f64,
}

impl DesTno {
    /// Perihelion distance q = a(1 − e) in AU.
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
    }

    /// Aphelion distance Q = a(1 + e) in AU.
    pub fn aphelion(&self) -> f64 {
        self.a * (1.0 + self.e)
    }

    /// Longitude of ascending node Ω (radians, wrapped to [0, 2π)).
    pub fn node(&self) -> f64 {
        (self.node_deg * DEG2RAD).rem_euclid(TWO_PI)
    }

    /// Argument of perihelion ω (radians, wrapped to [0, 2π)).
    pub fn argp(&self) -> f64 {
        (self.argp_deg * DEG2RAD).rem_euclid(TWO_PI)
    }

    /// Longitude of perihelion ϖ = ω + Ω (radians, wrapped to [0, 2π)).
    pub fn varpi(&self) -> f64 {
        ((self.argp_deg + self.node_deg) * DEG2RAD).rem_euclid(TWO_PI)
    }
}

/// The vetted DES extreme-TNO subset (a > 150 AU, q > 30 AU). Elements from
/// the DES eTNO table; H magnitudes from the DES discovery photometry.
pub const DES_EXTREME_TNOS: [DesTno; 8] = [
    DesTno {
        name: "2013 RA109",
        a: 463.3,
        e: 0.901,
        i_deg: 12.39,
        node_deg: 104.79,
        argp_deg: 262.91,
        h_mag: 6.5,
    },
    DesTno {
        name: "2015 BP519",
        a: 449.38,
        e: 0.922,
        i_deg: 54.11,
        node_deg: 135.21,
        argp_deg: 348.06,
        h_mag: 4.3,
    },
    DesTno {
        name: "2013 SL102",
        a: 314.31,
        e: 0.879,
        i_deg: 6.50,
        node_deg: 94.73,
        argp_deg: 265.49,
        h_mag: 6.4,
    },
    DesTno {
        name: "2014 WB556",
        a: 289.3,
        e: 0.853,
        i_deg: 24.15,
        node_deg: 114.89,
        argp_deg: 234.53,
        h_mag: 6.9,
    },
    DesTno {
        name: "2016 SG58",
        a: 233.0,
        e: 0.849,
        i_deg: 13.22,
        node_deg: 118.97,
        argp_deg: 296.29,
        h_mag: 7.0,
    },
    DesTno {
        name: "2016 QV89",
        a: 171.6,
        e: 0.767,
        i_deg: 21.38,
        node_deg: 173.21,
        argp_deg: 281.08,
        h_mag: 7.6,
    },
    DesTno {
        name: "508338 (2015 SO20)",
        a: 164.7,
        e: 0.799,
        i_deg: 23.41,
        node_deg: 33.63,
        argp_deg: 354.78,
        h_mag: 6.8,
    },
    DesTno {
        name: "2013 RF98",
        a: 358.2,
        e: 0.90,
        i_deg: 23.54,
        node_deg: 67.63,
        argp_deg: 312.05,
        h_mag: 8.7,
    },
];

/// The three orbital angles tested for clustering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Angle {
    /// Longitude of ascending node Ω.
    Node,
    /// Argument of perihelion ω.
    ArgPeri,
    /// Longitude of perihelion ϖ = ω + Ω (the Planet 9 alignment angle).
    Varpi,
}

impl Angle {
    pub fn all() -> [Angle; 3] {
        [Angle::Node, Angle::ArgPeri, Angle::Varpi]
    }

    pub fn label(&self) -> &'static str {
        match self {
            Angle::Node => "Omega",
            Angle::ArgPeri => "omega",
            Angle::Varpi => "varpi",
        }
    }

    /// Extract this angle (radians) from an object.
    pub fn of(&self, o: &DesTno) -> f64 {
        match self {
            Angle::Node => o.node(),
            Angle::ArgPeri => o.argp(),
            Angle::Varpi => o.varpi(),
        }
    }
}

/// Angle values (radians) for a sample.
pub fn angle_values(sample: &[DesTno], angle: Angle) -> Vec<f64> {
    sample.iter().map(|o| angle.of(o)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_objects_are_extreme_tnos() {
        for o in &DES_EXTREME_TNOS {
            assert!(o.a > 150.0, "{}: a = {}", o.name, o.a);
            assert!(
                o.perihelion() > 30.0,
                "{}: q = {:.1}",
                o.name,
                o.perihelion()
            );
        }
    }

    #[test]
    fn vetted_subset_within_paper_count() {
        // 8 vetted objects, a documented subset of the paper's 16 extreme TNOs.
        assert_eq!(DES_EXTREME_TNOS.len(), 8);
        assert!(DES_EXTREME_TNOS.len() <= crate::reference::N_EXTREME_TNOS);
    }

    #[test]
    fn varpi_matches_node_plus_argp() {
        // Spot-check 2013 RA109: 104.79 + 262.91 = 367.70 -> 7.70 deg.
        let o = DES_EXTREME_TNOS
            .iter()
            .find(|o| o.name == "2013 RA109")
            .unwrap();
        let varpi_deg = o.varpi() * p9_core::constants::RAD2DEG;
        assert!((varpi_deg - 7.70).abs() < 0.05, "varpi = {varpi_deg:.2}");
    }
}
