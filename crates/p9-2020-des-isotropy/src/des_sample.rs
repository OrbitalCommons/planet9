//! DES Y4 extreme trans-Neptunian object (eTNO) sample.
//!
//! Provenance: Bernardinelli et al. (2020), "Testing the isotropy of the Dark
//! Energy Survey's extreme trans-Neptunian objects" (arXiv:2003.08901), the
//! object table (the DES Y4 eTNOs plus 2013 RF98, recovered by DES). Elements
//! transcribed verbatim from the paper's table; parenthetical 1σ fit
//! uncertainties on the last digit(s) are dropped (we pin the central values).
//!
//! The paper defines four nested samples ("Cases") by semi-major-axis and
//! perihelion cuts:
//!   * Case 1: a > 150 AU, q > 30 AU  → 7 objects
//!   * Case 2: a > 250 AU, q > 30 AU  → 4 objects
//!   * Case 3: a > 150 AU, q > 37 AU  → 4 objects
//!   * Case 4: a > 250 AU, q > 37 AU  → 3 objects
//! and applies Kuiper's test to Ω, ω and ϖ in each Case (12 tests total). The
//! headline: the DES data alone are consistent with azimuthal isotropy.
//!
//! 2013 RF98 (the "⋆" object) is included in the table because DES recovered
//! it, but it sits outside the Y4-discovery Cases; it is carried here as a
//! separate entry and excluded from the Case membership predicates.

use p9_core::constants::{DEG2RAD, TWO_PI};

/// A DES eTNO with the orbital angles relevant to the isotropy tests.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DesEtno {
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
    /// Whether the object is a Y4-discovery (true) or only a DES recovery
    /// (false, the "⋆"-only 2013 RF98).
    pub y4_discovery: bool,
}

impl DesEtno {
    /// Perihelion distance q = a(1 − e) in AU.
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
    }

    /// Longitude of ascending node Ω in radians, wrapped to [0, 2π).
    pub fn node(&self) -> f64 {
        (self.node_deg * DEG2RAD).rem_euclid(TWO_PI)
    }

    /// Argument of perihelion ω in radians, wrapped to [0, 2π).
    pub fn argp(&self) -> f64 {
        (self.argp_deg * DEG2RAD).rem_euclid(TWO_PI)
    }

    /// Longitude of perihelion ϖ = ω + Ω in radians, wrapped to [0, 2π).
    pub fn varpi(&self) -> f64 {
        ((self.argp_deg + self.node_deg) * DEG2RAD).rem_euclid(TWO_PI)
    }
}

/// The DES eTNO table from Bernardinelli et al. (2020).
pub const DES_ETNOS: [DesEtno; 8] = [
    DesEtno {
        name: "2013 RA109",
        a: 463.3,
        e: 0.901,
        i_deg: 12.39,
        node_deg: 104.79,
        argp_deg: 262.91,
        y4_discovery: true,
    },
    DesEtno {
        name: "2015 BP519",
        a: 449.38,
        e: 0.922,
        i_deg: 54.11,
        node_deg: 135.21,
        argp_deg: 348.06,
        y4_discovery: true,
    },
    DesEtno {
        name: "2013 SL102",
        a: 314.31,
        e: 0.879,
        i_deg: 6.50,
        node_deg: 94.73,
        argp_deg: 265.49,
        y4_discovery: true,
    },
    DesEtno {
        name: "2014 WB556",
        a: 289.3,
        e: 0.853,
        i_deg: 24.15,
        node_deg: 114.89,
        argp_deg: 234.53,
        y4_discovery: true,
    },
    DesEtno {
        name: "2016 SG58",
        a: 233.0,
        e: 0.849,
        i_deg: 13.22,
        node_deg: 118.97,
        argp_deg: 296.29,
        y4_discovery: true,
    },
    DesEtno {
        name: "2016 QV89",
        a: 171.6,
        e: 0.767,
        i_deg: 21.38,
        node_deg: 173.21,
        argp_deg: 281.08,
        y4_discovery: true,
    },
    DesEtno {
        name: "508338 (2015 SO20)",
        a: 164.7,
        e: 0.799,
        i_deg: 23.41,
        node_deg: 33.63,
        argp_deg: 354.78,
        y4_discovery: true,
    },
    DesEtno {
        name: "2013 RF98",
        a: 358.2,
        e: 0.90,
        i_deg: 23.54,
        node_deg: 67.63,
        argp_deg: 312.05,
        y4_discovery: false,
    },
];

/// The four sample definitions of the paper.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SampleCase {
    /// a > 150 AU, q > 30 AU (7 objects).
    Case1,
    /// a > 250 AU, q > 30 AU (4 objects).
    Case2,
    /// a > 150 AU, q > 37 AU (4 objects).
    Case3,
    /// a > 250 AU, q > 37 AU (3 objects).
    Case4,
}

impl SampleCase {
    /// (a_min, q_min) cuts in AU for this Case.
    pub fn cuts(&self) -> (f64, f64) {
        match self {
            SampleCase::Case1 => (150.0, 30.0),
            SampleCase::Case2 => (250.0, 30.0),
            SampleCase::Case3 => (150.0, 37.0),
            SampleCase::Case4 => (250.0, 37.0),
        }
    }

    pub fn all() -> [SampleCase; 4] {
        [
            SampleCase::Case1,
            SampleCase::Case2,
            SampleCase::Case3,
            SampleCase::Case4,
        ]
    }
}

/// The Y4-discovery eTNOs passing the cuts of a given Case (the "⋆"-only
/// 2013 RF98 is never included, matching the paper's Case definitions).
pub fn case_sample(case: SampleCase) -> Vec<DesEtno> {
    let (a_min, q_min) = case.cuts();
    DES_ETNOS
        .iter()
        .filter(|o| o.y4_discovery && o.a > a_min && o.perihelion() > q_min)
        .copied()
        .collect()
}

/// The three orbital angles tested for isotropy in the paper.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Angle {
    /// Longitude of ascending node Ω.
    Node,
    /// Argument of perihelion ω.
    ArgPeri,
    /// Longitude of perihelion ϖ = ω + Ω.
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
    pub fn of(&self, o: &DesEtno) -> f64 {
        match self {
            Angle::Node => o.node(),
            Angle::ArgPeri => o.argp(),
            Angle::Varpi => o.varpi(),
        }
    }
}

/// Pull out a vector of one angle (radians) for a sample of objects.
pub fn angle_values(sample: &[DesEtno], angle: Angle) -> Vec<f64> {
    sample.iter().map(|o| angle.of(o)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_case_counts_match_paper() {
        assert_eq!(case_sample(SampleCase::Case1).len(), 7);
        assert_eq!(case_sample(SampleCase::Case2).len(), 4);
        assert_eq!(case_sample(SampleCase::Case3).len(), 4);
        assert_eq!(case_sample(SampleCase::Case4).len(), 3);
    }

    #[test]
    fn test_varpi_consistent_with_paper() {
        // Paper lists ϖ = Ω + ω; spot-check 2013 RA109 (104.79 + 262.91 = 367.70 → 7.70°).
        let ra109 = DES_ETNOS.iter().find(|o| o.name == "2013 RA109").unwrap();
        let varpi_deg = ra109.varpi() * p9_core::constants::RAD2DEG;
        assert!((varpi_deg - 7.70).abs() < 0.05, "varpi = {varpi_deg:.2}");
    }

    #[test]
    fn test_all_objects_are_etnos() {
        for o in &DES_ETNOS {
            assert!(o.a > 150.0, "{}: a = {}", o.name, o.a);
            assert!(
                o.perihelion() > 30.0,
                "{}: q = {:.1}",
                o.name,
                o.perihelion()
            );
        }
    }
}
