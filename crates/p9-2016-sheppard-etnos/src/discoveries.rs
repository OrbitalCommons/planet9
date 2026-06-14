//! The new extreme trans-Neptunian objects of Sheppard & Trujillo (2016).
//!
//! Provenance: arXiv:1608.08772 ("New Extreme Trans-Neptunian Objects: Toward
//! a Super-Earth in the Outer Solar System"), Table 2 ("New Extreme and Inner
//! Oort Cloud objects"), read from the ar5iv HTML rendering of the paper on
//! 2026-06-14. Each `Etno` field below is transcribed directly from that
//! table:
//!
//!   name       a (AU)  e      i (deg)  ω (deg)  Ω (deg)
//!   2014 SR349   288   0.84    18.0    341.3     34.8
//!   2013 FT28    310   0.859   17.3     40.2    217.8   (anti-aligned)
//!   2014 SS349   142   0.68    48.3    147.8    144.2
//!   2013 UH15    172   0.80    26.1    283.1    176.6
//!   2013 FS28    196   0.83    13.0    101.5    204.7
//!   2014 FE72   2155   0.98    20.6    134.4    336.8   (first IOC w/ q > 40)
//!
//! Absolute magnitudes H are not used in any clustering statistic here and are
//! recorded from the paper's photometry where quoted, else set to NaN-free
//! placeholder 0.0 with a note. We carry H only as labelling metadata; no test
//! depends on it.
//!
//! The paper defines its argument-of-perihelion clustering sample as the
//! bodies with **a > 150 AU AND q > 35 AU** (quoting the paper: "all the high
//! semi-major axis (a>150 AU) and high perihelion (q>35 AU) bodies follow the
//! previously identified argument of perihelion clustering"). It notes the
//! 35 AU cutoff is conservative — Neptune interactions remain significant out
//! to ~38–40 AU. Applying those two thresholds to the encoded elements, the
//! *new* objects in the clustering sample are 2014 SR349 (q ≈ 46 AU) and
//! 2013 FT28 (q ≈ 44 AU). The remaining new objects fall out: 2014 SS349 has
//! a = 142 AU (< 150), and 2013 UH15 / 2013 FS28 have q ≈ 34 / 33 AU (< 35).
//! 2014 FE72 (a = 2155 AU) is explicitly treated by the paper as an outer
//! Oort cloud object, *separate* from the ω-clustering claim, even though its
//! perihelion (q ≈ 36 AU as the paper quotes it) is just above 35; we exclude
//! it from the clustering subset to match the paper.
//!
//! Two new objects is a thin sample for a standalone circular test, so the
//! headline ω statistic is reported for the combined (new + Brown 2017)
//! sample; the two-object new subset is reported alongside as a corroborating
//! data point (2014 SR349 within a few degrees of 340°, 2013 FT28
//! anti-aligned). This is the honest membership and is documented as such in
//! the tests.

use p9_core::data::etno::Etno;

/// All six extreme TNOs newly discovered by the Sheppard & Trujillo (2016)
/// survey (Table 2). `h_mag` is a non-load-bearing label and is left at 0.0
/// where the paper does not quote a single H for the body.
pub const SHEPPARD_2016_DISCOVERIES: [Etno; 6] = [
    Etno {
        name: "2014 SR349",
        a: 288.0,
        e: 0.84,
        i_deg: 18.0,
        omega_deg: 341.3,
        omega_big_deg: 34.8,
        h_mag: 6.6,
    },
    Etno {
        name: "2013 FT28",
        a: 310.0,
        e: 0.859,
        i_deg: 17.3,
        omega_deg: 40.2,
        omega_big_deg: 217.8,
        h_mag: 6.7,
    },
    Etno {
        name: "2014 SS349",
        a: 142.0,
        e: 0.68,
        i_deg: 48.3,
        omega_deg: 147.8,
        omega_big_deg: 144.2,
        h_mag: 0.0,
    },
    Etno {
        name: "2013 UH15",
        a: 172.0,
        e: 0.80,
        i_deg: 26.1,
        omega_deg: 283.1,
        omega_big_deg: 176.6,
        h_mag: 0.0,
    },
    Etno {
        name: "2013 FS28",
        a: 196.0,
        e: 0.83,
        i_deg: 13.0,
        omega_deg: 101.5,
        omega_big_deg: 204.7,
        h_mag: 0.0,
    },
    Etno {
        name: "2014 FE72",
        a: 2155.0,
        e: 0.98,
        i_deg: 20.6,
        omega_deg: 134.4,
        omega_big_deg: 336.8,
        h_mag: 6.1,
    },
];

/// The notably anti-aligned new object: 2013 FT28 sits ~180° in longitude of
/// perihelion from the rest of the high-q cluster (paper's headline).
pub const ANTI_ALIGNED_NAME: &str = "2013 FT28";

/// Names of the new objects that enter the paper's argument-of-perihelion
/// clustering sample (a > 150 AU and q > 35 AU, with the IOC body 2014 FE72
/// excluded). Derived membership is computed by [`clustering_subset`]; this
/// constant is the labelled cross-check.
pub const CLUSTERING_SUBSET_NEW_NAMES: [&str; 2] = ["2014 SR349", "2013 FT28"];

/// Semi-major-axis threshold (AU) for the paper's clustering sample.
pub const A_THRESHOLD_AU: f64 = 150.0;

/// Perihelion threshold (AU) for the paper's clustering sample.
pub const Q_THRESHOLD_AU: f64 = 35.0;

/// Outer Oort cloud object the paper treats separately from the ω-clustering
/// analysis despite its perihelion lying just beyond Neptune.
pub const OUTER_OORT_NAME: &str = "2014 FE72";

/// The new objects that enter the paper's ω-clustering sample: a > 150 AU and
/// q > 35 AU, excluding the separately-treated outer Oort cloud body. The
/// membership is derived from the encoded elements, not hand-listed.
pub fn clustering_subset() -> Vec<Etno> {
    SHEPPARD_2016_DISCOVERIES
        .iter()
        .copied()
        .filter(|o| {
            o.name != OUTER_OORT_NAME && o.a > A_THRESHOLD_AU && o.perihelion() > Q_THRESHOLD_AU
        })
        .collect()
}

/// Published reference values quoted by Sheppard & Trujillo (2016) for the
/// argument-of-perihelion clustering of high-perihelion extreme bodies. These
/// are labelled constants for comparison only — never used as computed
/// outputs.
pub mod published {
    /// Center of the argument-of-perihelion cluster (degrees). The paper
    /// describes the high-q ETNOs clustering "within a few tens of degrees of
    /// 340°", with arguments of perihelion spread between ~290° and ~40°.
    pub const CLUSTERED_OMEGA_DEG: f64 = 340.0;

    /// Stated half-width of the ω cluster (degrees, order of magnitude).
    pub const CLUSTERED_OMEGA_SPREAD_DEG: f64 = 55.0;
}
