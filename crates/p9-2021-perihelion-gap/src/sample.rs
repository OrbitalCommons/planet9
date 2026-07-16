//! The observed distant-TNO perihelion sample.
//!
//! Oldroyd & Trujillo (2021) study the perihelion distribution of the distant
//! minor planets that straddle the scattering/detached divide: ETNOs with
//! q ≳ 40 AU, a ≳ 150 AU and the high-perihelion inner-Oort-cloud (IOC)
//! Sednoids with q ≳ 65 AU. To see *both* flanks of the gap we need the
//! high-q population (Sedna, 2012 VP113) and a broader low-q ETNO set than the
//! 10-object a ≳ 230 AU clustering sample alone.
//!
//! We therefore use:
//!   - the vetted `p9_core::data::etno::BROWN_2017_SAMPLE` (single workspace
//!     source for the a ≳ 230 AU objects), and
//!   - a documented `EXTENDED_DISTANT_TNOS` table for the additional distant,
//!     high-e objects (lower-a ETNOs and detached Sednoids) needed to populate
//!     the flanks. Objects already present in the Brown sample are NOT
//!     duplicated here.
//!
//! Provenance of the extended table: osculating elements re-transcribed from
//! JPL SBDB 2026-07 (a and e to map perihelion; angles are not used here so
//! are omitted). The previous table carried several rows whose elements
//! matched no JPL solution — most damagingly 2015 KH163 coded into the gap at
//! q ≈ 58 (real q ≈ 40) and the genuine gap object 2021 RR205 coded onto the
//! high flank at q ≈ 89 (real q ≈ 55.6) — which manufactured a deeper deficit
//! than the data support. 2014 SS349 (a ≈ 146) and 2010 ER65 (a ≈ 98) fail
//! the a ≳ 150 AU membership criterion outright and were removed.

use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::units::{au, Length};

/// A distant high-e TNO entered by (a, e). Only the perihelion q = a(1−e) and
/// the eccentricity are used by the gap analysis, so we keep the record
/// minimal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistantTno {
    pub name: &'static str,
    /// Semi-major axis (AU).
    pub a: f64,
    /// Eccentricity.
    pub e: f64,
}

impl DistantTno {
    /// Semi-major axis as a dimension-checked [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a)
    }

    /// Perihelion distance q = a(1 − e) as a dimension-checked [`Length`].
    pub fn perihelion_typed(&self) -> Length {
        au(self.a * (1.0 - self.e))
    }
}

/// Extended distant high-e TNO table (objects NOT already in
/// `BROWN_2017_SAMPLE`). Spans the low-q scattering ETNOs, the gap region, and
/// the high-q detached Sednoids so both flanks of the perihelion gap are
/// represented.
///
/// Membership criterion (matching the paper's distant sample): a ≳ 150 AU and
/// q ≳ 40 AU, plus the canonical Sednoids. The deliberate sparseness in
/// q ≈ 50–65 AU is the observed gap, not a curation choice — these are the
/// known distant high-e objects in this q range.
pub const EXTENDED_DISTANT_TNOS: [DistantTno; 10] = [
    // --- Low-q scattering-coupled ETNOs (q ~ 35-50 AU) ---
    DistantTno {
        name: "2014 FE72",
        a: 2690.0,
        e: 0.987,
    }, // q ~ 36 (scattering)
    DistantTno {
        name: "2015 KG163",
        a: 640.0,
        e: 0.937,
    }, // q ~ 40.5
    DistantTno {
        name: "2013 SY99",
        a: 813.0,
        e: 0.939,
    }, // q ~ 49.9 (top of the low flank)
    DistantTno {
        name: "2014 WB556",
        a: 285.0,
        e: 0.850,
    }, // q ~ 42.8
    DistantTno {
        name: "2016 SD106",
        a: 357.0,
        e: 0.881,
    }, // q ~ 42.7
    DistantTno {
        name: "2015 GT50",
        a: 314.0,
        e: 0.878,
    }, // q ~ 38.3 (just below floor, scattering edge)
    DistantTno {
        name: "2002 GB32",
        a: 205.0,
        e: 0.828,
    }, // q ~ 35.3 (scattering)
    DistantTno {
        name: "2015 KH163",
        a: 151.0,
        e: 0.735,
    }, // q ~ 39.9 (low flank — NOT a gap object; the old table miscoded it at q ~ 58)
    // --- Gap region (q ~ 50-65 AU): sparse, the observed deficit.
    DistantTno {
        name: "2021 RR205",
        a: 949.0,
        e: 0.941,
    }, // q ~ 55.6 (a genuine gap-region object)
    // --- High-q detached / inner-Oort-cloud Sednoids (q ~ 65+ AU) ---
    DistantTno {
        name: "2015 TG387",
        a: 1350.0,
        e: 0.952,
    }, // q ~ 64.7 (Leleakuhonua — at the gap's upper edge)
];

/// The full perihelion sample (AU): perihelia of the vetted Brown sample plus
/// the extended distant-TNO table. This is the observed q distribution the
/// gap analysis operates on.
pub fn observed_perihelia() -> Vec<f64> {
    BROWN_2017_SAMPLE
        .iter()
        .map(|k| k.perihelion())
        .chain(
            EXTENDED_DISTANT_TNOS
                .iter()
                .map(|t| (t.perihelion_typed() / au(1.0)).value),
        )
        .collect()
}

/// The perihelion sample as it stood at the paper's 2021 epoch: excludes
/// 2021 RR205 (announced 2022, after Oldroyd & Trujillo's submission), whose
/// q ≈ 55.6 AU sits squarely inside the published gap window. Use this for
/// reproducing the paper's statistic; use [`observed_perihelia`] for the
/// current state of knowledge (where the gap has partially filled in).
pub fn paper_epoch_perihelia() -> Vec<f64> {
    BROWN_2017_SAMPLE
        .iter()
        .map(|k| k.perihelion())
        .chain(
            EXTENDED_DISTANT_TNOS
                .iter()
                .filter(|t| t.name != "2021 RR205")
                .map(|t| t.a * (1.0 - t.e)),
        )
        .collect()
}

/// (perihelion, eccentricity) pairs for the full sample, for selecting the
/// high-e (e ≳ 0.65) subset in which the gap is cleanest.
pub fn observed_perihelia_with_e() -> Vec<(f64, f64)> {
    BROWN_2017_SAMPLE
        .iter()
        .map(|k| (k.perihelion(), k.e))
        .chain(
            EXTENDED_DISTANT_TNOS
                .iter()
                .map(|t| ((t.perihelion_typed() / au(1.0)).value, t.e)),
        )
        .collect()
}

/// Perihelia (AU) of the high-eccentricity subset (e ≥ `e_floor`) — the
/// population in which Oldroyd & Trujillo report the clearest gap.
pub fn high_e_perihelia(e_floor: f64) -> Vec<f64> {
    observed_perihelia_with_e()
        .into_iter()
        .filter(|(_, e)| *e >= e_floor)
        .map(|(q, _)| q)
        .collect()
}

/// The extended table as snapshot rows for `p9_core::data::refresh` (only a
/// and e are carried by this table; q is derived).
pub fn element_snapshots() -> Vec<p9_core::data::refresh::ElementSnapshot> {
    EXTENDED_DISTANT_TNOS
        .iter()
        .map(|t| p9_core::data::refresh::ElementSnapshot {
            name: t.name,
            designation: t.name,
            a: Some(t.a),
            e: Some(t.e),
            i_deg: None,
            omega_deg: None,
            omega_big_deg: None,
            h_mag: None,
        })
        .collect()
}

/// Diff the frozen extended table against the live JPL SBDB (network; see
/// the `#[ignore]`d test in `tests/sbdb_refresh_live.rs`). SBDB serves ~3
/// significant figures for the weakly-constrained members, and a/e drift
/// with the fit degeneracy, so the coarse `brown2017` allowances apply.
#[cfg(feature = "sbdb-refresh")]
pub fn refresh_from_sbdb(
    client: &p9_core::data::refresh::SbdbClient,
) -> Result<Vec<p9_core::data::refresh::EtnoDiff>, String> {
    p9_core::data::refresh::refresh_table_from_sbdb(
        client,
        &element_snapshots(),
        &p9_core::data::refresh::Tolerances::brown2017(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::published::{ETNO_Q_FLOOR_AU, GAP_ECCENTRICITY_FLOOR, IOC_Q_FLOOR_AU};
    use approx::assert_relative_eq;
    use uom::si::length::astronomical_unit;

    #[test]
    fn typed_accessors_match_f64() {
        let t = DistantTno {
            name: "test",
            a: 500.0,
            e: 0.9,
        };
        assert_relative_eq!(
            t.semi_major_axis().get::<astronomical_unit>(),
            t.a,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            t.perihelion_typed().get::<astronomical_unit>(),
            t.a * (1.0 - t.e),
            epsilon = 1e-9
        );
    }

    #[test]
    fn extended_table_has_no_brown_duplicates() {
        for t in &EXTENDED_DISTANT_TNOS {
            assert!(
                !BROWN_2017_SAMPLE.iter().any(|k| k.name == t.name),
                "{} duplicates a Brown-sample object",
                t.name
            );
        }
    }

    #[test]
    fn perihelia_are_distant_and_high_e() {
        for t in &EXTENDED_DISTANT_TNOS {
            assert!(t.a >= 150.0, "{}: a = {}", t.name, t.a);
            assert!(t.e > 0.5 && t.e < 1.0, "{}: e = {}", t.name, t.e);
            // All in the distant trans-Neptunian q range.
            let q = t.perihelion_typed().get::<astronomical_unit>();
            assert!(q > 30.0 && q < 120.0, "{}: q = {:.1}", t.name, q);
        }
    }

    #[test]
    fn sample_spans_both_flanks_of_the_gap() {
        let qs = observed_perihelia();
        // The combined sample has both a low-q scattering flank (q < 50)
        // and a high-q detached flank (q > 65). With SBDB-verified elements
        // the strict high flank is thin — Sedna and 2012 VP113 (2015 TG387's
        // current solution sits at q = 64.8, the gap window's upper edge).
        let low = qs.iter().filter(|&&q| q < 50.0).count();
        let high = qs.iter().filter(|&&q| q >= IOC_Q_FLOOR_AU).count();
        assert!(low >= 5, "scattering flank only {low} objects");
        assert!(high >= 2, "detached flank only {high} objects");
    }

    #[test]
    fn paper_epoch_sample_excludes_only_rr205() {
        let all = observed_perihelia();
        let paper = paper_epoch_perihelia();
        assert_eq!(all.len(), paper.len() + 1);
        // The excluded object is the in-gap 2021 RR205 (q ~ 55.6).
        assert!(all.iter().any(|&q| (q - 55.6).abs() < 0.5));
        assert!(!paper.iter().any(|&q| (q - 55.6).abs() < 0.5));
    }

    #[test]
    fn high_e_subset_is_a_subset() {
        let all = observed_perihelia();
        let hi = high_e_perihelia(GAP_ECCENTRICITY_FLOOR);
        assert!(hi.len() <= all.len());
        assert!(!hi.is_empty());
        // Sedna (e = 0.85, q ~ 76) must survive the e cut.
        assert!(hi.iter().any(|&q| (q - 75.9).abs() < 1.0));
    }

    #[test]
    fn etno_q_floor_documented() {
        // Most of the sample sits above the paper's q ~ 40 AU ETNO floor.
        let qs = observed_perihelia();
        let above = qs.iter().filter(|&&q| q >= ETNO_Q_FLOOR_AU).count();
        assert!(above >= qs.len() / 2);
    }
}
