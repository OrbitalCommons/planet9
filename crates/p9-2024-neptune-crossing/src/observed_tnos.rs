//! The 17 well-characterized, multi-opposition Neptune-crossing TNOs of
//! Batygin, Morbidelli, Brown & Nesvorny (2024), ApJL 966, L8
//! (arXiv:2404.11594), Figure 1.
//!
//! Selection criteria (paper Sec. 1 and Fig. 1 caption): semi-major axis
//! a > 100 AU, perihelion q < 30 AU (Neptune-crossing), inclination
//! i < 40 deg, restricted to objects "whose orbits have been quantified
//! through multi-opposition observations" (17 of the 29 MPC objects meeting
//! the orbital cuts at the time of writing).
//!
//! Provenance: orbital elements queried from the JPL Small-Body Database
//! (https://ssd-api.jpl.nasa.gov/sbdb_query.api) with the constraints
//! a > 100 AU, q < 30 AU, i < 40 deg (June 2026). The 17 paper objects were
//! identified by cross-matching each labelled inclination in the paper's
//! Figure 1 (5.7, 6.1, 9.2, 10.8, 11.5, 12.3, 16.5, 19.4, 20.0, 22.4, 26.1,
//! 26.2, 26.3, 26.5, 28.0, 30.8, 32.1 deg) against the SBDB result set; all
//! 17 match uniquely. The previous version of this table contained
//! fabricated entries ("2016 QV89" is a ~30 m NEA, "2014 LU28" is a
//! retrograde centaur); both are gone.

use p9_core::units::{au, degrees, Angle, Length};
use serde::{Deserialize, Serialize};

/// Provenance statement for the table below (tested non-empty and SBDB-sourced).
pub const PROVENANCE: &str = "JPL SBDB query (sbdb_query.api, June 2026): \
a > 100 AU, q < 30 AU, i < 40 deg, asteroid kind; multi-opposition subset \
cross-matched object-by-object against the inclination labels of Figure 1 \
of Batygin, Morbidelli, Brown & Nesvorny (2024), ApJL 966, L8.";

/// A Neptune-crossing TNO with basic orbital parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NeptuneCrossingTno {
    pub name: &'static str,
    /// Semi-major axis in AU
    pub a: f64,
    /// Eccentricity
    pub e: f64,
    /// Inclination in degrees
    pub i: f64,
    /// Perihelion distance in AU (q = a(1-e))
    pub q: f64,
}

impl NeptuneCrossingTno {
    /// Semi-major axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a)
    }
    /// Inclination as a typed [`Angle`].
    pub fn inclination(&self) -> Angle {
        degrees(self.i)
    }
    /// Perihelion distance as a typed [`Length`].
    pub fn perihelion(&self) -> Length {
        au(self.q)
    }
}

/// Selection criteria for the Neptune-crossing TNO sample.
pub struct SelectionCriteria {
    pub a_min: f64,
    pub i_max: f64,
    pub q_max: f64,
}

impl SelectionCriteria {
    /// Minimum semi-major axis as a typed [`Length`].
    pub fn a_min(&self) -> Length {
        au(self.a_min)
    }
    /// Maximum inclination as a typed [`Angle`].
    pub fn i_max(&self) -> Angle {
        degrees(self.i_max)
    }
    /// Maximum perihelion as a typed [`Length`].
    pub fn q_max(&self) -> Length {
        au(self.q_max)
    }
}

/// Returns the selection criteria used in the paper.
pub fn selection_criteria() -> SelectionCriteria {
    SelectionCriteria {
        a_min: 100.0,
        i_max: 40.0,
        q_max: 30.0,
    }
}

/// The 17 multi-opposition TNOs of the paper sample (JPL SBDB osculating
/// heliocentric elements; see module docs for the query and cross-match).
pub fn observed_sample() -> Vec<NeptuneCrossingTno> {
    vec![
        NeptuneCrossingTno {
            name: "54520 (2000 PJ30)",
            a: 121.9,
            e: 0.7656,
            i: 5.72,
            q: 28.560,
        },
        NeptuneCrossingTno {
            name: "87269 (2000 OO67)",
            a: 597.1,
            e: 0.9651,
            i: 20.05,
            q: 20.846,
        },
        NeptuneCrossingTno {
            name: "308933 (2006 SQ372)",
            a: 805.4,
            e: 0.9699,
            i: 19.47,
            q: 24.214,
        },
        NeptuneCrossingTno {
            name: "353222 (2009 YD7)",
            a: 126.1,
            e: 0.8938,
            i: 30.76,
            q: 13.390,
        },
        NeptuneCrossingTno {
            name: "469750 (2005 PU21)",
            a: 179.2,
            e: 0.8369,
            i: 6.17,
            q: 29.234,
        },
        NeptuneCrossingTno {
            name: "523771 (2014 XP40)",
            a: 103.6,
            e: 0.7313,
            i: 11.53,
            q: 27.840,
        },
        NeptuneCrossingTno {
            name: "600217 (2011 QY100)",
            a: 177.8,
            e: 0.8899,
            i: 26.22,
            q: 19.581,
        },
        NeptuneCrossingTno {
            name: "767044 (2014 QW510)",
            a: 217.2,
            e: 0.8804,
            i: 26.12,
            q: 25.977,
        },
        NeptuneCrossingTno {
            name: "2012 GU11",
            a: 175.8,
            e: 0.8967,
            i: 10.76,
            q: 18.162,
        },
        NeptuneCrossingTno {
            name: "2013 AZ60",
            a: 428.1,
            e: 0.9816,
            i: 16.54,
            q: 7.883,
        },
        NeptuneCrossingTno {
            name: "2013 GJ138",
            a: 123.0,
            e: 0.7877,
            i: 27.96,
            q: 26.107,
        },
        NeptuneCrossingTno {
            name: "2013 GW141",
            a: 566.4,
            e: 0.9585,
            i: 32.08,
            q: 23.527,
        },
        NeptuneCrossingTno {
            name: "2014 NV65",
            a: 113.0,
            e: 0.8033,
            i: 12.33,
            q: 22.232,
        },
        NeptuneCrossingTno {
            name: "2014 RK86",
            a: 372.7,
            e: 0.9410,
            i: 26.39,
            q: 22.005,
        },
        NeptuneCrossingTno {
            name: "2014 UY224",
            a: 133.5,
            e: 0.8473,
            i: 26.44,
            q: 20.390,
        },
        NeptuneCrossingTno {
            name: "2015 VD168",
            a: 116.5,
            e: 0.7786,
            i: 22.45,
            q: 25.789,
        },
        NeptuneCrossingTno {
            name: "2021 CP5",
            a: 414.8,
            e: 0.9745,
            i: 9.23,
            q: 10.564,
        },
    ]
}

/// SBDB search string (`sstr`) for a table entry: the bare number for
/// numbered objects ("54520 (2000 PJ30)" → "54520"), the designation itself
/// otherwise.
pub fn sbdb_designation(name: &'static str) -> &'static str {
    name.split(" (").next().unwrap_or(name)
}

/// The sample as snapshot rows for `p9_core::data::refresh` (only a, e, i
/// are vetted columns of this table; q is derived and the angles are not
/// carried).
pub fn element_snapshots() -> Vec<p9_core::data::refresh::ElementSnapshot> {
    observed_sample()
        .iter()
        .map(|t| p9_core::data::refresh::ElementSnapshot {
            name: t.name,
            designation: sbdb_designation(t.name),
            a: Some(t.a),
            e: Some(t.e),
            i_deg: Some(t.i),
            omega_deg: None,
            omega_big_deg: None,
            h_mag: None,
        })
        .collect()
}

/// Diff the frozen sample against the live JPL SBDB (network; see the
/// `#[ignore]`d test in `tests/sbdb_refresh_live.rs`). Empty result = the
/// table still matches JPL within `Tolerances::full_precision()`.
#[cfg(feature = "sbdb-refresh")]
pub fn refresh_from_sbdb(
    client: &p9_core::data::refresh::SbdbClient,
) -> Result<Vec<p9_core::data::refresh::EtnoDiff>, String> {
    use p9_core::data::refresh::refresh_table_from_sbdb;
    refresh_table_from_sbdb(
        client,
        &element_snapshots(),
        &p9_core::data::refresh::Tolerances::full_precision(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sample_size() {
        assert_eq!(observed_sample().len(), 17);
    }

    #[test]
    fn typed_accessors_match_f64() {
        use approx::assert_relative_eq;
        use uom::si::angle::degree;
        use uom::si::length::astronomical_unit;
        let t = &observed_sample()[0];
        assert_relative_eq!(
            t.semi_major_axis().get::<astronomical_unit>(),
            t.a,
            epsilon = 1e-12
        );
        assert_relative_eq!(t.inclination().get::<degree>(), t.i, epsilon = 1e-12);
        assert_relative_eq!(
            t.perihelion().get::<astronomical_unit>(),
            t.q,
            epsilon = 1e-12
        );
        let c = selection_criteria();
        assert_relative_eq!(
            c.a_min().get::<astronomical_unit>(),
            c.a_min,
            epsilon = 1e-12
        );
        assert_relative_eq!(c.i_max().get::<degree>(), c.i_max, epsilon = 1e-12);
        assert_relative_eq!(
            c.q_max().get::<astronomical_unit>(),
            c.q_max,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_sbdb_designations() {
        assert_eq!(sbdb_designation("54520 (2000 PJ30)"), "54520");
        assert_eq!(sbdb_designation("2012 GU11"), "2012 GU11");
        let snaps = element_snapshots();
        assert_eq!(snaps.len(), 17);
        assert!(snaps
            .iter()
            .all(|s| !s.designation.is_empty() && !s.designation.contains('(')));
    }

    #[test]
    fn test_provenance_documented() {
        assert!(PROVENANCE.contains("SBDB"), "provenance must cite JPL SBDB");
        assert!(
            PROVENANCE.contains("2024"),
            "provenance must cite the paper"
        );
    }

    #[test]
    fn test_fabricated_objects_removed() {
        // "2016 QV89" is a ~30 m near-Earth asteroid and "2014 LU28" a
        // retrograde centaur; neither belongs in this sample.
        for tno in observed_sample() {
            assert!(!tno.name.contains("2016 QV89"), "fabricated entry");
            assert!(!tno.name.contains("2014 LU28"), "fabricated entry");
        }
    }

    #[test]
    fn test_all_above_100_au() {
        for tno in observed_sample() {
            assert!(
                tno.a > 100.0,
                "{} has a = {} AU, expected > 100",
                tno.name,
                tno.a,
            );
        }
    }

    #[test]
    fn test_all_low_inclination() {
        let criteria = selection_criteria();
        for tno in observed_sample() {
            assert!(
                tno.i < criteria.i_max,
                "{} has i = {} deg, expected < {}",
                tno.name,
                tno.i,
                criteria.i_max,
            );
        }
    }

    #[test]
    fn test_all_neptune_crossing() {
        let criteria = selection_criteria();
        for tno in observed_sample() {
            assert!(
                tno.q < criteria.q_max,
                "{} has q = {} AU, expected < {}",
                tno.name,
                tno.q,
                criteria.q_max,
            );
        }
    }

    #[test]
    fn test_perihelion_consistency() {
        // q must equal a(1-e) to within the rounding of the published
        // elements (a to 0.1 AU, e to 1e-4 => |dq| <~ 0.07 AU).
        for tno in observed_sample() {
            let q_computed = tno.a * (1.0 - tno.e);
            assert!(
                (tno.q - q_computed).abs() < 0.1,
                "{}: listed q = {}, computed q = {:.3}",
                tno.name,
                tno.q,
                q_computed,
            );
        }
    }

    #[test]
    fn test_unique_names() {
        let sample = observed_sample();
        let mut names: Vec<&str> = sample.iter().map(|t| t.name).collect();
        names.sort();
        names.dedup();
        assert_eq!(names.len(), 17);
    }
}
