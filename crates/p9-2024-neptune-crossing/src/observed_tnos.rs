//! The 17 well-characterized multi-opposition Neptune-crossing TNOs from the paper.
//!
//! Selection criteria: a > 100 AU, i < 40 deg, q < 30 AU.
//! Orbital elements are approximate barycentric osculating values.

use serde::{Deserialize, Serialize};

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

/// Selection criteria for the Neptune-crossing TNO sample.
pub struct SelectionCriteria {
    pub a_min: f64,
    pub i_max: f64,
    pub q_max: f64,
}

/// Returns the selection criteria used in the paper.
pub fn selection_criteria() -> SelectionCriteria {
    SelectionCriteria {
        a_min: 100.0,
        i_max: 40.0,
        q_max: 30.0,
    }
}

/// Returns the 17 well-characterized multi-opposition TNOs from the paper sample.
///
/// All objects satisfy a > 100 AU, i < 40 deg, q < 30 AU.
/// Orbital elements are approximate values from MPC/JPL.
pub fn observed_sample() -> Vec<NeptuneCrossingTno> {
    vec![
        NeptuneCrossingTno {
            name: "2014 SS349",
            a: 129.2,
            e: 0.793,
            i: 19.8,
            q: 26.7,
        },
        NeptuneCrossingTno {
            name: "2013 UL10",
            a: 109.2,
            e: 0.781,
            i: 7.5,
            q: 23.9,
        },
        NeptuneCrossingTno {
            name: "2015 KH163",
            a: 153.0,
            e: 0.846,
            i: 27.1,
            q: 23.6,
        },
        NeptuneCrossingTno {
            name: "2013 FS28",
            a: 193.6,
            e: 0.854,
            i: 13.1,
            q: 28.3,
        },
        NeptuneCrossingTno {
            name: "2017 FO161",
            a: 225.6,
            e: 0.870,
            i: 10.9,
            q: 29.3,
        },
        NeptuneCrossingTno {
            name: "2010 ER65",
            a: 147.0,
            e: 0.808,
            i: 22.0,
            q: 28.2,
        },
        NeptuneCrossingTno {
            name: "2014 OE394",
            a: 165.1,
            e: 0.838,
            i: 11.7,
            q: 26.7,
        },
        NeptuneCrossingTno {
            name: "2015 RY245",
            a: 130.3,
            e: 0.859,
            i: 4.8,
            q: 18.4,
        },
        NeptuneCrossingTno {
            name: "2013 JO64",
            a: 102.5,
            e: 0.718,
            i: 31.8,
            q: 28.9,
        },
        NeptuneCrossingTno {
            name: "2015 GA58",
            a: 115.5,
            e: 0.771,
            i: 20.2,
            q: 26.5,
        },
        NeptuneCrossingTno {
            name: "2014 QR441",
            a: 139.7,
            e: 0.810,
            i: 16.0,
            q: 26.5,
        },
        NeptuneCrossingTno {
            name: "2016 QV89",
            a: 110.2,
            e: 0.760,
            i: 12.4,
            q: 26.5,
        },
        NeptuneCrossingTno {
            name: "2018 AD39",
            a: 134.0,
            e: 0.810,
            i: 8.6,
            q: 25.5,
        },
        NeptuneCrossingTno {
            name: "2014 LU28",
            a: 101.6,
            e: 0.707,
            i: 24.3,
            q: 29.8,
        },
        NeptuneCrossingTno {
            name: "2012 GA32",
            a: 105.1,
            e: 0.746,
            i: 14.0,
            q: 26.7,
        },
        NeptuneCrossingTno {
            name: "2013 AT183",
            a: 120.0,
            e: 0.784,
            i: 25.5,
            q: 25.9,
        },
        NeptuneCrossingTno {
            name: "2015 DW224",
            a: 140.8,
            e: 0.828,
            i: 6.2,
            q: 24.2,
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sample_size() {
        assert_eq!(observed_sample().len(), 17);
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
        for tno in observed_sample() {
            let q_computed = tno.a * (1.0 - tno.e);
            assert!(
                (tno.q - q_computed).abs() < 1.0,
                "{}: listed q = {}, computed q = {:.1}",
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
