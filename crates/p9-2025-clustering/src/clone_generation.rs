//! TNO clone generation from observational uncertainties.
//!
//! 25 clones per TNO, orbital elements sampled from normal distributions
//! around measured values. Dataset: a > 150 AU (or a + sigma_a > 150 AU),
//! excluding sigma_a > 100 AU.

use rand::Rng;
use rand_distr::{Distribution, Normal};
use serde::{Deserialize, Serialize};

/// Observed TNO with uncertainties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObservedTno {
    pub name: &'static str,
    /// Semi-major axis (AU)
    pub a: f64,
    /// Uncertainty in a (AU)
    pub sigma_a: f64,
    /// Eccentricity
    pub e: f64,
    /// Uncertainty in e
    pub sigma_e: f64,
    /// Inclination (rad)
    pub i: f64,
    /// Uncertainty in i (rad)
    pub sigma_i: f64,
    /// Longitude of perihelion (rad)
    pub varpi: f64,
    /// Longitude of ascending node (rad)
    pub omega_big: f64,
}

/// A clone of a TNO with sampled orbital elements.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TnoClone {
    pub parent_name: &'static str,
    pub clone_index: usize,
    pub a: f64,
    pub e: f64,
    pub i: f64,
    pub varpi: f64,
    pub omega_big: f64,
}

impl TnoClone {
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
    }
}

/// Generate clones for a single TNO.
pub fn generate_clones<R: Rng>(tno: &ObservedTno, n_clones: usize, rng: &mut R) -> Vec<TnoClone> {
    let a_dist = Normal::new(tno.a, tno.sigma_a.max(0.01)).unwrap();
    let e_dist = Normal::new(tno.e, tno.sigma_e.max(0.001)).unwrap();
    let i_dist = Normal::new(tno.i, tno.sigma_i.max(0.001)).unwrap();

    let mut clones = Vec::with_capacity(n_clones);

    for idx in 0..n_clones {
        let a = a_dist.sample(rng).max(50.0);
        let e = e_dist.sample(rng).clamp(0.01, 0.999);
        let i = i_dist.sample(rng).max(0.0);

        clones.push(TnoClone {
            parent_name: tno.name,
            clone_index: idx,
            a,
            e,
            i,
            varpi: tno.varpi,
            omega_big: tno.omega_big,
        });
    }

    clones
}

/// Selection criteria from the paper.
pub fn passes_selection(tno: &ObservedTno) -> bool {
    (tno.a > 150.0 || tno.a + tno.sigma_a > 150.0) && tno.sigma_a <= 100.0
}

/// Sample of distant TNOs used in the analysis.
pub fn distant_tno_sample() -> Vec<ObservedTno> {
    use p9_core::constants::DEG2RAD;
    vec![
        ObservedTno {
            name: "Sedna",
            a: 506.0,
            sigma_a: 6.0,
            e: 0.843,
            sigma_e: 0.001,
            i: 11.93 * DEG2RAD,
            sigma_i: 0.01 * DEG2RAD,
            varpi: 311.5 * DEG2RAD,
            omega_big: 144.5 * DEG2RAD,
        },
        ObservedTno {
            name: "2012 VP113",
            a: 261.0,
            sigma_a: 8.0,
            e: 0.693,
            sigma_e: 0.002,
            i: 24.05 * DEG2RAD,
            sigma_i: 0.01 * DEG2RAD,
            varpi: 293.8 * DEG2RAD,
            omega_big: 90.8 * DEG2RAD,
        },
        ObservedTno {
            name: "2015 TG387",
            a: 1094.0,
            sigma_a: 75.0,
            e: 0.945,
            sigma_e: 0.003,
            i: 11.67 * DEG2RAD,
            sigma_i: 0.01 * DEG2RAD,
            varpi: 118.2 * DEG2RAD,
            omega_big: 300.8 * DEG2RAD,
        },
        ObservedTno {
            name: "2013 SY99",
            a: 735.0,
            sigma_a: 50.0,
            e: 0.932,
            sigma_e: 0.003,
            i: 4.23 * DEG2RAD,
            sigma_i: 0.01 * DEG2RAD,
            varpi: 32.5 * DEG2RAD,
            omega_big: 29.5 * DEG2RAD,
        },
        ObservedTno {
            name: "2014 SR349",
            a: 289.0,
            sigma_a: 15.0,
            e: 0.840,
            sigma_e: 0.005,
            i: 18.0 * DEG2RAD,
            sigma_i: 0.02 * DEG2RAD,
            varpi: 341.5 * DEG2RAD,
            omega_big: 34.8 * DEG2RAD,
        },
        ObservedTno {
            name: "2004 VN112",
            a: 319.0,
            sigma_a: 4.0,
            e: 0.854,
            sigma_e: 0.001,
            i: 25.56 * DEG2RAD,
            sigma_i: 0.01 * DEG2RAD,
            varpi: 327.1 * DEG2RAD,
            omega_big: 66.0 * DEG2RAD,
        },
        ObservedTno {
            name: "2010 GB174",
            a: 351.0,
            sigma_a: 9.0,
            e: 0.862,
            sigma_e: 0.002,
            i: 21.54 * DEG2RAD,
            sigma_i: 0.01 * DEG2RAD,
            varpi: 347.8 * DEG2RAD,
            omega_big: 130.6 * DEG2RAD,
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn test_sample_count() {
        let sample = distant_tno_sample();
        assert!(sample.len() >= 7, "should have at least 7 TNOs");
    }

    #[test]
    fn test_all_pass_selection() {
        for tno in distant_tno_sample() {
            assert!(
                passes_selection(&tno),
                "{} should pass selection criteria",
                tno.name
            );
        }
    }

    #[test]
    fn test_generate_clones_count() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let tno = &distant_tno_sample()[0];
        let clones = generate_clones(tno, 25, &mut rng);
        assert_eq!(clones.len(), 25);
    }

    #[test]
    fn test_clones_near_parent() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let tno = &distant_tno_sample()[0]; // Sedna
        let clones = generate_clones(tno, 25, &mut rng);

        for c in &clones {
            assert!(
                (c.a - tno.a).abs() < 5.0 * tno.sigma_a,
                "clone a = {} too far from parent a = {}",
                c.a,
                tno.a
            );
        }
    }

    #[test]
    fn test_clone_eccentricity_valid() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(99);
        let tno = &distant_tno_sample()[0];
        let clones = generate_clones(tno, 100, &mut rng);
        for c in &clones {
            assert!(c.e > 0.0 && c.e < 1.0, "e = {} invalid", c.e);
        }
    }
}
