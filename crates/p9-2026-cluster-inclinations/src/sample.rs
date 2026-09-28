//! The observed high-perihelion sample: the 19 TNOs with q ∈ (40, 80) AU and
//! a ∈ (200, 2000) AU that the paper (Section 3.1, footnote) fits for the
//! intrinsic inclination width.
//!
//! Provenance: JPL Small-Body Database osculating elements at epoch
//! JD 2461200.5, fetched 2026-09-26, with each orbit's `first_obs` date (the
//! first observation of the fitted arc, which stands in for the discovery
//! epoch; where precovery was later linked it is earlier than the discovery
//! announcement). Angles are in degrees as served. Alicanto (2012 VP113)
//! sits at q = 80.6 AU in the current JPL solution, just outside the paper's
//! nominal q < 80 AU cut, but it is in the paper's list and is kept.
//!
//! The Brown (2001) debiasing needs, per object, the ecliptic latitude at
//! which it was discovered. That is not served directly, so it is
//! recomputed here from the object's own two-body orbit at the `first_obs`
//! epoch (heliocentric latitude; at ≥ 40 AU the geocentric correction is
//! < 1.5°).

use nalgebra::Vector3;
use p9_core::constants::{DEG2RAD, GM_SUN, TWO_PI};
use p9_core::data::refresh::parse_sbdb_date;
use p9_core::types::{elements_to_cartesian, OrbitalElements, StateVector};

/// Epoch of the pinned osculating elements (JD TDB).
pub const SAMPLE_EPOCH_JD: f64 = 2461200.5;

/// Paper's sample cuts (Section 3.1): q ∈ (40, 80) AU, a ∈ (200, 2000) AU,
/// all with i < 40°.
pub const Q_MIN_AU: f64 = 40.0;
pub const Q_MAX_AU: f64 = 80.0;
pub const A_MIN_AU: f64 = 200.0;
pub const A_MAX_AU: f64 = 2000.0;
pub const I_MAX_DEG: f64 = 40.0;

/// One object of the observed high-q sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HighQTno {
    pub name: &'static str,
    /// Semi-major axis (AU)
    pub a: f64,
    pub e: f64,
    /// Inclination (deg)
    pub i_deg: f64,
    /// Longitude of ascending node (deg)
    pub omega_big_deg: f64,
    /// Argument of perihelion (deg)
    pub omega_deg: f64,
    /// Mean anomaly at [`SAMPLE_EPOCH_JD`] (deg)
    pub mean_anomaly_deg: f64,
    /// SBDB `first_obs` date (YYYY-MM-DD)
    pub first_obs: &'static str,
}

#[allow(clippy::too_many_arguments)]
const fn tno(
    name: &'static str,
    a: f64,
    e: f64,
    i_deg: f64,
    omega_big_deg: f64,
    omega_deg: f64,
    mean_anomaly_deg: f64,
    first_obs: &'static str,
) -> HighQTno {
    HighQTno {
        name,
        a,
        e,
        i_deg,
        omega_big_deg,
        omega_deg,
        mean_anomaly_deg,
        first_obs,
    }
}

/// The 19 high-q TNOs of the paper's footnote, JPL SBDB elements.
pub const HIGH_Q_SAMPLE: [HighQTno; 19] = [
    tno(
        "2000 CR105",
        229.2563,
        0.807337,
        22.70799,
        128.20375,
        317.04090,
        6.31803,
        "2000-02-06",
    ),
    tno(
        "Alicanto (2012 VP113)",
        267.5210,
        0.698559,
        24.00179,
        90.89262,
        294.28523,
        3.83428,
        "2007-09-19",
    ),
    tno(
        "2013 UT15",
        205.3661,
        0.786232,
        10.65007,
        191.95588,
        252.08833,
        355.05344,
        "2004-09-21",
    ),
    tno(
        "Leleakuhonua (2015 TG387)",
        1345.9747,
        0.951945,
        11.67890,
        301.13467,
        118.23173,
        359.61901,
        "2005-10-05",
    ),
    tno(
        "2013 RA109",
        480.1107,
        0.904050,
        12.39043,
        104.89786,
        263.11434,
        0.64366,
        "2013-09-12",
    ),
    tno(
        "Sedna",
        543.7195,
        0.859882,
        11.92528,
        144.50617,
        311.09877,
        358.59569,
        "1990-09-25",
    ),
    tno(
        "2010 GB174",
        360.2421,
        0.865318,
        21.53976,
        130.57233,
        347.42133,
        3.92341,
        "2009-06-26",
    ),
    tno(
        "2013 FT28",
        288.3328,
        0.848932,
        17.39338,
        217.70895,
        40.51974,
        357.53196,
        "2013-03-16",
    ),
    tno(
        "2013 SY99",
        812.8892,
        0.938607,
        4.22556,
        29.51424,
        32.28018,
        359.55577,
        "2013-09-05",
    ),
    tno(
        "2014 SR349",
        304.1123,
        0.844164,
        17.96288,
        34.91767,
        341.66526,
        358.14682,
        "2014-09-19",
    ),
    tno(
        "2014 WB556",
        284.5402,
        0.849682,
        24.14837,
        115.00751,
        235.52627,
        2.02567,
        "2014-11-21",
    ),
    tno(
        "2015 KG163",
        640.1613,
        0.936762,
        14.00425,
        219.13302,
        31.94344,
        0.08953,
        "2015-05-17",
    ),
    tno(
        "2015 RX245",
        446.1028,
        0.897988,
        12.14273,
        8.61523,
        65.14447,
        358.49437,
        "2015-06-23",
    ),
    tno(
        "2016 SD106",
        357.4513,
        0.880582,
        4.80735,
        219.46157,
        163.08346,
        359.49247,
        "2013-09-08",
    ),
    tno(
        "2017 OF201",
        823.6372,
        0.945188,
        16.22096,
        328.76035,
        337.83543,
        1.45795,
        "2004-09-21",
    ),
    tno(
        "2018 VM35",
        306.9899,
        0.854731,
        8.47759,
        192.38284,
        303.35025,
        357.89975,
        "2018-11-06",
    ),
    tno(
        "2021 RR205",
        949.0345,
        0.941410,
        7.64896,
        108.45123,
        208.63449,
        0.42689,
        "2017-07-24",
    ),
    tno(
        "2023 KQ14",
        247.4285,
        0.733566,
        11.00096,
        72.07903,
        198.70375,
        356.57051,
        "2005-04-11",
    ),
    tno(
        "2024 FX26",
        253.7267,
        0.838749,
        15.66536,
        158.33936,
        34.54236,
        357.85351,
        "2013-04-01",
    ),
];

impl HighQTno {
    /// Osculating elements at [`SAMPLE_EPOCH_JD`] (radians / AU).
    pub fn elements(&self) -> OrbitalElements {
        self.elements_at(SAMPLE_EPOCH_JD)
    }

    /// Two-body elements at Julian date `jd`: only the mean anomaly moves,
    /// M(t) = M₀ + n (t − t₀).
    pub fn elements_at(&self, jd: f64) -> OrbitalElements {
        let n = (GM_SUN / self.a.powi(3)).sqrt();
        OrbitalElements {
            a: self.a,
            e: self.e,
            i: self.i_deg * DEG2RAD,
            omega_big: self.omega_big_deg * DEG2RAD,
            omega: self.omega_deg * DEG2RAD,
            mean_anomaly: (self.mean_anomaly_deg * DEG2RAD + n * (jd - SAMPLE_EPOCH_JD))
                .rem_euclid(TWO_PI),
        }
    }

    /// Perihelion distance q = a(1 − e) (AU).
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
    }

    /// Julian date of the first observation of the fitted arc.
    pub fn first_obs_jd(&self) -> f64 {
        parse_sbdb_date(self.first_obs).expect("pinned first_obs dates are well-formed")
    }

    /// Heliocentric ecliptic state vector at the first observation.
    pub fn discovery_state(&self) -> StateVector {
        elements_to_cartesian(&self.elements_at(self.first_obs_jd()), GM_SUN)
    }

    /// Heliocentric ecliptic latitude at the first observation (radians).
    pub fn discovery_latitude(&self) -> f64 {
        latitude(&self.discovery_state().pos)
    }
}

/// Latitude of a position vector relative to the frame's z-axis (radians).
pub fn latitude(pos: &Vector3<f64>) -> f64 {
    (pos.z / pos.norm()).clamp(-1.0, 1.0).asin()
}

/// Osculating elements of the whole sample.
pub fn sample_elements() -> Vec<OrbitalElements> {
    HIGH_Q_SAMPLE.iter().map(HighQTno::elements).collect()
}

/// Discovery-epoch heliocentric positions of the whole sample (AU).
pub fn discovery_positions() -> Vec<Vector3<f64>> {
    HIGH_Q_SAMPLE
        .iter()
        .map(|t| t.discovery_state().pos)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::RAD2DEG;

    #[test]
    fn sample_satisfies_paper_cuts() {
        for t in &HIGH_Q_SAMPLE {
            assert!(t.a > A_MIN_AU && t.a < A_MAX_AU, "{}: a = {}", t.name, t.a);
            // Alicanto's current JPL solution is 0.6 AU past the nominal cut.
            assert!(
                t.perihelion() > Q_MIN_AU && t.perihelion() < Q_MAX_AU + 1.0,
                "{}: q = {}",
                t.name,
                t.perihelion()
            );
            assert!(t.i_deg < I_MAX_DEG, "{}: i = {}", t.name, t.i_deg);
        }
    }

    #[test]
    fn discovery_latitude_never_exceeds_inclination() {
        // A body on an orbit of inclination i is confined to |β| ≤ i.
        for t in &HIGH_Q_SAMPLE {
            let beta = t.discovery_latitude().abs() * RAD2DEG;
            assert!(
                beta <= t.i_deg + 1e-6,
                "{}: |β| = {beta:.2}° > i = {}°",
                t.name,
                t.i_deg
            );
        }
    }

    #[test]
    fn objects_were_discovered_near_perihelion() {
        // Flux ∝ r⁻⁴: every object in this sample was first seen well inside
        // twice its perihelion distance.
        for t in &HIGH_Q_SAMPLE {
            let r = t.discovery_state().pos.norm();
            assert!(
                r < 2.0 * t.perihelion(),
                "{}: r_disc = {r:.1} AU vs q = {:.1}",
                t.name,
                t.perihelion()
            );
        }
    }
}
