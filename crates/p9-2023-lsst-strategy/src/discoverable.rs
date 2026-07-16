//! Discoverable fraction of a reference Planet Nine population under an LSST
//! observing strategy.
//!
//! For each synthetic Planet Nine realization (mass, distance, sky position,
//! reflected-light `V`), the per-orbit *discovery* probability is
//!
//!   P(discover) = footprint(δ, b) · coverage · P(link | V)
//!
//! where
//!   * `footprint(δ, b)` is the Rubin/LSST declination band and Galactic-plane
//!     cut (reusing `p9_core::analysis::surveys::Footprint::accepts`), with the
//!     equatorial declination and Galactic latitude of the current sky
//!     position computed by `p9_core::coords::sky` from the orbit's
//!     heliocentric ecliptic state;
//!   * `coverage` is the fraction of the in-band sky imaged to depth;
//!   * `P(link | V)` is the binomial-survival linking probability from
//!     [`crate::strategy::LsstStrategy::linking_probability`] (≥ N detections
//!     on independent nights at the single-visit recovery efficiency).
//!
//! The **discoverable fraction** is the mean of these per-orbit probabilities
//! over the population. It is, by construction, a probability in `(0, 1)`.
//! Nothing is hard-coded; the population is the seeded reference draw and the
//! strategy is the labelled-constant baseline (or a varied knob).

use serde::{Deserialize, Serialize};

use p9_core::analysis::surveys::Footprint;
use p9_core::constants::{DEG2RAD, GM_SUN, RAD2DEG};
use p9_core::coords::sky::{ecliptic_vec_to_equatorial_deg, equatorial_to_galactic};
use p9_core::types::OrbitalElements;

use crate::strategy::LsstStrategy;
use p9_core::analysis::photometry::SOLAR_V_MINUS_R;

/// One synthetic Planet Nine realization: its orbital elements (for the sky
/// position) and its reflected-light apparent `V` magnitude at that position.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct P9Sample {
    pub elements: OrbitalElements,
    pub v_magnitude: f64,
}

/// Result of a discoverable-fraction computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiscoverableResult {
    /// Mean per-orbit discovery probability over the population (in `(0, 1)`).
    pub fraction: f64,
    /// Population size evaluated.
    pub n_total: usize,
    /// Expected number of discovered objects (`fraction · n_total`).
    pub expected_discovered: f64,
    /// Fraction passing the footprint (geometry + coverage) gate alone.
    pub fraction_in_footprint: f64,
}

impl LsstStrategy {
    /// The footprint of this strategy as a reusable `p9_core` `Footprint`.
    pub fn footprint(&self) -> Footprint {
        Footprint {
            dec_min_deg: self.dec_min_deg,
            dec_max_deg: self.dec_max_deg,
            galactic_lat_min_deg: self.galactic_lat_min_deg,
            coverage_fraction: self.coverage_fraction,
        }
    }
}

/// Equatorial declination and Galactic latitude (degrees) of an orbit's
/// current heliocentric sky position. Parallax between heliocentric and
/// geocentric directions is < 0.2° at hundreds of AU and is neglected, as in
/// the other workspace footprint models.
pub fn dec_and_galactic_lat_deg(elem: &OrbitalElements) -> (f64, f64) {
    let state = elem.to_state_vector(GM_SUN);
    let (ra_deg, dec_deg) = ecliptic_vec_to_equatorial_deg(&state.pos);
    let (_, b) = equatorial_to_galactic(ra_deg * DEG2RAD, dec_deg * DEG2RAD);
    (dec_deg, b * RAD2DEG)
}

/// Per-orbit discovery probability under `strategy`.
pub fn discovery_probability(strategy: &LsstStrategy, sample: &P9Sample) -> f64 {
    let (dec_deg, gal_lat_deg) = dec_and_galactic_lat_deg(&sample.elements);
    let footprint = strategy.footprint();
    if !footprint.accepts(dec_deg, gal_lat_deg) {
        return 0.0;
    }
    footprint.coverage_fraction * strategy.linking_probability(r_magnitude(sample))
}

/// r-band apparent magnitude of a sample: the population carries V, the
/// LSST depths are r — convert with the labeled neutral-reflector color.
fn r_magnitude(sample: &P9Sample) -> f64 {
    sample.v_magnitude - SOLAR_V_MINUS_R
}

/// Discoverable fraction of `population` under `strategy`: the mean per-orbit
/// discovery probability (accumulated deterministically — no RNG).
pub fn discoverable_fraction(
    strategy: &LsstStrategy,
    population: &[P9Sample],
) -> DiscoverableResult {
    let n_total = population.len();
    if n_total == 0 {
        return DiscoverableResult {
            fraction: 0.0,
            n_total: 0,
            expected_discovered: 0.0,
            fraction_in_footprint: 0.0,
        };
    }
    let footprint = strategy.footprint();
    let mut expected = 0.0;
    let mut in_fp = 0.0;
    for s in population {
        let (dec_deg, gal_lat_deg) = dec_and_galactic_lat_deg(&s.elements);
        if footprint.accepts(dec_deg, gal_lat_deg) {
            in_fp += footprint.coverage_fraction;
            expected += footprint.coverage_fraction * strategy.linking_probability(r_magnitude(s));
        }
    }
    DiscoverableResult {
        fraction: expected / n_total as f64,
        n_total,
        expected_discovered: expected,
        fraction_in_footprint: in_fp / n_total as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn elem(i_deg: f64, omega_big_deg: f64, m_deg: f64) -> OrbitalElements {
        OrbitalElements {
            a: 460.0,
            e: 0.2,
            i: i_deg * DEG2RAD,
            omega: 150.0 * DEG2RAD,
            omega_big: omega_big_deg * DEG2RAD,
            mean_anomaly: m_deg * DEG2RAD,
        }
    }

    /// A handful of bright, southern, off-plane objects: every per-orbit
    /// probability is in (0, 1) and the mean is in (0, 1).
    #[test]
    fn discoverable_fraction_is_a_probability() {
        let strategy = LsstStrategy::baseline();
        let pop: Vec<P9Sample> = (0..60)
            .map(|k| P9Sample {
                elements: elem(20.0, 90.0, k as f64 * 6.0),
                v_magnitude: 22.0 + (k % 5) as f64 * 0.4,
            })
            .collect();
        let r = discoverable_fraction(&strategy, &pop);
        assert!(
            r.fraction > 0.0 && r.fraction < 1.0,
            "frac = {}",
            r.fraction
        );
        for s in &pop {
            let p = discovery_probability(&strategy, s);
            assert!((0.0..=1.0).contains(&p));
        }
    }

    /// A bright object inside the southern footprint is discovered with high
    /// probability; the same object pushed north of the dec_max is not.
    #[test]
    fn footprint_gate_zeroes_out_north_of_dec_max() {
        let strategy = LsstStrategy::baseline();
        // Find a southern in-footprint position.
        let south = (0..360)
            .map(|m| elem(15.0, 90.0, m as f64))
            .find(|e| {
                let (dec, b) = dec_and_galactic_lat_deg(e);
                dec < strategy.dec_max_deg - 10.0 && b.abs() > strategy.galactic_lat_min_deg
            })
            .expect("a southern off-plane position exists");
        let bright = P9Sample {
            elements: south,
            v_magnitude: 21.0,
        };
        assert!(discovery_probability(&strategy, &bright) > 0.5);

        // Shrinking dec_max below this object's declination removes it.
        let (dec, _) = dec_and_galactic_lat_deg(&bright.elements);
        let shrunk = strategy.clone().with_dec_max(dec - 5.0);
        assert_eq!(discovery_probability(&shrunk, &bright), 0.0);
    }

    #[test]
    fn empty_population_is_zero() {
        let r = discoverable_fraction(&LsstStrategy::baseline(), &[]);
        assert_relative_eq!(r.fraction, 0.0);
        assert_eq!(r.n_total, 0);
    }
}
