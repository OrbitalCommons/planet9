//! Mass–distance exclusion map: which `(mass, heliocentric-distance)` Planet
//! Nines does the Cassini range precision rule out?
//!
//! Holman & Payne (2016) confront the induced Earth–Saturn range signal with
//! the Cassini range precision (order tens of metres). A P9 whose peak range
//! residual over the orbit exceeds that precision would have shown up in the
//! data and is EXCLUDED; one below it is ALLOWED. Sweeping mass and distance
//! gives the paper's mass–distance constraint: close, massive planets are
//! excluded, distant or light ones survive.

use p9_core::types::P9Params;
use p9_core::units::{au, earth_masses, meters, Length, Mass};

use crate::published::CASSINI_RANGE_PRECISION_M;
use crate::signal::range_residual_m;

/// Number of true-anomaly samples used to find a P9's peak range residual.
const PERIHELION_ARC_SAMPLES: usize = 180;

/// The strongest range residual (METRES) P9 can produce anywhere on its orbit.
///
/// The signal is largest on the perihelion-facing arc (closest approach to
/// Saturn), so we scan true anomaly and take the maximum. This is the quantity
/// compared against the Cassini precision: if even the worst-case phase stays
/// below precision, the planet is unconstrained.
pub fn peak_range_residual_m(params: &P9Params) -> f64 {
    let mut peak = 0.0_f64;
    for k in 0..PERIHELION_ARC_SAMPLES {
        let nu = std::f64::consts::TAU * (k as f64) / (PERIHELION_ARC_SAMPLES as f64);
        peak = peak.max(range_residual_m(params, nu));
    }
    peak
}

/// The strongest range residual P9 can produce anywhere on its orbit, as a
/// typed [`Length`] — the [`peak_range_residual_m`] value in metres.
pub fn peak_range_residual(params: &P9Params) -> Length {
    meters(peak_range_residual_m(params))
}

/// Verdict for a candidate Planet Nine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExclusionVerdict {
    /// Peak range residual exceeds the Cassini precision — ruled out.
    Excluded,
    /// Peak range residual is below the Cassini precision — survives.
    Allowed,
}

/// Whether `params` is excluded by the Cassini range precision.
pub fn excluded(params: &P9Params) -> ExclusionVerdict {
    if peak_range_residual_m(params) > CASSINI_RANGE_PRECISION_M {
        ExclusionVerdict::Excluded
    } else {
        ExclusionVerdict::Allowed
    }
}

/// Largest P9 mass (Earth masses) that remains ALLOWED at the given semi-major
/// axis `a_au`, holding the other B&B orbital elements fixed.
///
/// Because the residual is exactly linear in mass, the boundary is found in
/// closed form from a single unit-mass evaluation: scale the peak residual at
/// `m = 1 M⊕` up to the Cassini precision. Returns the mass at which the peak
/// residual equals the precision floor.
pub fn max_allowed_mass_earth(template: &P9Params, a_au: f64) -> f64 {
    let mut p = *template;
    p.a = a_au;
    p.mass_earth = 1.0;
    let peak_per_earth_mass = peak_range_residual_m(&p);
    CASSINI_RANGE_PRECISION_M / peak_per_earth_mass
}

/// Largest allowed P9 mass at semi-major axis `a_au` as a typed [`Mass`] — the
/// [`max_allowed_mass_earth`] value in Earth masses.
pub fn max_allowed_mass(template: &P9Params, a_au: f64) -> Mass {
    earth_masses(max_allowed_mass_earth(template, a_au))
}

/// A sampled mass–distance exclusion map over a grid of semi-major axes.
#[derive(Debug, Clone)]
pub struct ExclusionMap {
    /// Sampled semi-major axes (AU).
    pub a_au: Vec<f64>,
    /// Largest allowed mass (Earth masses) at each `a` — the exclusion boundary.
    pub max_allowed_mass_earth: Vec<f64>,
}

impl ExclusionMap {
    /// Build the map over `n` semi-major axes uniform in `[a_min, a_max]`.
    pub fn build(template: &P9Params, a_min: f64, a_max: f64, n: usize) -> Self {
        assert!(n >= 2, "need at least two samples");
        let mut a_au = Vec::with_capacity(n);
        let mut max_allowed_mass_earth = Vec::with_capacity(n);
        for k in 0..n {
            let frac = (k as f64) / ((n - 1) as f64);
            let a = a_min + frac * (a_max - a_min);
            a_au.push(a);
            max_allowed_mass_earth.push(super::max_allowed_mass_earth(template, a));
        }
        Self {
            a_au,
            max_allowed_mass_earth,
        }
    }

    /// Sampled semi-major axes as typed [`Length`]s.
    pub fn semi_major_axes(&self) -> Vec<Length> {
        self.a_au.iter().map(|&a| au(a)).collect()
    }

    /// The exclusion-boundary masses as typed [`Mass`]es.
    pub fn max_allowed_masses(&self) -> Vec<Mass> {
        self.max_allowed_mass_earth
            .iter()
            .map(|&m| earth_masses(m))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::data::ephemeris_constraint::brown_batygin_orbit;

    #[test]
    fn close_massive_p9_is_excluded() {
        // A heavy P9 brought in close produces a residual far above precision.
        let mut p = brown_batygin_orbit();
        p.mass_earth = 20.0;
        p.a = 300.0;
        assert_eq!(excluded(&p), ExclusionVerdict::Excluded);
    }

    #[test]
    fn distant_light_p9_is_allowed() {
        // A light P9 pushed far out drops below the Cassini precision floor.
        let mut p = brown_batygin_orbit();
        p.mass_earth = 1.0;
        p.a = 2000.0;
        assert_eq!(excluded(&p), ExclusionVerdict::Allowed);
    }

    #[test]
    fn allowed_mass_boundary_is_consistent() {
        // At the boundary mass the peak residual equals precision; just above it
        // the planet flips to excluded.
        let template = brown_batygin_orbit();
        let a = 500.0;
        let m_max = max_allowed_mass_earth(&template, a);
        let mut at_boundary = template;
        at_boundary.a = a;
        at_boundary.mass_earth = m_max * 0.99;
        let mut over_boundary = template;
        over_boundary.a = a;
        over_boundary.mass_earth = m_max * 1.01;
        assert_eq!(excluded(&at_boundary), ExclusionVerdict::Allowed);
        assert_eq!(excluded(&over_boundary), ExclusionVerdict::Excluded);
    }

    #[test]
    fn exclusion_boundary_rises_with_distance() {
        // Farther planets need to be more massive to be detectable, so the
        // max-allowed-mass boundary increases monotonically with distance.
        let template = brown_batygin_orbit();
        let map = ExclusionMap::build(&template, 300.0, 1500.0, 12);
        for w in map.max_allowed_mass_earth.windows(2) {
            assert!(
                w[1] > w[0],
                "allowed-mass boundary not monotone in distance: {:?}",
                map.max_allowed_mass_earth
            );
        }
        // Steep distance dependence: the allowed mass should grow by a large
        // factor across a 5× distance span (tidal ≈ 1/d³).
        let lo = map.max_allowed_mass_earth.first().copied().unwrap();
        let hi = map.max_allowed_mass_earth.last().copied().unwrap();
        assert!(hi / lo > 10.0, "boundary growth too weak: {}", hi / lo);
    }

    #[test]
    fn typed_exclusion_accessors_match_f64() {
        use approx::assert_relative_eq;
        let template = brown_batygin_orbit();
        assert_relative_eq!(
            (peak_range_residual(&template) / meters(1.0)).value,
            peak_range_residual_m(&template),
            max_relative = 1e-12
        );
        assert_relative_eq!(
            (max_allowed_mass(&template, 500.0) / earth_masses(1.0)).value,
            max_allowed_mass_earth(&template, 500.0),
            max_relative = 1e-12
        );
        let map = ExclusionMap::build(&template, 300.0, 1500.0, 5);
        let axes = map.semi_major_axes();
        let masses = map.max_allowed_masses();
        for (k, (&a, &m)) in map
            .a_au
            .iter()
            .zip(map.max_allowed_mass_earth.iter())
            .enumerate()
        {
            assert_relative_eq!((axes[k] / au(1.0)).value, a, max_relative = 1e-12);
            assert_relative_eq!(
                (masses[k] / earth_masses(1.0)).value,
                m,
                max_relative = 1e-12
            );
        }
    }

    #[test]
    fn nominal_brown_batygin_p9_is_excluded() {
        // The 10 M⊕, a = 700 AU nominal B&B P9 sits where Cassini has leverage:
        // somewhere on its orbit it would exceed the range precision (which is
        // exactly why the data localize it rather than ignore it).
        let p = brown_batygin_orbit();
        assert_eq!(excluded(&p), ExclusionVerdict::Excluded);
    }
}
