//! The scattering/detached divide in the (a, q) plane.
//!
//! An actively-scattering TNO has its perihelion deep enough in the Neptune
//! region that the overlapping 2:j mean-motion resonances drive a chaotic
//! random walk in semi-major axis (Kaib et al. 2019, OSSOS XV). As perihelion
//! rises, the resonances cease to overlap and the object *detaches*: it stops
//! exchanging energy with Neptune and freezes its semi-major axis.
//!
//! The divide is the Chirikov resonance-overlap criterion, sourced entirely
//! from `p9_core::analysis::resonance`: an object at (a, q) is actively
//! scattering iff the overlap parameter K(a, q) > 1, i.e. iff q < q_crit(a),
//! the K = 1 root. We do not reimplement any of that machinery here; this
//! module maps the boundary curve and classifies orbits with it.

use p9_core::analysis::resonance::{chirikov_overlap_parameter, critical_perihelion, is_chaotic};

/// Dynamical state of a TNO relative to the Neptune-scattering region.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DynamicalClass {
    /// Resonances overlap (K > 1): Neptune-coupled chaotic random walk in a.
    Scattering,
    /// Resonances do not overlap (K < 1): semi-major axis frozen.
    Detached,
}

/// Classify an orbit at semi-major axis `a_au` and perihelion `q_au`.
///
/// Thin wrapper over `p9_core::analysis::resonance::is_chaotic` — the single
/// source of the overlap criterion.
pub fn classify(a_au: f64, q_au: f64) -> DynamicalClass {
    if is_chaotic(a_au, q_au) {
        DynamicalClass::Scattering
    } else {
        DynamicalClass::Detached
    }
}

/// The critical perihelion q_crit(a) (AU) separating scattering from detached
/// at semi-major axis `a_au`. Re-exported from p9-core so callers map the
/// boundary through the same root that `classify` uses.
pub fn boundary_perihelion(a_au: f64) -> f64 {
    critical_perihelion(a_au)
}

/// A point on the scattering/detached boundary curve.
#[derive(Debug, Clone, Copy)]
pub struct BoundaryPoint {
    pub a_au: f64,
    pub q_crit_au: f64,
}

/// Sample the scattering/detached boundary q_crit(a) on a logarithmic grid in
/// semi-major axis from `a_min` to `a_max` (AU), `n` points inclusive.
///
/// Covers the OSSOS XV scattering range a ≈ 50–1000 AU. Points where
/// q_crit = 0 (no chaos at any q) are still returned so the caller sees the
/// full curve.
pub fn boundary_curve(a_min: f64, a_max: f64, n: usize) -> Vec<BoundaryPoint> {
    assert!(n >= 2 && a_min > 0.0 && a_max > a_min);
    let log_lo = a_min.ln();
    let log_hi = a_max.ln();
    (0..n)
        .map(|i| {
            let a = (log_lo + (log_hi - log_lo) * i as f64 / (n as f64 - 1.0)).exp();
            BoundaryPoint {
                a_au: a,
                q_crit_au: critical_perihelion(a),
            }
        })
        .collect()
}

/// Overlap parameter K(a, q), re-exported for direct inspection of how deep
/// into (or out of) the chaotic regime an orbit sits.
pub fn overlap(a_au: f64, q_au: f64) -> f64 {
    chirikov_overlap_parameter(a_au, q_au)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn boundary_matches_analytic_critical_perihelion() {
        // q_crit(a) here must be exactly p9-core's analytic relation, and
        // the overlap parameter must equal 1 there (cross-check at several a).
        // The 2:j chain only overlaps (q_crit > 0) above a ~ 226 AU; below
        // that the overlap argument < 1 and there is no chaos at any q.
        for &a in &[250.0, 500.0, 1000.0] {
            let q = boundary_perihelion(a);
            assert!(q > 0.0, "expected nonzero q_crit at a = {a}");
            assert_relative_eq!(q, critical_perihelion(a), epsilon = 1e-12);
            assert_relative_eq!(overlap(a, q), 1.0, epsilon = 1e-10);
        }
        // Below the overlap threshold the boundary collapses to q_crit = 0.
        assert_eq!(boundary_perihelion(100.0), 0.0);
    }

    #[test]
    fn scattering_region_is_chaotic_detached_is_not() {
        for &a in &[300.0, 500.0, 1000.0] {
            let q = boundary_perihelion(a);
            // Just below the boundary: scattering (chaotic).
            assert_eq!(classify(a, q - 1.0), DynamicalClass::Scattering);
            assert!(overlap(a, q - 1.0) > 1.0);
            // Just above: detached (stable).
            assert_eq!(classify(a, q + 1.0), DynamicalClass::Detached);
            assert!(overlap(a, q + 1.0) < 1.0);
        }
    }

    #[test]
    fn boundary_scale_near_published_neptune_region() {
        use crate::published::SCATTERING_Q_BOUNDARY_AU;
        // OSSOS XV: the scattering population detaches near the Neptune region,
        // q ~ 37 AU. The 2:j overlap boundary q_crit(a) rises monotonically
        // through this scale; it crosses the published ~37 AU somewhere in the
        // large-a scattering range and sits in the trans-Neptunian 30-60 AU
        // band there. (q_crit keeps climbing past 37 at the largest a, which is
        // the resolved physics, not a tuned cutoff.)
        let q400 = boundary_perihelion(400.0);
        let q500 = boundary_perihelion(500.0);
        assert!(
            q400 < SCATTERING_Q_BOUNDARY_AU && q500 > SCATTERING_Q_BOUNDARY_AU,
            "published 37 AU should fall between q_crit(400)={q400} and q_crit(500)={q500}"
        );
        for &a in &[300.0, 500.0, 800.0] {
            let q = boundary_perihelion(a);
            assert!((20.0..60.0).contains(&q), "q_crit({a}) = {q}");
        }
    }

    #[test]
    fn critical_perihelion_rises_with_semi_major_axis() {
        // Larger a => resonances are wider in the (a) direction relative to
        // their spacing, so chaos persists to higher q.
        let curve = boundary_curve(50.0, 1000.0, 40);
        let nonzero: Vec<f64> = curve
            .iter()
            .filter(|p| p.q_crit_au > 0.0)
            .map(|p| p.q_crit_au)
            .collect();
        // Roughly half the log-grid lies above the a ~ 226 AU overlap
        // threshold where q_crit becomes nonzero.
        assert!(nonzero.len() > 15, "only {} nonzero", nonzero.len());
        for w in nonzero.windows(2) {
            assert!(w[1] >= w[0] - 1e-9, "q_crit not monotone: {w:?}");
        }
    }

    #[test]
    fn boundary_curve_spans_requested_range() {
        let curve = boundary_curve(50.0, 1000.0, 10);
        assert_eq!(curve.len(), 10);
        assert_relative_eq!(curve[0].a_au, 50.0, epsilon = 1e-9);
        assert_relative_eq!(curve[9].a_au, 1000.0, epsilon = 1e-9);
    }
}
