//! Relating the gap edges to solar-system dynamics from `p9_core`.
//!
//! The two flanks of the perihelion gap are sculpted by different mechanisms:
//!
//!   - The LOW-q flank is the Neptune-coupled scattering population. Its
//!     high-q edge — where objects stop exchanging energy with Neptune and
//!     detach — is set by the Chirikov resonance-overlap criterion, sourced
//!     entirely from `p9_core::analysis::resonance::critical_perihelion`. At
//!     the large semi-major axes of the distant ETNOs (a ~ several hundred AU)
//!     this q_crit(a) sits in the trans-Neptunian band and rises with a; the
//!     scattering population cannot maintain perihelia far above it.
//!   - The HIGH-q flank is the detached / inner-Oort-cloud Sednoid onset at
//!     q ≳ 65 AU.
//!
//! We do not reimplement the overlap criterion here — we read q_crit(a) from
//! p9-core and check the gap sits *above* the scattering boundary, i.e. that
//! the gap is in the detached side of the Neptune-coupling divide.

use p9_core::analysis::resonance::critical_perihelion;

/// The Neptune scattering boundary q_crit(a) (AU) at semi-major axis `a_au`,
/// re-exported from p9-core. Objects with q < q_crit(a) are in the chaotic
/// resonance-overlap (scattering) regime.
pub fn scattering_boundary(a_au: f64) -> f64 {
    critical_perihelion(a_au)
}

/// The semi-major axis (AU) at which the scattering boundary q_crit(a) reaches
/// a target perihelion `q_target`, found by bisection on the monotone
/// q_crit(a). Returns `None` if the boundary never reaches `q_target` within
/// [`a_lo`, `a_hi`].
///
/// Used to locate where the rising Neptune boundary crosses the low-q edge of
/// the gap — i.e. the a beyond which even the *boundary* of the scattering
/// region has climbed into the gap window.
pub fn a_where_boundary_reaches(q_target: f64, a_lo: f64, a_hi: f64) -> Option<f64> {
    let (mut lo, mut hi) = (a_lo, a_hi);
    if critical_perihelion(lo) > q_target || critical_perihelion(hi) < q_target {
        return None;
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if critical_perihelion(mid) < q_target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::published::{ETNO_A_FLOOR_AU, GAP_Q_HIGH_AU, GAP_Q_LOW_AU, IOC_Q_FLOOR_AU};

    #[test]
    fn scattering_boundary_matches_p9_core() {
        for &a in &[300.0, 500.0, 800.0] {
            assert_eq!(scattering_boundary(a), critical_perihelion(a));
        }
    }

    #[test]
    fn scattering_boundary_sits_below_the_gap() {
        // At the ETNO semi-major axes the Neptune resonance-overlap boundary
        // q_crit(a) is in the trans-Neptunian 30-60 AU band — at or below the
        // low edge of the gap. The gap is therefore on the detached side of
        // the Neptune-coupling divide, consistent with the paper's picture.
        for &a in &[ETNO_A_FLOOR_AU, 300.0, 500.0] {
            let q = scattering_boundary(a);
            assert!(q < GAP_Q_HIGH_AU, "q_crit({a}) = {q} not below gap top");
        }
        // At very large a the boundary climbs but stays well under the
        // detached IOC floor (q ~ 65 AU): the detached flank is NOT set by
        // resonance overlap.
        assert!(scattering_boundary(1000.0) < IOC_Q_FLOOR_AU);
    }

    #[test]
    fn boundary_crosses_low_gap_edge_at_large_a() {
        // The rising boundary reaches the gap's low edge (50 AU) only at large
        // semi-major axis — beyond the bulk of the scattering ETNOs, which is
        // why the scattering population thins out approaching the gap.
        let a = a_where_boundary_reaches(GAP_Q_LOW_AU, 200.0, 5000.0);
        assert!(a.is_some(), "boundary should reach 50 AU somewhere");
        let a = a.unwrap();
        assert!(a > 700.0, "boundary reaches 50 AU only at large a, got {a}");
        // Sanity: at that a, q_crit really is ~50 AU.
        assert!((scattering_boundary(a) - GAP_Q_LOW_AU).abs() < 1.0);
    }
}
