//! Labelled reference constants from Khain et al. (2020), arXiv:2010.02234,
//! and the Planet Nine orbit they assume (Batygin & Brown nominal). These are
//! kept ONLY as documentation / cross-checks; nothing in the computation reads
//! them as an answer.

/// Planet Nine semimajor axis assumed by the resonance-hopping study (AU).
/// Khain et al. (2020) integrate distant TNOs against the nominal Batygin &
/// Brown Planet Nine; a₉ ≈ 500 AU is the canonical value of that era. The
/// resonance LOCATIONS scale with this, so it is exposed as the model input.
pub const A_P9_AU: f64 = 500.0;

/// Planet Nine eccentricity (nominal Batygin & Brown). Used to set the
/// perihelion distance of the synthetic distant population (objects detached
/// to comparable perihelia interact secularly/resonantly with P9).
pub const E_P9: f64 = 0.25;

/// Planet Nine perihelion q₉ = a₉(1 − e₉) ≈ 375 AU.
pub const Q_P9_AU: f64 = A_P9_AU * (1.0 - E_P9);

/// The paper's qualitative headline, encoded as a documented ordering rather
/// than a single pinned scalar (the paper reports classes per-object, not a
/// single population fraction at our reduced scale): resonance STICKING
/// dominates near low-order n:1 resonances at smaller a, and resonance
/// HOPPING dominates at larger a where the resonances overlap.
///
/// Anti-aligned objects survive without the resonances — the dynamics are
/// "predominantly driven by secular, rather than resonant, interactions"
/// (abstract). So a substantial NON-resonant / scattering class is expected
/// alongside the resonant and hopping classes.
pub const HOPPING_RISES_WITH_A: bool = true;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::population::{
        class_fractions, default_landscape, default_m_p9_solar, hopping_fraction_by_bin,
        synthetic_population,
    };

    #[test]
    fn perihelion_identity_holds() {
        assert!((Q_P9_AU - A_P9_AU * (1.0 - E_P9)).abs() < 1e-12);
        // Distant population perihelia bracket the giant-planet decoupling
        // scale and sit well interior to P9's perihelion.
        assert!(Q_P9_AU > 300.0 && Q_P9_AU < A_P9_AU);
    }

    #[test]
    fn headline_ordering_matches_published_claim() {
        // The paper's qualitative headline (encoded in HOPPING_RISES_WITH_A):
        // hopping rises with a. Confirm the computed bin trend agrees with the
        // documented reference flag, so the constant is not a stale label.
        let (a_min, a_max) = (150.0, 499.0);
        let l = default_landscape(a_min, a_max);
        let m = default_m_p9_solar();
        let pop = synthetic_population(2020, 6000, a_min, a_max, 30.0, 120.0);
        let bins = hopping_fraction_by_bin(&pop, &l, m, 1.2, a_min, a_max, 5);
        let rises = bins.last().unwrap().1 > bins.first().unwrap().1;
        assert_eq!(rises, HOPPING_RISES_WITH_A);

        // And a non-trivial non-resonant (secular/scattering) class exists.
        let f = class_fractions(&pop, &l, m, 1.2);
        assert!(
            f.non_resonant > 0.1,
            "expected a substantial secular/scattering class, got {:.3}",
            f.non_resonant
        );
    }
}
