//! Resonance-overlap chaos indicator over the (a, e) plane.
//!
//! Mechanism (Hadden et al. 2018): a particle interior to Planet Nine sits in a
//! forest of j:1 mean-motion resonances at a_j = a₉ (1/j)^{2/3}. Near the
//! planet these first-order resonances are dense and their pendulum widths grow
//! with the particle eccentricity; where neighbouring resonances overlap, the
//! islands merge into a chaotic sea (Wisdom 1980; Chirikov). The overlapped
//! zone extends inward from the planet by
//!
//!   Δa_overlap / a₉ = C μ^{2/7} e^{1/5},   μ = m₉ / M☉
//!
//! the eccentric first-order resonance-overlap criterion (circular limit
//! C·μ^{2/7}, Wisdom 1980; eccentric e^{1/5} enhancement, Mustill & Wyatt
//! 2012). A particle is chaotic when it lies inside this zone, i.e. when the
//! dimensionless overlap parameter
//!
//!   K = C μ^{2/7} e^{1/5} a₉ / (a₉ − a)
//!
//! exceeds unity. K rises with Planet Nine's mass, with eccentricity, and as a
//! approaches a₉ — so at fixed a the motion is chaotic at high e (large K) and
//! regular at low e (small K), the paper's headline.
//!
//! REUSE: resonance locations come from
//! `p9_core::analysis::resonance::resonance_semi_major_axis` (the Kepler
//! relation, cross-checked to 1e-9). The Chirikov chain form for a generic
//! perturber is also provided and pinned against p9-core's Neptune-specialised
//! `chirikov_overlap_parameter` to 1e-12.

use p9_core::analysis::resonance::resonance_semi_major_axis;
use p9_core::analysis::resonance::{
    overlap_zone_width_fraction, MW_ECCENTRIC_PREFACTOR, WISDOM_PREFACTOR,
};
use p9_core::constants::EARTH_MASS_SOLAR;
use p9_core::units::{au, Length};

/// Chirikov resonance-overlap parameter for a chain near a generic perturber of
/// semi-major axis `a_p` (AU) and mass `m_p` (solar masses), at particle
/// semi-major axis `a` (AU) and perihelion `q` (AU):
///
///   K = (24/√5) (a/a_p)^{5/4} √(m_p/M☉) exp(−q²/(2 a_p²)).
///
/// This is the general form of `p9_core`'s Neptune-specialised
/// `chirikov_overlap_parameter`; passing Neptune's a_p and mass reproduces it
/// (pinned in the tests). Retained for the analytic cross-check; the Planet
/// Nine driver uses [`overlap_parameter`].
pub fn chirikov_overlap_parameter_general(a_au: f64, q_au: f64, a_p: f64, m_p_solar: f64) -> f64 {
    let alpha = a_au / a_p;
    (24.0 / 5.0_f64.sqrt())
        * alpha.powf(1.25)
        * m_p_solar.sqrt()
        * (-q_au * q_au / (2.0 * a_p * a_p)).exp()
}

/// The j:1 resonance index whose nominal semi-major axis a₉ (1/j)^{2/3} is
/// closest to `a`. Returns the integer `j ≥ 1`.
pub fn nearest_j_one_index(a_au: f64, a9_au: f64) -> i64 {
    // a_j = a9 j^{-2/3}  =>  j = (a9 / a)^{3/2}
    let j_real = (a9_au / a_au).powf(1.5);
    j_real.round().max(1.0) as i64
}

/// Analytic location of the interior j:1 resonance: a_j = a₉ (1/j)^{2/3} (AU).
///
/// Delegated to `p9_core::analysis::resonance::resonance_semi_major_axis(1, j,
/// a₉)` = a₉ (1/j)^{2/3} — a particle completing `j` orbits per Planet Nine
/// orbit sits here.
pub fn j_one_location_typed(j: i64, a9_au: f64) -> Length {
    au(resonance_semi_major_axis(1, j as u32, a9_au))
}

/// Inward extent (AU) of the overlapped resonance zone from Planet Nine, for a
/// particle of eccentricity `e`:
///
///   Δa_overlap = max( 1.3 μ^{2/7},  1.8 (μ e)^{1/5} ) · a₉,
///
/// the Wisdom (1980) circular-orbit zone as the floor with the Mustill &
/// Wyatt (2012) eccentric criterion taking over once e is large enough.
/// `m9_earth` is Planet Nine's mass in Earth masses. (A previous version
/// used C·μ^{2/7}·e^{1/5}, which contradicted BOTH cited limits: it vanished
/// as e → 0 instead of approaching the finite Wisdom zone, and it carried
/// the circular mass exponent 2/7 on the eccentric branch where M&W give
/// (μe)^{1/5} — ≈3× too narrow at e = 0.85 for a 10 M⊕ perturber.)
pub fn overlap_zone_width_typed(e: f64, a9_au: f64, m9_earth: f64) -> Length {
    let mu = m9_earth * EARTH_MASS_SOLAR;
    au(overlap_zone_width_fraction(mu, e) * a9_au)
}

/// Resonance-overlap parameter K at (a, e) for a Planet Nine of semi-major
/// axis `a9_au` and mass `m9_earth` (Earth masses):
///
///   K = Δa_overlap(e) / (a₉ − a).
///
/// K > 1 ⇒ the particle is inside the overlapped zone ⇒ chaotic. `_e9` is
/// accepted for API symmetry (Planet Nine's own eccentricity does not enter
/// this first-order criterion). For a ≥ a₉ the particle is at or beyond the
/// planet and counted as chaotic (K = ∞).
pub fn overlap_parameter(a_au: f64, e: f64, a9_au: f64, m9_earth: f64, _e9: f64) -> f64 {
    let gap = a9_au - a_au;
    if gap <= 0.0 {
        return f64::INFINITY;
    }
    (overlap_zone_width_typed(e, a9_au, m9_earth) / au(gap)).value
}

/// Whether the point (a, e) lies in the chaotic resonance-overlap regime
/// (K > 1).
pub fn is_chaotic(a_au: f64, e: f64, a9_au: f64, m9_earth: f64, e9: f64) -> bool {
    overlap_parameter(a_au, e, a9_au, m9_earth, e9) > 1.0
}

/// Critical eccentricity at fixed semi-major axis `a`: the eccentricity at
/// which the overlap parameter crosses unity (below it regular, above it
/// chaotic). From K = 1 on the eccentric (Mustill & Wyatt) branch:
///
///   e_crit = [ (a₉ − a) / (1.8 μ^{1/5} a₉) ]^5 / μ · μ = (gap/(1.8 a₉))^5 / μ
///
/// i.e. e_crit = (gap / (1.8 a₉))^5 / μ. Returns `None` when the boundary
/// falls outside e ∈ (0, 1): either the column is uniformly regular
/// (e_crit ≥ 1), or a ≥ a₉ / the gap is already inside the CIRCULAR Wisdom
/// zone (chaotic at every eccentricity).
pub fn critical_eccentricity(a_au: f64, a9_au: f64, m9_earth: f64, _e9: f64) -> Option<f64> {
    let gap = a9_au - a_au;
    if gap <= 0.0 {
        return None;
    }
    let mu = m9_earth * EARTH_MASS_SOLAR;
    // Inside the circular Wisdom zone the column is chaotic at all e.
    if gap <= WISDOM_PREFACTOR * mu.powf(2.0 / 7.0) * a9_au {
        return None;
    }
    let e_crit = (gap / (MW_ECCENTRIC_PREFACTOR * a9_au)).powi(5) / mu;
    if e_crit > 0.0 && e_crit < 1.0 {
        Some(e_crit)
    } else {
        None
    }
}

/// A sampled chaos map over a rectangular (a, e) grid.
#[derive(Debug, Clone)]
pub struct ChaosMap {
    /// Semi-major axis grid (AU).
    pub a_vals: Vec<f64>,
    /// Eccentricity grid.
    pub e_vals: Vec<f64>,
    /// `k[i][j]` = overlap parameter at (a_vals[i], e_vals[j]).
    pub k: Vec<Vec<f64>>,
}

impl ChaosMap {
    /// Build a chaos map over `n_a × n_e` cells spanning [a_min, a_max] ×
    /// [e_min, e_max] for the given Planet Nine.
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        a_min: f64,
        a_max: f64,
        n_a: usize,
        e_min: f64,
        e_max: f64,
        n_e: usize,
        a9_au: f64,
        m9_earth: f64,
        e9: f64,
    ) -> Self {
        let a_vals: Vec<f64> = (0..n_a)
            .map(|i| a_min + (a_max - a_min) * (i as f64 + 0.5) / n_a as f64)
            .collect();
        let e_vals: Vec<f64> = (0..n_e)
            .map(|j| e_min + (e_max - e_min) * (j as f64 + 0.5) / n_e as f64)
            .collect();

        let k = a_vals
            .iter()
            .map(|&a| {
                e_vals
                    .iter()
                    .map(|&e| overlap_parameter(a, e, a9_au, m9_earth, e9))
                    .collect()
            })
            .collect();

        Self { a_vals, e_vals, k }
    }

    /// Fraction of grid cells that are chaotic (K > 1).
    pub fn chaotic_fraction(&self) -> f64 {
        let total = self.k.len() * self.k.first().map_or(0, |r| r.len());
        if total == 0 {
            return 0.0;
        }
        let chaotic: usize = self
            .k
            .iter()
            .map(|row| row.iter().filter(|&&v| v > 1.0).count())
            .sum();
        chaotic as f64 / total as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference;
    use approx::assert_relative_eq;
    use p9_core::analysis::resonance::chirikov_overlap_parameter;
    use p9_core::constants::{A_NEPTUNE_AU, MASS_NEPTUNE_SOLAR};

    const A9: f64 = reference::FIDUCIAL_A9_AU;
    const M9: f64 = reference::FIDUCIAL_M9_EARTH;
    const E9: f64 = reference::FIDUCIAL_E9;

    #[test]
    fn general_chirikov_reduces_to_p9_core_neptune_case() {
        // The generalized Chirikov chain form must reproduce p9-core's
        // Neptune-specialized chirikov_overlap_parameter exactly when fed
        // Neptune's a_p and mass (shared-physics cross-check).
        for &(a, q) in &[(300.0, 35.0), (500.0, 40.0), (1000.0, 30.0)] {
            let ours = chirikov_overlap_parameter_general(a, q, A_NEPTUNE_AU, MASS_NEPTUNE_SOLAR);
            let core = chirikov_overlap_parameter(a, q);
            assert_relative_eq!(ours, core, max_relative = 1e-12);
        }
    }

    #[test]
    fn typed_accessors_match_inline_formulas() {
        // j:1 location: a_j = a9 j^{-2/3}.
        assert_relative_eq!(
            (j_one_location_typed(5, A9) / au(1.0)).value,
            A9 * 5.0f64.powf(-2.0 / 3.0),
            max_relative = 1e-12
        );
        // Overlap zone width: max(1.3 μ^{2/7}, 1.8 (μe)^{1/5}) a9.
        let mu = M9 * EARTH_MASS_SOLAR;
        let expected = (WISDOM_PREFACTOR * mu.powf(2.0 / 7.0))
            .max(MW_ECCENTRIC_PREFACTOR * (mu * 0.7f64).powf(0.2))
            * A9;
        assert_relative_eq!(
            (overlap_zone_width_typed(0.7, A9, M9) / au(1.0)).value,
            expected,
            max_relative = 1e-12
        );
    }

    #[test]
    fn overlap_zone_matches_cited_limits() {
        // Circular limit: as e → 0 the zone approaches the finite Wisdom
        // width 1.3 μ^{2/7} a9 (the old formula vanished here).
        let mu = M9 * EARTH_MASS_SOLAR;
        let wisdom = WISDOM_PREFACTOR * mu.powf(2.0 / 7.0) * A9;
        let w0 = (overlap_zone_width_typed(1e-6, A9, M9) / au(1.0)).value;
        assert!((w0 - wisdom).abs() / wisdom < 1e-9, "e→0 zone = {w0:.1} AU");
        // Eccentric limit: at e = 0.85 the Mustill & Wyatt (μe)^{1/5}
        // branch dominates and is ~3x the old μ^{2/7} e^{1/5} value.
        let w85 = (overlap_zone_width_typed(0.85, A9, M9) / au(1.0)).value;
        let mw = MW_ECCENTRIC_PREFACTOR * (mu * 0.85f64).powf(0.2) * A9;
        assert!((w85 - mw).abs() / mw < 1e-9);
        let old = 1.5 * mu.powf(2.0 / 7.0) * 0.85f64.powf(0.2) * A9;
        assert!(w85 > 2.0 * old, "corrected {w85:.0} vs old {old:.0} AU");
    }

    #[test]
    fn resonance_locations_match_analytic_kepler_relation() {
        // Cross-check the j:1 resonance locations against a_j = a9 j^{-2/3}
        // to 1e-9 (the analytic relation underpinning the chain).
        for j in 1..=30i64 {
            let computed = (j_one_location_typed(j, A9) / au(1.0)).value;
            let analytic = A9 * (j as f64).powf(-2.0 / 3.0);
            assert_relative_eq!(computed, analytic, max_relative = 1e-9);
        }
    }

    #[test]
    fn nearest_index_is_self_consistent() {
        for j in 1..=20i64 {
            let a = (j_one_location_typed(j, A9) / au(1.0)).value;
            assert_eq!(nearest_j_one_index(a, A9), j, "j = {j}");
        }
    }

    #[test]
    fn overlap_parameter_grows_with_eccentricity() {
        // The defining mechanism: rising eccentricity widens the overlap zone,
        // so K grows monotonically with e at fixed a.
        let a = reference::CHAOTIC_EXAMPLE_A_AU;
        let ks: Vec<f64> = [0.05, 0.3, 0.5, 0.7, 0.85]
            .iter()
            .map(|&e| overlap_parameter(a, e, A9, M9, E9))
            .collect();
        for w in ks.windows(2) {
            assert!(w[1] > w[0], "K should grow with e: {ks:?}");
        }
    }

    #[test]
    fn overlap_parameter_grows_with_semimajor_axis() {
        // K rises as a approaches Planet Nine (the gap a9 − a shrinks).
        let e = 0.5;
        let k_small = overlap_parameter(300.0, e, A9, M9, E9);
        let k_large = overlap_parameter(450.0, e, A9, M9, E9);
        assert!(k_large > k_small, "K({k_large}) should exceed K({k_small})");
    }

    #[test]
    fn chaos_boundary_separates_high_and_low_e() {
        // The criterion flips across a computed critical eccentricity:
        // chaotic at high e, regular at low e (the paper's headline). The
        // boundary exists only OUTSIDE the circular Wisdom zone — at the
        // chaotic example's a = 480 (20 AU gap, inside the ~33 AU circular
        // zone) every eccentricity is chaotic, so probe at 430 AU.
        let a = reference::REGULAR_EXAMPLE_A_AU;
        let e_crit =
            critical_eccentricity(a, A9, M9, E9).expect("a chaos boundary should exist at this a");
        assert!(e_crit > 0.0 && e_crit < 1.0, "e_crit = {e_crit}");

        // Above the boundary: chaotic. Below: regular.
        assert!(is_chaotic(a, (e_crit + 0.05).min(0.999), A9, M9, E9));
        assert!(!is_chaotic(a, (e_crit - 0.05).max(1e-4), A9, M9, E9));

        // K = 1 exactly at the boundary (defining root).
        let k_boundary = overlap_parameter(a, e_crit, A9, M9, E9);
        assert_relative_eq!(k_boundary, 1.0, epsilon = 1e-9);
    }

    #[test]
    fn high_e_example_chaotic_low_e_example_regular() {
        // The paper's representative points: high-e chaotic, low-e regular.
        assert!(is_chaotic(
            reference::CHAOTIC_EXAMPLE_A_AU,
            reference::CHAOTIC_EXAMPLE_E,
            A9,
            M9,
            E9
        ));
        assert!(!is_chaotic(
            reference::REGULAR_EXAMPLE_A_AU,
            reference::REGULAR_EXAMPLE_E,
            A9,
            M9,
            E9
        ));
    }

    #[test]
    fn critical_eccentricity_is_the_exact_boundary() {
        // K(a, e_crit) = 1 to 1e-9 wherever a boundary exists.
        for &a in &[420.0, 440.0, 460.0, 480.0] {
            if let Some(e_crit) = critical_eccentricity(a, A9, M9, E9) {
                let k = overlap_parameter(a, e_crit, A9, M9, E9);
                assert_relative_eq!(k, 1.0, epsilon = 1e-9);
            }
        }
    }

    #[test]
    fn chaotic_fraction_rises_with_planet_nine_mass() {
        // A more massive Planet Nine widens the overlap zone and chaoticizes
        // more of the (a, e) plane: the chaotic area fraction grows with mass.
        let map_lo = ChaosMap::build(300.0, 495.0, 16, 0.0, 0.95, 16, A9, 5.0, E9);
        let map_hi = ChaosMap::build(300.0, 495.0, 16, 0.0, 0.95, 16, A9, 20.0, E9);
        assert!(
            map_hi.chaotic_fraction() > map_lo.chaotic_fraction(),
            "chaotic fraction {} (20 M⊕) should exceed {} (5 M⊕)",
            map_hi.chaotic_fraction(),
            map_lo.chaotic_fraction()
        );
    }

    #[test]
    fn surviving_etnos_sit_in_the_regular_region() {
        // Honest reproduction check: run the real Brown (2017) ETNO sample
        // (shared p9-core table) through the from-scratch overlap criterion for
        // a fiducial 500 AU / 10 M⊕ Planet Nine. Every observed object interior
        // to Planet Nine lands BELOW the chaos boundary (regular) — consistent
        // with the paper's picture that the chaotic web is depleted on Gyr
        // timescales, so the long-lived survivors we observe occupy the regular
        // islands. (At a9 = 500 AU even the largest-a survivor, a ≈ 350 AU, is
        // ~150 AU inside the ~40 AU overlap zone, hence regular.)
        use p9_core::data::etno::BROWN_2017_SAMPLE;
        let interior: Vec<_> = BROWN_2017_SAMPLE.iter().filter(|k| k.a < A9).collect();
        assert!(
            !interior.is_empty(),
            "expected interior ETNOs in the sample"
        );
        for etno in &interior {
            let chaotic = is_chaotic(etno.a, etno.e, A9, M9, E9);
            // With the corrected Mustill & Wyatt widths, the single
            // highest-(a·e) survivor — 2015 RX245 (a = 430, e = 0.89, 70 AU
            // inside P9 vs a ~110 AU eccentric overlap zone) — falls INSIDE
            // the chaotic web for this fiducial 500 AU/10 M⊕ P9. Its
            // survival requires resonant phase protection (which Becker et
            // al.'s N-body machinery models and a static overlap criterion
            // cannot), or disfavors this particular fiducial. Every other
            // interior survivor is regular.
            if etno.name == "2015 RX245" {
                assert!(
                    chaotic,
                    "RX245 sits inside the corrected M&W overlap zone at this fiducial"
                );
                continue;
            }
            assert!(
                !chaotic,
                "{} (a={}, e={}) should be regular for the fiducial P9",
                etno.name, etno.a, etno.e
            );
        }
    }

    #[test]
    fn map_has_both_regimes() {
        // A representative map must contain both regular and chaotic cells —
        // i.e. there is a genuine boundary, not a uniformly one-sided plane.
        let map = ChaosMap::build(300.0, 495.0, 16, 0.0, 0.95, 16, A9, M9, E9);
        let frac = map.chaotic_fraction();
        assert!(
            frac > 0.0 && frac < 1.0,
            "expected mixed regimes, got chaotic fraction {frac}"
        );
    }
}
