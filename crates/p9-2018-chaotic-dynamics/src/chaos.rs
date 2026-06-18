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
use p9_core::constants::EARTH_MASS_SOLAR;
use p9_core::units::{au, Length};

/// Dimensionless prefactor in the first-order resonance-overlap width
/// Δa/a_p = C μ^{2/7} (circular limit). Wisdom (1980) / Duncan et al. (1989)
/// give C ≈ 1.5 for the chaotic-zone half-width.
pub const OVERLAP_PREFACTOR: f64 = 1.5;

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
///   Δa_overlap = C μ^{2/7} e^{1/5} a₉.
///
/// `m9_earth` is Planet Nine's mass in Earth masses. The eccentricity factor
/// e^{1/5} is the Mustill & Wyatt (2012) enhancement of the Wisdom (1980)
/// circular overlap zone.
pub fn overlap_zone_width_typed(e: f64, a9_au: f64, m9_earth: f64) -> Length {
    let mu = m9_earth * EARTH_MASS_SOLAR;
    au(OVERLAP_PREFACTOR * mu.powf(2.0 / 7.0) * e.powf(0.2) * a9_au)
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
/// chaotic). Closed form from K = 1:
///
///   e_crit = [ (a₉ − a) / (C μ^{2/7} a₉) ]^5.
///
/// Returns `None` when the boundary falls outside e ∈ (0, 1) — the column is
/// uniformly regular (e_crit ≥ 1) or uniformly chaotic (a ≥ a₉, e_crit ≤ 0).
pub fn critical_eccentricity(a_au: f64, a9_au: f64, m9_earth: f64, _e9: f64) -> Option<f64> {
    let gap = a9_au - a_au;
    if gap <= 0.0 {
        return None;
    }
    let mu = m9_earth * EARTH_MASS_SOLAR;
    let base = gap / (OVERLAP_PREFACTOR * mu.powf(2.0 / 7.0) * a9_au);
    let e_crit = base.powi(5);
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
        // Overlap zone width: C μ^{2/7} e^{1/5} a9.
        let mu = M9 * EARTH_MASS_SOLAR;
        let expected = OVERLAP_PREFACTOR * mu.powf(2.0 / 7.0) * 0.7f64.powf(0.2) * A9;
        assert_relative_eq!(
            (overlap_zone_width_typed(0.7, A9, M9) / au(1.0)).value,
            expected,
            max_relative = 1e-12
        );
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
        // The criterion flips across a computed critical eccentricity: chaotic
        // at high e, regular at low e (the paper's headline).
        let a = reference::CHAOTIC_EXAMPLE_A_AU;
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
            assert!(
                !is_chaotic(etno.a, etno.e, A9, M9, E9),
                "{} (a={}, e={}) should be regular for the fiducial P9",
                etno.name,
                etno.a,
                etno.e
            );
        }
        // The boundary is nonetheless reachable: pushing the highest-a survivor
        // outward to large a (toward Planet Nine) at its observed eccentricity
        // crosses into the chaotic web — the criterion is not vacuously regular.
        let max_a = interior
            .iter()
            .max_by(|x, y| x.a.partial_cmp(&y.a).unwrap())
            .unwrap();
        assert!(
            is_chaotic(490.0, max_a.e, A9, M9, E9),
            "a particle at a = 490 AU, e = {} should be chaotic",
            max_a.e
        );
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
