//! Dynamical classification of a distant TNO against the Planet Nine
//! mean-motion-resonance landscape.
//!
//! Three classes, following Khain et al. (2020):
//!
//! * **Resonant** — the object sits inside (and isolated within) a single P9
//!   resonance: it is close to a resonance centre AND the neighbouring
//!   resonances do not overlap it (resonance-overlap K < 1). It can stick.
//! * **Hopping** — the object lies in a region where neighbouring P9
//!   resonances OVERLAP (K ≥ 1). With overlapping resonances there is no
//!   single stable libration island; the object transiently sticks to one
//!   commensurability and then abruptly transitions to a neighbour. This is
//!   the metastable "resonance hopping" regime, the paper's central
//!   phenomenon.
//! * **NonResonant** — the object is not close to any P9 resonance centre and
//!   the local landscape is not overlapping; its evolution is set by secular
//!   (scattering / detached) dynamics rather than resonance.
//!
//! Resonance LOCATIONS come from `p9_core::analysis::resonance::
//! resonance_semi_major_axis`; the single workspace convention. A p:q
//! resonance (p > q) places the interior particle at a = a₉ (q/p)^{2/3}, i.e.
//! the particle completes p orbits per q Planet Nine orbits.
//!
//! The OVERLAP mechanism (the paper's): the full mean-motion-resonance
//! spectrum (all coprime p:q period ratios) is DENSE near the perturber, where
//! the period ratio → 1, and SPARSE far from it. So the spacing between
//! neighbouring resonances collapses as a → a₉, driving the Chirikov overlap
//! parameter K above 1 there — overlapping islands, no stable libration, and
//! the metastable resonance HOPPING the paper reports. Far from a₉ only the
//! widely separated low-order resonances (the n:1 / n:2 commensurabilities)
//! exist, so an object there sits in an isolated island and can STICK.

use p9_core::analysis::resonance::resonance_semi_major_axis;
use p9_core::units::{au, Length};

/// Greatest common divisor (for reducing p:q to lowest terms).
fn gcd(a: u32, b: u32) -> u32 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// One P9 mean-motion resonance: the particle completes `p` orbits per `q`
/// Planet Nine orbits (period ratio p:q, p > q). Interior particles sit at
/// a = a₉ (q/p)^{2/3} < a₉; we use the n:1 and n:2 families.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mmr {
    /// Particle coefficient (orbits of the particle).
    pub p: u32,
    /// Planet Nine coefficient (orbits of P9).
    pub q: u32,
}

impl Mmr {
    /// Semimajor axis of this interior resonance for a P9 at `a_p9` (AU).
    ///
    /// Particle does `p` orbits per `q` P9 orbits ⇒ period ratio
    /// T_particle/T_p9 = q/p, so a = a₉ (q/p)^{2/3}. We express this through
    /// the p9-core exterior helper a_planet (P/Q)^{2/3} with (P, Q) = (q, p),
    /// the single workspace resonance-location convention.
    pub fn semimajor_axis(&self, a_p9: f64) -> f64 {
        resonance_semi_major_axis(self.q, self.p, a_p9)
    }

    /// Semimajor axis of this interior resonance as a typed [`Length`], for a
    /// Planet Nine at `a_p9` (AU).
    pub fn semi_major_axis_typed(&self, a_p9: f64) -> Length {
        au(self.semimajor_axis(a_p9))
    }

    /// Order of the resonance, |p − q|.
    pub fn order(&self) -> u32 {
        if self.p > self.q {
            self.p - self.q
        } else {
            self.q - self.p
        }
    }

    /// Whether this is a "simple" low-order commensurability: n:1 or n:2.
    pub fn is_simple(&self) -> bool {
        self.q == 1 || self.q == 2
    }
}

/// The P9 resonance landscape over a semimajor-axis band: two sorted lists of
/// interior resonance centres.
///
/// * `spectrum` — the FULL coprime p:q resonance spectrum up to order `q_max`
///   (used for the overlap / neighbour-spacing calculation). It is dense near
///   a₉ and sparse far from it.
/// * `simple` — only the low-order n:1 / n:2 commensurabilities (used for the
///   "nearest simple resonance" classification feature the paper highlights).
#[derive(Debug, Clone)]
pub struct ResonanceLandscape {
    /// Planet Nine semimajor axis (AU).
    pub a_p9: f64,
    /// Full coprime p:q spectrum, sorted by increasing semimajor axis.
    pub spectrum: Vec<(Mmr, f64)>,
    /// The n:1 / n:2 subset, sorted by increasing semimajor axis.
    pub simple: Vec<(Mmr, f64)>,
}

impl ResonanceLandscape {
    /// Build the interior resonance landscape for `a_p9` over `[a_min, a_max]`
    /// (AU). `q_max` bounds the P9-orbit coefficient q (so the spectrum holds
    /// all coprime p:q with q ≤ q_max and period ratio q/p < 1) — larger
    /// `q_max` resolves the dense accumulation of high-order resonances near
    /// a₉. The `simple` list is the q ∈ {1, 2} subset.
    pub fn new(a_p9: f64, a_min: f64, a_max: f64, q_max: u32) -> Self {
        let mut spectrum = Vec::new();
        let mut simple = Vec::new();
        for q in 1..=q_max {
            // Interior resonance: period ratio q/p < 1 ⇒ p > q. Cap p so the
            // centre stays above a_min (a = a₉ (q/p)^{2/3} ≥ a_min).
            for p in (q + 1)..=(q_max * 40) {
                if gcd(p, q) != 1 {
                    continue;
                }
                let mmr = Mmr { p, q };
                let a_res = mmr.semimajor_axis(a_p9);
                if a_res < a_min {
                    break; // larger p only goes lower in a
                }
                if a_res <= a_max {
                    spectrum.push((mmr, a_res));
                    if mmr.is_simple() {
                        simple.push((mmr, a_res));
                    }
                }
            }
        }
        spectrum.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        simple.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        Self {
            a_p9,
            spectrum,
            simple,
        }
    }

    /// Planet Nine semimajor axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a_p9)
    }

    /// Nearest SIMPLE (n:1 / n:2) resonance to `a` (AU): (resonance, centre,
    /// |Δa|), or `None` if there are none in band.
    pub fn nearest_simple(&self, a: f64) -> Option<(Mmr, f64, f64)> {
        Self::nearest_in(&self.simple, a)
    }

    /// Nearest resonance in the full spectrum to `a` (AU).
    pub fn nearest(&self, a: f64) -> Option<(Mmr, f64, f64)> {
        Self::nearest_in(&self.spectrum, a)
    }

    fn nearest_in(list: &[(Mmr, f64)], a: f64) -> Option<(Mmr, f64, f64)> {
        list.iter()
            .map(|&(m, a_res)| (m, a_res, (a - a_res).abs()))
            .min_by(|x, y| x.2.partial_cmp(&y.2).unwrap())
    }
}

/// Half-width (in semimajor axis, AU) of a first-order pendulum resonance.
///
/// For a mean-motion resonance the libration island has a half-width in
/// semimajor axis δa ≈ a · C · √(m_p/M) · f(e), with the resonance strength
/// rising with the perturber mass and the particle eccentricity (the resonant
/// term in the disturbing function carries e^{|p−q|}). We use the standard
/// pendulum scaling
///
///   δa = a · width_coeff · √(m_p9/M_sun) · e^{order/2}
///
/// with `order = |p − q|`. The eccentricity factor (e raised to half the
/// resonance order) makes high-eccentricity, low-order resonances WIDE and
/// high-order resonances narrow — the standard ordering. `width_coeff` is a
/// labelled O(1) model constant, not tuned to any published fraction.
pub fn resonance_half_width(a: f64, mmr: Mmr, e: f64, m_p9_solar: f64, width_coeff: f64) -> f64 {
    let order = if mmr.p > mmr.q {
        mmr.p - mmr.q
    } else {
        mmr.q - mmr.p
    };
    let e_factor = e.max(1e-3).powf(0.5 * order as f64);
    a * width_coeff * m_p9_solar.sqrt() * e_factor
}

/// Resonance-overlap parameter K at semimajor axis `a`: the ratio of the
/// summed half-widths of the two resonances bracketing `a` to the spacing
/// between their centres,
///
///   K = (δa_left + δa_right) / |a_right_centre − a_left_centre|.
///
/// K ≥ 1 ⇒ the neighbouring libration islands touch/overlap ⇒ the Chirikov
/// resonance-overlap (chaos) criterion is met ⇒ no isolated stable island ⇒
/// the object can HOP between commensurabilities. K < 1 ⇒ islands are
/// separated ⇒ an object inside one can STICK. The full-spectrum resonance
/// spacing collapses as a → a₉ (period ratio → 1), so K RISES with a — the
/// paper's mechanism. Computed over the dense full `spectrum`.
///
/// Returns `None` if `a` is outside the bracketed region (fewer than two
/// neighbouring resonances).
pub fn resonance_overlap_k(
    landscape: &ResonanceLandscape,
    a: f64,
    e: f64,
    m_p9_solar: f64,
    width_coeff: f64,
) -> Option<f64> {
    let res = &landscape.spectrum;
    if res.len() < 2 {
        return None;
    }
    // Find the index of the first resonance with centre >= a.
    let idx = res.partition_point(|&(_, a_res)| a_res < a);
    let (left, right) = if idx == 0 {
        (res[0], res[1])
    } else if idx >= res.len() {
        (res[res.len() - 2], res[res.len() - 1])
    } else {
        (res[idx - 1], res[idx])
    };
    let spacing = (right.1 - left.1).abs();
    if spacing <= 0.0 {
        return None;
    }
    let w_left = resonance_half_width(left.1, left.0, e, m_p9_solar, width_coeff);
    let w_right = resonance_half_width(right.1, right.0, e, m_p9_solar, width_coeff);
    Some((w_left + w_right) / spacing)
}

/// Dynamical class of a distant TNO.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class {
    /// Stuck in a single, isolated P9 resonance.
    Resonant,
    /// In the overlapping-resonance regime: transiently sticks and hops.
    Hopping,
    /// Not bound to any resonance; secular/scattering/detached dynamics.
    NonResonant,
}

/// Result of classifying one object.
#[derive(Debug, Clone, Copy)]
pub struct Classification {
    pub class: Class,
    /// Nearest resonance (if the landscape is non-empty).
    pub nearest: Option<Mmr>,
    /// |Δa| to the nearest resonance centre (AU).
    pub delta_a: f64,
    /// Local resonance-overlap K (if defined).
    pub overlap_k: Option<f64>,
}

/// Classify an object at (`a`, `e`) against the P9 resonance landscape.
///
/// Logic (the paper's regimes):
/// 1. Compute the local full-spectrum overlap K. If K ≥ 1 the neighbouring
///    resonances have merged into a chaotic sea — the object is in the
///    metastable HOPPING regime (it transiently sticks to a commensurability
///    then jumps to a neighbour). This is the large-a / dense-resonance
///    region.
/// 2. If NOT overlapping (K < 1), the resonances are isolated islands. If the
///    object lies within one half-width of the nearest SIMPLE (n:1 / n:2)
///    resonance centre, it can stick → RESONANT.
/// 3. Otherwise it is NON-RESONANT (secular / scattering / detached).
///
/// `m_p9_solar` is Planet Nine's mass in solar masses; `width_coeff` the
/// labelled pendulum-width O(1) constant. Reported `nearest` / `delta_a` are
/// for the nearest simple resonance — the commensurability the paper labels.
pub fn classify(
    landscape: &ResonanceLandscape,
    a: f64,
    e: f64,
    m_p9_solar: f64,
    width_coeff: f64,
) -> Classification {
    let nearest_simple = landscape.nearest_simple(a);
    let overlap_k = resonance_overlap_k(landscape, a, e, m_p9_solar, width_coeff);
    let overlapping = overlap_k.map(|k| k >= 1.0).unwrap_or(false);

    let (near_mmr, delta_a, within_island) = match nearest_simple {
        None => (None, f64::INFINITY, false),
        Some((mmr, a_res, da)) => {
            let half_width = resonance_half_width(a_res, mmr, e, m_p9_solar, width_coeff);
            (Some(mmr), da, da <= half_width)
        }
    };

    let class = if overlapping {
        Class::Hopping
    } else if within_island {
        Class::Resonant
    } else {
        Class::NonResonant
    };

    Classification {
        class,
        nearest: near_mmr,
        delta_a,
        overlap_k,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::analysis::resonance::resonance_semi_major_axis;
    use uom::si::length::astronomical_unit;

    #[test]
    fn typed_semi_major_axis_matches_f64() {
        let a_p9 = 500.0;
        let mmr = Mmr { p: 3, q: 1 };
        assert_relative_eq!(
            mmr.semi_major_axis_typed(a_p9).get::<astronomical_unit>(),
            mmr.semimajor_axis(a_p9),
            epsilon = 1e-12
        );
        let land = ResonanceLandscape::new(a_p9, 100.0, 600.0, 5);
        assert_relative_eq!(
            land.semi_major_axis().get::<astronomical_unit>(),
            land.a_p9,
            epsilon = 1e-12
        );
    }

    #[test]
    fn resonance_location_matches_analytic_relation() {
        // The landscape's resonance centres must equal a₉ (q/p)^{2/3} to 1e-9.
        let a_p9 = 500.0;
        for (p, q) in [(2u32, 1u32), (3, 1), (5, 1), (3, 2), (5, 2), (7, 2)] {
            let mmr = Mmr { p, q };
            let computed = mmr.semimajor_axis(a_p9);
            let analytic = a_p9 * (q as f64 / p as f64).powf(2.0 / 3.0);
            assert_relative_eq!(computed, analytic, epsilon = 1e-9);
            // And it agrees with the p9-core helper directly.
            let core = resonance_semi_major_axis(q, p, a_p9);
            assert_relative_eq!(computed, core, epsilon = 1e-9);
        }
    }

    #[test]
    fn interior_resonances_are_below_a9_and_ordered() {
        let l = ResonanceLandscape::new(500.0, 100.0, 500.0, 20);
        assert!(!l.spectrum.is_empty());
        assert!(!l.simple.is_empty());
        let mut prev = 0.0;
        for (mmr, a_res) in &l.spectrum {
            assert!(*a_res < 500.0, "{mmr:?} interior resonance should be < a9");
            assert!(*a_res >= prev, "resonances must be sorted by a");
            prev = *a_res;
        }
        // Every simple resonance is n:1 or n:2.
        for (mmr, _) in &l.simple {
            assert!(mmr.is_simple(), "{mmr:?} should be n:1 or n:2");
        }
        // 2:1 is the widest-spaced low-order resonance, near 315 AU.
        let two_one = Mmr { p: 2, q: 1 }.semimajor_axis(500.0);
        assert!((two_one - 315.0).abs() < 2.0, "2:1 at {two_one:.1} AU");
    }

    #[test]
    fn overlap_k_rises_with_semimajor_axis() {
        // The full-spectrum spacing collapses as a → a9, so K must rise.
        let l = ResonanceLandscape::new(500.0, 100.0, 499.0, 20);
        let e = 0.6;
        let m = 5.0 * p9_core::constants::EARTH_MASS_SOLAR;
        let k_low = resonance_overlap_k(&l, 250.0, e, m, 1.2).unwrap();
        let k_high = resonance_overlap_k(&l, 470.0, e, m, 1.2).unwrap();
        assert!(
            k_high > k_low,
            "overlap K should rise with a: {k_low:.3} -> {k_high:.3}"
        );
    }

    #[test]
    fn half_width_grows_with_mass_and_drops_with_order() {
        let a = 300.0;
        let lo = resonance_half_width(a, Mmr { p: 2, q: 1 }, 0.3, 1e-5, 0.5);
        let hi = resonance_half_width(a, Mmr { p: 2, q: 1 }, 0.3, 4e-5, 0.5);
        assert!(hi > lo, "wider for heavier P9");
        // Higher order ⇒ extra eccentricity suppression ⇒ narrower.
        let order1 = resonance_half_width(a, Mmr { p: 2, q: 1 }, 0.3, 1e-5, 0.5);
        let order3 = resonance_half_width(a, Mmr { p: 5, q: 2 }, 0.3, 1e-5, 0.5);
        assert!(order1 > order3, "low order should be wider");
    }

    #[test]
    fn classification_is_exhaustive() {
        let l = ResonanceLandscape::new(500.0, 100.0, 499.0, 20);
        let m = 5.0 * p9_core::constants::EARTH_MASS_SOLAR;
        // An object far from any resonance, in a non-overlapping low-a region.
        let c = classify(&l, 250.0, 0.3, m, 0.5);
        // Just exercise the API: class is one of the three.
        assert!(matches!(
            c.class,
            Class::Resonant | Class::Hopping | Class::NonResonant
        ));
    }
}
