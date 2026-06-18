//! Exterior-perturber (Planet Nine) mean-motion resonance chains near a
//! distant TNO: the n:1 and n:2 resonance semimajor axes, the libration width
//! of each, and the Chirikov resonance-overlap parameter between neighbours.
//!
//! ## Physics (Becker, Adams et al. 2017, arXiv:1706.06609)
//!
//! Becker+ integrate distant TNOs in the presence of the nominal Planet Nine
//! and find the objects are *not* locked in a single P9 mean-motion resonance
//! for their lifetimes. They HOP between adjacent resonances (metastable
//! resonance sticking) while keeping approximate apsidal anti-alignment,
//! because neighbouring P9 resonances OVERLAP. Remaining in one resonance is
//! not a requirement for stability.
//!
//! The analytic backbone reproduced here:
//!
//! 1. **Resonance locations.** A particle interior to P9 in the p:q resonance
//!    (p > q: the particle completes q orbits while P9 completes p) sits at
//!    `a_res = a9 (q/p)^{2/3} < a9`. We delegate to the single workspace Kepler
//!    relation `p9_core::analysis::resonance::resonance_semi_major_axis`, whose
//!    `(p, q)` arguments map to `(q, p)` here (it returns `a_planet (p/q)^{2/3}`
//!    for the *exterior* case, so the interior location is `…(q, p, a9)`).
//!
//! 2. **Libration width.** Each resonance has a finite half-width in semimajor
//!    axis from the pendulum model (Murray & Dermott 1999, Ch. 8):
//!    `δa = a_res · sqrt( (16/3) · μ · S )`, μ = m9 / M_sun, with a
//!    dimensionless strength `S = α · α/(1−α) · e` (`resonance_strength`).
//!    The factor `α/(1−α)` (α = a_res/a9 < 1) is the leading growth of the
//!    Laplace-coefficient amplitude as the resonance crowds toward the
//!    perturber; the linear `e` is the eccentric-encounter scaling appropriate
//!    to the high-e distant TNOs (e ≈ 0.9) for which the d'Alembert `e^{order}`
//!    truncation does not converge. No free constant is tuned to a target a.
//!
//! 3. **Overlap.** Two neighbouring resonances overlap when the sum of their
//!    half-widths exceeds their separation — the Chirikov criterion
//!    `K = (δa_lo + δa_hi)/|a_hi − a_lo|`. K < 1: isolated, an object can stay
//!    put. K ≥ 1: the resonances merge into a chaotic hopping band.
//!
//! The widely-spaced n:1 / n:2 resonances out in the TNO region are isolated
//! for a nominal-mass P9; it is the **first-order p:(p−1) chain crowding toward
//! P9** whose neighbour overlap parameter rises monotonically and crosses
//! K = 1 at a computed `a_hop ≈ 590 AU` (inside a9 = 700 AU) — the inner edge
//! of the Wisdom (1980) resonance-overlap zone. The eccentric TNOs whose
//! perihelia penetrate this zone are the ones Becker+ find hopping; see the
//! `tno` module and REPRODUCTION_NOTES for the honest scale caveat.

use p9_core::analysis::resonance::resonance_semi_major_axis;
use p9_core::units::{au, Length};

/// Published / reference constants for Becker, Adams et al. (2017).
pub mod published {
    /// Nominal Planet Nine semimajor axis (AU); the Batygin & Brown orbit.
    pub const A9_AU: f64 = 700.0;
    /// Nominal Planet Nine mass in Earth masses.
    pub const M9_EARTH: f64 = 10.0;
    /// Nominal Planet Nine eccentricity.
    pub const E9: f64 = 0.6;
    /// Representative distant-TNO eccentricity used for resonance widths
    /// (Becker+'s hoppers — 2007 TG422, 2013 RF98 — have e ≈ 0.9).
    pub const TNO_E: f64 = 0.9;
    /// Becker+ headline: distant TNOs hop between adjacent resonances rather
    /// than staying locked in one.
    pub const HOPS_BETWEEN_RESONANCES: bool = true;
}

/// Earth mass in solar masses (Earth alone; matches p9-core's constant).
const EARTH_MASS_SOLAR: f64 = 3.003_489e-6;

/// A single interior p:q mean-motion resonance of Planet Nine.
///
/// Convention (the workspace single convention): an interior particle in the
/// p:q resonance completes q orbits while P9 completes p (p > q), sitting at
/// `a_res = a9 (q/p)^{2/3} < a9`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct P9Resonance {
    /// Planet Nine's integer (orbits of P9 per `q` orbits of the particle).
    pub p: u32,
    /// Particle's integer.
    pub q: u32,
}

impl P9Resonance {
    pub fn new(p: u32, q: u32) -> Self {
        Self { p, q }
    }

    /// Order of the resonance |p − q|.
    pub fn order(&self) -> u32 {
        self.p.abs_diff(self.q)
    }

    /// Interior resonance semimajor axis as a typed [`Length`] for a perturber
    /// at `a9`: `a_res = a9 (q/p)^{2/3}`. Delegates to the single workspace
    /// Kepler relation (with arguments swapped to select the interior branch).
    pub fn semi_major_axis_typed(&self, a9_au: f64) -> Length {
        au(resonance_semi_major_axis(self.q, self.p, a9_au))
    }

    /// Dimensionless ratio α = a_res / a9 (< 1 for interior resonances).
    pub fn alpha(&self) -> f64 {
        (self.q as f64 / self.p as f64).powf(2.0 / 3.0)
    }

    /// Pendulum libration half-width in semimajor axis (AU) for a particle of
    /// eccentricity `e` perturbed by a planet of mass ratio `mu` at `a9`:
    /// `δa = a_res · sqrt( (16/3) · mu · S )`, S from `resonance_strength`.
    pub fn libration_half_width_au(&self, a9_au: f64, mu: f64, e: f64) -> f64 {
        let a_res = (self.semi_major_axis_typed(a9_au) / au(1.0)).value;
        let strength = resonance_strength(self.alpha(), e);
        a_res * ((16.0 / 3.0) * mu * strength).sqrt()
    }

    /// Pendulum libration half-width as a typed [`Length`] (see
    /// [`Self::libration_half_width_au`]).
    pub fn libration_half_width(&self, a9_au: f64, mu: f64, e: f64) -> Length {
        au(self.libration_half_width_au(a9_au, mu, e))
    }
}

/// Dimensionless resonance strength `S(α, e) = α · α/(1−α) · e` (see module
/// docs). Grows toward the perturber (α → 1) and with eccentricity. No tuned
/// constant.
pub fn resonance_strength(alpha: f64, e: f64) -> f64 {
    alpha * (alpha / (1.0 - alpha)) * e
}

/// Convert a Planet Nine mass in Earth masses to a solar mass ratio μ.
pub fn mass_ratio(m9_earth: f64) -> f64 {
    m9_earth * EARTH_MASS_SOLAR
}

/// The interior n:1 resonance chain `p:1` for `p ∈ [p_min, p_max]`, sorted by
/// increasing semimajor axis (which is decreasing p — larger p sits closer in).
pub fn n_over_1_chain(p_min: u32, p_max: u32) -> Vec<P9Resonance> {
    sorted_by_a((p_min..=p_max).map(|p| P9Resonance::new(p, 1)))
}

/// The interior n:2 resonance chain `p:2` (odd p in lowest terms), sorted by
/// increasing semimajor axis.
pub fn n_over_2_chain(p_min: u32, p_max: u32) -> Vec<P9Resonance> {
    sorted_by_a(
        (p_min..=p_max)
            .filter(|p| p % 2 == 1)
            .map(|p| P9Resonance::new(p, 2)),
    )
}

/// The first-order p:(p−1) resonance chain — the resonances that crowd toward
/// the perturber and whose neighbour overlap parameter crosses K = 1 (the
/// hopping-band inner edge). Sorted by increasing semimajor axis.
pub fn first_order_chain(p_min: u32, p_max: u32) -> Vec<P9Resonance> {
    sorted_by_a((p_min..=p_max).map(|p| P9Resonance::new(p, p - 1)))
}

fn sorted_by_a(iter: impl Iterator<Item = P9Resonance>) -> Vec<P9Resonance> {
    let mut v: Vec<P9Resonance> = iter.collect();
    v.sort_by(|x, y| x.alpha().partial_cmp(&y.alpha()).unwrap());
    v
}

/// Chirikov overlap parameter between two neighbouring resonances:
/// `K = (δa_lo + δa_hi) / |a_hi − a_lo|`. K ≥ 1 ⇒ overlap (hopping).
pub fn neighbour_overlap(lo: P9Resonance, hi: P9Resonance, a9_au: f64, mu: f64, e: f64) -> f64 {
    let a_lo = (lo.semi_major_axis_typed(a9_au) / au(1.0)).value;
    let a_hi = (hi.semi_major_axis_typed(a9_au) / au(1.0)).value;
    let sep = (a_hi - a_lo).abs();
    let w = lo.libration_half_width_au(a9_au, mu, e) + hi.libration_half_width_au(a9_au, mu, e);
    w / sep
}

/// One link in the overlap analysis of an ordered resonance chain.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct OverlapLink {
    pub lo: P9Resonance,
    pub hi: P9Resonance,
    /// Mean semimajor axis of the pair (AU).
    pub a_mid_au: f64,
    /// Separation between the two resonance centres (AU).
    pub separation_au: f64,
    /// Sum of the two libration half-widths (AU).
    pub width_sum_au: f64,
    /// Chirikov overlap parameter K.
    pub k: f64,
}

impl OverlapLink {
    pub fn overlaps(&self) -> bool {
        self.k >= 1.0
    }

    /// Mean semimajor axis of the pair as a typed [`Length`].
    pub fn a_mid(&self) -> Length {
        au(self.a_mid_au)
    }

    /// Separation between the two resonance centres as a typed [`Length`].
    pub fn separation(&self) -> Length {
        au(self.separation_au)
    }

    /// Sum of the two libration half-widths as a typed [`Length`].
    pub fn width_sum(&self) -> Length {
        au(self.width_sum_au)
    }
}

/// Build the overlap links for an ordered (increasing-a) chain.
pub fn overlap_profile(chain: &[P9Resonance], a9_au: f64, mu: f64, e: f64) -> Vec<OverlapLink> {
    chain
        .windows(2)
        .map(|w| {
            let (lo, hi) = (w[0], w[1]);
            let a_lo = (lo.semi_major_axis_typed(a9_au) / au(1.0)).value;
            let a_hi = (hi.semi_major_axis_typed(a9_au) / au(1.0)).value;
            OverlapLink {
                lo,
                hi,
                a_mid_au: 0.5 * (a_lo + a_hi),
                separation_au: (a_hi - a_lo).abs(),
                width_sum_au: lo.libration_half_width_au(a9_au, mu, e)
                    + hi.libration_half_width_au(a9_au, mu, e),
                k: neighbour_overlap(lo, hi, a9_au, mu, e),
            }
        })
        .collect()
}

/// The semimajor axis (AU) at which a chain crosses into the hopping regime:
/// the mid-pair a of the first overlap link with K ≥ 1 (K rises monotonically
/// toward P9). `None` if no link reaches K = 1 in the supplied chain.
pub fn hopping_threshold_au(chain: &[P9Resonance], a9_au: f64, mu: f64, e: f64) -> Option<f64> {
    overlap_profile(chain, a9_au, mu, e)
        .into_iter()
        .find(|link| link.overlaps())
        .map(|link| link.a_mid_au)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::analysis::resonance::chirikov_overlap_parameter;
    use uom::si::length::astronomical_unit;

    const A9: f64 = published::A9_AU;

    #[test]
    fn typed_accessors_match_f64_sources() {
        let res = P9Resonance::new(5, 4);
        let mu = mass_ratio(published::M9_EARTH);
        assert_relative_eq!(
            res.semi_major_axis_typed(A9).get::<astronomical_unit>(),
            resonance_semi_major_axis(res.q, res.p, A9),
            epsilon = 1e-9
        );
        assert_relative_eq!(
            res.libration_half_width(A9, mu, 0.9)
                .get::<astronomical_unit>(),
            res.libration_half_width_au(A9, mu, 0.9),
            epsilon = 1e-9
        );

        let chain = first_order_chain(2, 20);
        let link = overlap_profile(&chain, A9, mu, published::TNO_E)[0];
        assert_relative_eq!(
            link.a_mid().get::<astronomical_unit>(),
            link.a_mid_au,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            link.separation().get::<astronomical_unit>(),
            link.separation_au,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            link.width_sum().get::<astronomical_unit>(),
            link.width_sum_au,
            epsilon = 1e-9
        );
    }

    #[test]
    fn resonance_locations_match_kepler_relation() {
        // a_res = a9 (q/p)^{2/3} to 1e-9 — the load-bearing analytic relation.
        for &(p, q) in &[(2u32, 1u32), (3, 1), (5, 1), (3, 2), (5, 2), (4, 3)] {
            let res = P9Resonance::new(p, q);
            let expected = A9 * (q as f64 / p as f64).powf(2.0 / 3.0);
            assert_relative_eq!(
                res.semi_major_axis_typed(A9).get::<astronomical_unit>(),
                expected,
                epsilon = 1e-9
            );
            // All interior resonances sit inside P9.
            assert!(res.semi_major_axis_typed(A9) < au(A9));
        }
    }

    #[test]
    fn n_over_1_chain_is_ordered_and_crowds_toward_p9() {
        let chain = n_over_1_chain(2, 30);
        for w in chain.windows(2) {
            assert!(w[1].semi_major_axis_typed(A9) > w[0].semi_major_axis_typed(A9));
        }
        // Spacing shrinks toward smaller a (higher p crowds together): the
        // first (small-a) gap is smaller than the last (large-a, near 2:1).
        let seps: Vec<Length> = chain
            .windows(2)
            .map(|w| w[1].semi_major_axis_typed(A9) - w[0].semi_major_axis_typed(A9))
            .collect();
        assert!(seps.first().unwrap() < seps.last().unwrap());
    }

    #[test]
    fn n_over_2_interleaves_n_over_1() {
        // A (2p+1):2 resonance sits between the (p+1):1 and p:1 resonances.
        let p = 5u32;
        let inner = P9Resonance::new(p + 1, 1).semi_major_axis_typed(A9); // smaller a
        let outer = P9Resonance::new(p, 1).semi_major_axis_typed(A9); // larger a
        let half = P9Resonance::new(2 * p + 1, 2).semi_major_axis_typed(A9);
        assert!(
            inner < half && half < outer,
            "{inner:?} < {half:?} < {outer:?}"
        );
        // n:2 chain is odd-p and ordered.
        let chain = n_over_2_chain(3, 21);
        assert!(chain.iter().all(|r| r.p % 2 == 1));
        for w in chain.windows(2) {
            assert!(w[1].semi_major_axis_typed(A9) > w[0].semi_major_axis_typed(A9));
        }
    }

    #[test]
    fn overlap_rises_monotonically_and_crosses_hopping_threshold() {
        // The first-order chain crowds toward P9; its overlap parameter rises
        // monotonically and crosses K = 1 — Becker+'s hopping band edge.
        let chain = first_order_chain(2, 60);
        let mu = mass_ratio(published::M9_EARTH);
        let profile = overlap_profile(&chain, A9, mu, published::TNO_E);

        for w in profile.windows(2) {
            assert!(w[1].k > w[0].k, "K not rising: {} -> {}", w[0].k, w[1].k);
        }
        // Isolated at the small-a end, overlapping at the large-a (near-P9) end.
        assert!(
            profile.first().unwrap().k < 1.0,
            "innermost should be isolated"
        );
        assert!(profile.last().unwrap().k > 1.0, "near-P9 should overlap");

        let a_cross =
            hopping_threshold_au(&chain, A9, mu, published::TNO_E).expect("chain should cross K=1");
        // Inner edge of the overlap zone sits inside P9, in the distant band.
        assert!(
            a_cross > 400.0 && a_cross < A9,
            "hopping threshold a = {a_cross}"
        );
    }

    #[test]
    fn below_threshold_isolated_above_overlapping() {
        let chain = first_order_chain(2, 60);
        let mu = mass_ratio(published::M9_EARTH);
        let a_cross = hopping_threshold_au(&chain, A9, mu, published::TNO_E).unwrap();
        let profile = overlap_profile(&chain, A9, mu, published::TNO_E);
        for link in &profile {
            if link.a_mid_au < a_cross {
                assert!(
                    link.k < 1.0,
                    "K<1 expected below threshold at {}",
                    link.a_mid_au
                );
            }
        }
        // Monotone K ⇒ every link at or above the threshold overlaps.
        assert!(profile
            .iter()
            .filter(|l| l.a_mid_au >= a_cross)
            .all(|l| l.k >= 1.0));
    }

    #[test]
    fn n_over_1_resonances_isolated_for_nominal_p9() {
        // The widely-spaced n:1 resonances out in the TNO region do NOT overlap
        // for a nominal-mass P9 — they are the resonances a locked object sits
        // in. The honest scale finding (REPRODUCTION_NOTES): only the crowded
        // first-order chain near P9 overlaps.
        let chain = n_over_1_chain(2, 40);
        let mu = mass_ratio(published::M9_EARTH);
        let profile = overlap_profile(&chain, A9, mu, published::TNO_E);
        assert!(
            profile.iter().all(|l| l.k < 1.0),
            "nominal-P9 n:1 chain should be isolated everywhere"
        );
    }

    #[test]
    fn libration_width_and_overlap_scale_as_sqrt_mu() {
        // δa ∝ √μ, separation is μ-free ⇒ K ∝ √μ exactly.
        let res = P9Resonance::new(5, 4);
        let w1 = res.libration_half_width_au(A9, mass_ratio(10.0), 0.9);
        let w4 = res.libration_half_width_au(A9, mass_ratio(40.0), 0.9);
        assert_relative_eq!(w4 / w1, 2.0, epsilon = 1e-12);

        let (lo, hi) = (P9Resonance::new(5, 4), P9Resonance::new(6, 5));
        let k1 = neighbour_overlap(lo, hi, A9, mass_ratio(10.0), 0.9);
        let k4 = neighbour_overlap(lo, hi, A9, mass_ratio(40.0), 0.9);
        assert_relative_eq!(k4 / k1, 2.0, epsilon = 1e-12);
    }

    #[test]
    fn heavier_planet_and_higher_e_do_not_raise_threshold() {
        // Wider resonances (heavier P9, higher e) overlap at the same or
        // smaller a (the crossing is discretised by the chain).
        let chain = first_order_chain(2, 80);
        let a_light = hopping_threshold_au(&chain, A9, mass_ratio(5.0), 0.9).unwrap();
        let a_heavy = hopping_threshold_au(&chain, A9, mass_ratio(20.0), 0.9).unwrap();
        assert!(a_heavy <= a_light, "heavier: {a_heavy} vs {a_light}");

        let a_lo_e = hopping_threshold_au(&chain, A9, mass_ratio(10.0), 0.5).unwrap();
        let a_hi_e = hopping_threshold_au(&chain, A9, mass_ratio(10.0), 0.9).unwrap();
        assert!(a_hi_e <= a_lo_e, "higher e: {a_hi_e} vs {a_lo_e}");
    }

    #[test]
    fn overlap_a_exponent_brackets_validated_chirikov_slope() {
        // Cross-check the a-dependence of our overlap parameter against the
        // validated Neptune-form Chirikov K ∝ (a/a_pert)^{5/4} in
        // p9_core::analysis::resonance::chirikov_overlap_parameter. Our K rises
        // even faster toward P9 (the α/(1−α) crowding diverges at a9), so its
        // local log-slope should equal or exceed 5/4 over a mid-range window.
        let chain = first_order_chain(2, 60);
        let mu = mass_ratio(published::M9_EARTH);
        let profile = overlap_profile(&chain, A9, mu, published::TNO_E);
        let l1 = profile.iter().find(|l| l.a_mid_au > 480.0).unwrap();
        let l2 = profile.iter().find(|l| l.a_mid_au > 560.0).unwrap();
        let slope = (l2.k.ln() - l1.k.ln()) / (l2.a_mid_au.ln() - l1.a_mid_au.ln());
        assert!(
            slope >= 1.25,
            "overlap a-exponent {slope} should be >= Chirikov 5/4"
        );

        // Anchor: the validated core function rises with a, as ours does.
        assert!(chirikov_overlap_parameter(500.0, 40.0) > chirikov_overlap_parameter(300.0, 40.0));
    }
}
