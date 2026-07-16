//! Placing a representative distant TNO within the Planet Nine resonance
//! structure and deciding whether it lives in the isolated (locked) or the
//! overlapping (hopping) regime.
//!
//! Becker+ (2017) single out two behavioural classes: objects like Sedna and
//! 2012 VP113 that "do not migrate significantly" and "tend to stay in the same
//! resonance", versus objects like 2007 TG422 and 2013 RF98 that "may both
//! migrate and transition between different resonances". An eccentric TNO
//! whose perihelion penetrates the first-order resonance-overlap zone near P9
//! (inner edge `a_hop ≈ 590 AU`, computed in `chain`) samples the chaotic web
//! during its perihelion passage and can hop; an object whose orbit stays
//! clear of that zone stays locked. We use aphelion `Q = a(1+e)` as the
//! penetration depth: a high-e TNO with a few hundred AU reaches well past
//! a_hop at aphelion.

use crate::chain::{
    first_order_chain, hopping_threshold_au, mass_ratio, n_over_1_chain, published, P9Resonance,
};
use p9_core::data::etno::{Etno, BROWN_2017_SAMPLE};
use p9_core::units::au;

/// Classification of a TNO's dynamical state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum HoppingState {
    /// The orbit stays clear of the overlap zone: resonances are isolated and
    /// the object tends to stay in one resonance.
    Locked,
    /// The orbit penetrates the overlap zone: resonances overlap and the
    /// object hops between adjacent resonances.
    Hopping,
}

/// Inner edge (AU) of the first-order resonance-overlap zone for the nominal
/// Planet Nine and a representative high-e TNO — the computed hopping
/// threshold `a_hop`.
pub fn hopping_threshold_au_nominal() -> f64 {
    let chain = first_order_chain(2, 200);
    let mu = mass_ratio(published::M9_EARTH);
    hopping_threshold_au(&chain, published::A9_AU, mu, published::TNO_E)
        .expect("nominal P9 first-order chain crosses K=1")
}

/// Classify a TNO by whether its orbit (aphelion `Q = a(1+e)`) reaches into
/// the overlap zone.
pub fn classify(a_au: f64, e: f64) -> HoppingState {
    let aphelion = a_au * (1.0 + e);
    if aphelion >= hopping_threshold_au_nominal() {
        HoppingState::Hopping
    } else {
        HoppingState::Locked
    }
}

/// The closest interior n:1 resonance (by semimajor axis) to a given a — the
/// resonance a locked object is most likely sitting in.
pub fn nearest_n1_resonance(a_au: f64) -> P9Resonance {
    let chain = n_over_1_chain(2, 200);
    *chain
        .iter()
        .min_by(|x, y| {
            (x.semi_major_axis_typed(published::A9_AU) - au(a_au))
                .abs()
                .partial_cmp(&(y.semi_major_axis_typed(published::A9_AU) - au(a_au)).abs())
                .unwrap()
        })
        .unwrap()
}

/// Classify every ETNO in the Brown (2017) sample.
pub fn classify_sample() -> Vec<(&'static str, f64, f64, HoppingState)> {
    BROWN_2017_SAMPLE
        .iter()
        .map(|o: &Etno| (o.name, o.a, o.e, classify(o.a, o.e)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obj(name: &str) -> &'static Etno {
        BROWN_2017_SAMPLE.iter().find(|o| o.name == name).unwrap()
    }

    #[test]
    fn threshold_is_inside_p9_in_the_distant_band() {
        let t = hopping_threshold_au_nominal();
        assert!(t > 400.0 && t < published::A9_AU, "a_hop = {t}");
    }

    #[test]
    fn classification_splits_on_aphelion_penetration() {
        let t = hopping_threshold_au_nominal();
        // A circular object well inside the zone: Locked.
        assert_eq!(classify(0.3 * t, 0.0), HoppingState::Locked);
        // The same a but high e reaches the zone at aphelion: Hopping.
        assert_eq!(classify(0.7 * t, 0.9), HoppingState::Hopping);
    }

    #[test]
    fn eccentric_distant_tnos_penetrate_overlap_zone() {
        // Becker+'s MIGRATING objects are high-e and distant; their aphelia
        // reach into the overlap zone, so they are classified as hoppers.
        for name in ["2007 TG422", "2013 RF98"] {
            let o = obj(name);
            assert_eq!(
                classify(o.a, o.e),
                HoppingState::Hopping,
                "{name} (a={}, e={}, Q={})",
                o.a,
                o.e,
                o.a * (1.0 + o.e)
            );
        }
    }

    #[test]
    fn sedna_is_a_known_false_positive_of_the_aphelion_proxy() {
        // Becker et al. report Sedna (with 2012 VP113) as an object that
        // does NOT migrate — it "tends to stay in the same resonance" — yet
        // its aphelion (Q ≈ 936 AU ≥ a_hop ≈ 590) reaches the overlap zone,
        // so the Q-based proxy classifies it Hopping. The proxy is a
        // NECESSARY-not-sufficient criterion: reaching the chaotic web at
        // aphelion does not force migration when the secular phase keeps
        // perihelion passages resonantly protected (the paper's N-body
        // captures this; a static aphelion cut cannot). Documented as a
        // known false positive rather than asserted as reproducing the
        // paper's class.
        let o = obj("Sedna");
        assert_eq!(
            classify(o.a, o.e),
            HoppingState::Hopping,
            "the proxy's false positive is itself pinned (Q = {:.0})",
            o.a * (1.0 + o.e)
        );
    }

    #[test]
    fn nearest_n1_resonance_is_within_local_spacing() {
        // The n:1 resonances are sparse out near P9 (2:1 at 441 AU, then 1:1
        // at a9 = 700), so the nearest can be many tens of AU away — the very
        // reason these distant objects are not cleanly locked in one n:1
        // resonance. We just check the lookup returns a sensible neighbour:
        // within the local 2:1→1:1 gap.
        let o = obj("2007 TG422"); // a = 501 AU
        let res = nearest_n1_resonance(o.a);
        let local_gap =
            au(published::A9_AU) - P9Resonance::new(2, 1).semi_major_axis_typed(published::A9_AU);
        let da = (res.semi_major_axis_typed(published::A9_AU) - au(o.a)).abs();
        assert!(
            da < local_gap,
            "nearest n:1 resonance {da:?} away (gap {local_gap:?})"
        );
        // And the nearest n:1 to 2007 TG422 is the 2:1.
        assert_eq!(res, P9Resonance::new(2, 1));
    }

    #[test]
    fn whole_sample_classifies_without_panic() {
        let rows = classify_sample();
        assert_eq!(rows.len(), BROWN_2017_SAMPLE.len());
        assert!(rows.iter().any(|(_, _, _, s)| *s == HoppingState::Hopping));
    }
}
