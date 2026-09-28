//! Film export for `p9-2017-resonance-hopping`: the numbers its scene and ledger entry draw.
//!
//! Becker et al. (2017) find that distant objects hop between neighbouring
//! Planet Nine resonances instead of staying locked in one. The reproduction
//! crate gives the analytic reason: the libration width of each resonance, the
//! Chirikov overlap parameter between neighbours, the semimajor axis where the
//! first-order chain starts to overlap, and which observed objects reach it.

use p9_2017_resonance_hopping::chain::resonance_strength;
use p9_2017_resonance_hopping::{
    P9Resonance, classify_sample, first_order_chain, hopping_threshold_au_nominal, mass_ratio,
    n_over_1_chain, n_over_2_chain, overlap_profile, published,
};
use p9_core::analysis::resonance::overlap_zone_width_fraction;
use p9_core::units::au;
use serde_json::{Value, json};

/// Longest first-order chain used for the overlap profile.
const P_MAX: u32 = 200;

fn resonance_row(res: &P9Resonance, mu: f64) -> Value {
    let a9 = published::A9_AU;
    json!({
        "p": res.p,
        "q": res.q,
        "a_au": (res.semi_major_axis_typed(a9) / au(1.0)).value,
        "half_width_au": res.libration_half_width_au(a9, mu, published::TNO_E),
        "strength": resonance_strength(res.alpha(), published::TNO_E),
    })
}

pub fn export() -> Value {
    let a9 = published::A9_AU;
    let e = published::TNO_E;
    let mu = mass_ratio(published::M9_EARTH);

    // The widely spaced n:1 and n:2 resonances out in the ETNO region.
    let simple: Vec<Value> = n_over_1_chain(2, 6)
        .iter()
        .chain(n_over_2_chain(5, 11).iter())
        .map(|r| resonance_row(r, mu))
        .collect();
    let simple_k_max = overlap_profile(&n_over_1_chain(2, 40), a9, mu, e)
        .iter()
        .map(|l| l.k)
        .fold(0.0_f64, f64::max);

    // The first-order p:(p-1) chain crowding toward Planet Nine.
    let chain = first_order_chain(2, P_MAX);
    let first_order: Vec<Value> = first_order_chain(2, 40)
        .iter()
        .map(|r| resonance_row(r, mu))
        .collect();
    let profile = overlap_profile(&chain, a9, mu, e);
    let a_mid: Vec<f64> = profile.iter().map(|l| l.a_mid_au).collect();

    let a_hop = hopping_threshold_au_nominal();
    let a_zone = a9 * (1.0 - overlap_zone_width_fraction(mu, e));

    let sample: Vec<Value> = classify_sample()
        .iter()
        .map(|(name, a, ecc, state)| {
            json!({
                "name": name,
                "a_au": a,
                "e": ecc,
                "q_au": a * (1.0 - ecc),
                "big_q_au": a * (1.0 + ecc),
                "state": format!("{state:?}"),
            })
        })
        .collect();
    let n_hopping = sample.iter().filter(|s| s["state"] == "Hopping").count();

    json!({
        "a9_au": a9,
        "m9_earth": published::M9_EARTH,
        "e9": published::E9,
        "tno_e": e,
        "simple": simple,
        "simple_k_max": simple_k_max,
        "first_order": first_order,
        "profile": {"a_mid_au": a_mid, "k": profile.iter().map(|l| l.k).collect::<Vec<f64>>()},
        "a_hop": a_hop,
        "a_zone_sourced": a_zone,
        "sample": sample,
        "n_sample": sample.len(),
        "n_hopping": n_hopping,
    })
}
