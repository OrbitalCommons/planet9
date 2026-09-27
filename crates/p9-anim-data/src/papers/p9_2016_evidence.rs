//! Film export for `p9-2016-evidence`: the numbers its scene and ledger entry draw.

use p9_2016_evidence::kbo_elements::joint_clustering_significance;
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::data::stable_kbos::stable_kbos;
use p9_core::types::P9Params;
use serde_json::{Value, json};

pub fn export() -> Value {
    let joint = joint_clustering_significance(2_000_000, 2016);
    let p9 = P9Params::nominal_2016();
    json!({
        "p_joint": joint.p_joint,
        "sigma": p_value_to_sigma(joint.p_joint),
        "n_sample": stable_kbos().len(),
        "mass": p9.mass_earth,
        "a": p9.a,
        "e": p9.e,
        "i": p9.i.to_degrees(),
    })
}
