//! Preface and finale film exports: one module per scene file, each returning the
//! values its scenes read from `anim.json -> preface -> <module>`.

use std::collections::BTreeMap;

use serde_json::Value;

pub mod a_foundations;
pub mod b_elements;
pub mod c_dynamics;
pub mod d_mechanism;
pub mod e_detection;
pub mod f_indirect_map;
pub mod finale;

pub fn all() -> BTreeMap<&'static str, Value> {
    BTreeMap::from([
        ("a_foundations", a_foundations::export()),
        ("b_elements", b_elements::export()),
        ("c_dynamics", c_dynamics::export()),
        ("d_mechanism", d_mechanism::export()),
        ("e_detection", e_detection::export()),
        ("f_indirect_map", f_indirect_map::export()),
        ("finale", finale::export()),
    ])
}
