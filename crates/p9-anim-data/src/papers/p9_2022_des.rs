//! Film export for `p9-2022-des`: the numbers its scene and ledger entry draw.

use p9_2024_panstarrs::combined_exclusion::compute_combined_from_population;
use serde_json::{Value, json};

pub fn export() -> Value {
    let ex = compute_combined_from_population(3000, 2024);
    json!({
        "des_unique": ex.des_unique,
        "cumulative": ex.ztf_frac + ex.des_unique,
    })
}
