//! Film export for `p9-2024-panstarrs`: the numbers its scene and ledger entry draw.

use p9_2024_panstarrs::combined_exclusion::compute_combined_from_population;
use serde_json::{Value, json};

pub fn export() -> Value {
    let ex = compute_combined_from_population(3000, 2024);
    json!({
        "ps1_unique": ex.ps1_unique,
        "cumulative": ex.combined,
        "remaining": 1.0 - ex.combined,
    })
}
