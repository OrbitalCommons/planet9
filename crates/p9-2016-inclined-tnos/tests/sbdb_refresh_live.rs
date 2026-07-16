//! Live SBDB drift check for the frozen paper-TNO table.
//!
//! Network-gated: `cargo test -p p9-2016-inclined-tnos --features sbdb-refresh -- --ignored`
#![cfg(feature = "sbdb-refresh")]

use p9_2016_inclined_tnos::known_objects::refresh_from_sbdb;
use p9_core::data::refresh::SbdbClient;

#[test]
#[ignore = "hits the live JPL SBDB API"]
fn paper_tno_table_matches_live_sbdb() {
    let client = SbdbClient::new().expect("HTTP client");
    let diffs = refresh_from_sbdb(&client).expect("live SBDB refresh");
    for d in &diffs {
        eprintln!("paper-tno drift: {d}");
    }
    assert!(
        diffs.is_empty(),
        "{} element(s) of the paper TNO table drifted vs live SBDB (see stderr)",
        diffs.len()
    );
}
