//! Live SBDB drift check for the frozen extended distant-TNO table.
//!
//! Network-gated: `cargo test -p p9-2021-perihelion-gap --features sbdb-refresh -- --ignored`
#![cfg(feature = "sbdb-refresh")]

use p9_2021_perihelion_gap::sample::refresh_from_sbdb;
use p9_core::data::refresh::SbdbClient;

#[test]
#[ignore = "hits the live JPL SBDB API"]
fn extended_table_matches_live_sbdb() {
    let client = SbdbClient::new().expect("HTTP client");
    let diffs = refresh_from_sbdb(&client).expect("live SBDB refresh");
    for d in &diffs {
        eprintln!("extended-table drift: {d}");
    }
    assert!(
        diffs.is_empty(),
        "{} element(s) of EXTENDED_DISTANT_TNOS drifted vs live SBDB (see stderr)",
        diffs.len()
    );
}
