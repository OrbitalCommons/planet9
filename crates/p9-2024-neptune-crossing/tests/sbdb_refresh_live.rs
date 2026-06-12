//! Live SBDB checks for the frozen 17-object sample table and the
//! first-observation distance table. `#[ignore]`d and feature-gated; run:
//!
//! ```text
//! cargo test -p p9-2024-neptune-crossing --features sbdb-refresh -- --ignored
//! ```
#![cfg(feature = "sbdb-refresh")]

use p9_2024_neptune_crossing::hypothesis_test::{fetch_first_obs_distances, FIRST_OBS_DISTANCES};
use p9_2024_neptune_crossing::observed_tnos::refresh_from_sbdb;
use p9_core::data::refresh::SbdbClient;

#[test]
#[ignore = "hits the live JPL SBDB API"]
fn sample_table_matches_live_sbdb() {
    let client = SbdbClient::new().expect("HTTP client");
    let diffs = refresh_from_sbdb(&client).expect("live SBDB refresh");
    for d in &diffs {
        eprintln!("DRIFT: {d}");
    }
    assert!(
        diffs.is_empty(),
        "{} element(s) of the 17-object sample drifted out of tolerance vs live SBDB \
         (transcription error or orbit-solution update; see stderr)",
        diffs.len()
    );
}

#[test]
#[ignore = "hits the live JPL SBDB API"]
fn first_obs_distances_match_live_sbdb() {
    let client = SbdbClient::new().expect("HTTP client");
    let live = fetch_first_obs_distances(&client).expect("live SBDB fetch");
    assert_eq!(live.len(), FIRST_OBS_DISTANCES.len());

    let mut mismatches = 0;
    for (name, date, r) in &live {
        eprintln!("    (\"{name}\", \"{date}\", {r:.3}),");
        let Some((_, snap_date, snap_r)) = FIRST_OBS_DISTANCES
            .iter()
            .find(|(snap_name, _, _)| snap_name == name)
        else {
            eprintln!("UNEXPECTED: {name} missing from FIRST_OBS_DISTANCES");
            mismatches += 1;
            continue;
        };
        // The distance moves with the orbit solution; 1% covers normal a-e
        // refits without hiding a wrong date or a wrong row.
        if snap_date != date || (snap_r - r).abs() > 0.01 * r.abs().max(1.0) {
            eprintln!(
                "UNEXPECTED: {name}: const table says ({snap_date}, {snap_r}), live says \
                 ({date}, {r:.3})"
            );
            mismatches += 1;
        }
    }
    assert_eq!(
        mismatches, 0,
        "FIRST_OBS_DISTANCES no longer matches the live SBDB recomputation (see stderr; \
         the rows printed above are paste-ready replacements)"
    );
}
