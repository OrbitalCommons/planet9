//! Live SBDB drift checks for the frozen ETNO table.
//!
//! These hit the JPL SBDB API, so they are `#[ignore]`d and additionally
//! gated behind the `sbdb-refresh` feature; default `cargo test` runs
//! nothing from this file. Run with:
//!
//! ```text
//! cargo test -p p9-core --features sbdb-refresh -- --ignored
//! ```
#![cfg(feature = "sbdb-refresh")]

use p9_core::data::refresh::{fetch_live_elements, refresh_etno_from_sbdb, Element};
use starfield::sbdb::SbdbClient;

/// Orbit-solution drift documented at the 2026-06-12 live check (see the
/// `data::etno` module docs): the table pins the Brown (2017) paper-epoch
/// solutions, and JPL's a–e fit degeneracy has since moved a/e on 8 of the
/// 10 objects (plus ω of the short-arc 2013 RF98) while q, i, Ω and H stayed
/// put. These are *expected* diffs, not table bugs; anything outside this
/// list is a new drift or a transcription error and fails the test.
const KNOWN_DRIFT: &[(&str, Element)] = &[
    ("Sedna", Element::SemiMajorAxis),
    ("Sedna", Element::Eccentricity),
    ("2012 VP113", Element::SemiMajorAxis),
    ("2013 RF98", Element::SemiMajorAxis),
    ("2013 RF98", Element::Eccentricity),
    ("2013 RF98", Element::ArgPerihelion),
    ("2004 VN112", Element::SemiMajorAxis),
    ("2004 VN112", Element::Eccentricity),
    ("2010 GB174", Element::SemiMajorAxis),
    ("2010 GB174", Element::ArgPerihelion),
    ("2007 TG422", Element::SemiMajorAxis),
    ("2013 FT28", Element::SemiMajorAxis),
    ("2013 FT28", Element::Eccentricity),
    ("2014 SR349", Element::SemiMajorAxis),
    ("2014 SR349", Element::ArgPerihelion),
    ("2015 RX245", Element::SemiMajorAxis),
];

#[test]
#[ignore = "hits the live JPL SBDB API"]
fn etno_table_matches_live_sbdb_modulo_documented_drift() {
    let client = SbdbClient::new().expect("HTTP client");
    let diffs = refresh_etno_from_sbdb(&client).expect("live SBDB refresh");
    let mut unexpected = 0;
    for d in &diffs {
        let known = KNOWN_DRIFT
            .iter()
            .any(|(obj, el)| *obj == d.object && *el == d.element);
        eprintln!("{}: {d}", if known { "known drift" } else { "UNEXPECTED" });
        if !known {
            unexpected += 1;
        }
    }
    assert_eq!(
        unexpected, 0,
        "{unexpected} element(s) of BROWN_2017_SAMPLE drifted out of tolerance vs live \
         SBDB beyond the drift documented 2026-06-12 (transcription error or new \
         orbit-solution update; see stderr)"
    );
}

#[test]
#[ignore = "hits the live JPL SBDB API"]
fn sbdb_lookup_returns_full_element_set_for_sedna() {
    let client = SbdbClient::new().expect("HTTP client");
    let live = fetch_live_elements(&client, "Sedna").expect("Sedna lookup");
    assert!(live.a.is_some());
    assert!(live.e.is_some());
    assert!(live.i_deg.is_some());
    assert!(live.omega_deg.is_some());
    assert!(live.omega_big_deg.is_some());
    assert!(live.h_mag.is_some());
    assert!(live.epoch_jd.is_some());
    assert!(live.mean_anomaly_deg.is_some());
    assert!(live.mean_motion_deg_per_day.is_some());
    assert!(live.first_obs.is_some());
}

/// Solution drift documented at the 2026-07-15 live check: `stable_kbos`
/// deliberately pins the paper-epoch (~2016) solutions (see the module
/// docs), and JPL has since republished a/e along the fit degeneracy for
/// the long-period members (live SBDB also serves only ~3 significant
/// figures for these), plus a sub-degree ω update for the short-arc
/// 2010 GB174. Anything outside this list is a transcription error or a
/// new solution worth a human look.
const STABLE_KBO_KNOWN_DRIFT: &[(&str, Element)] = &[
    ("Sedna", Element::SemiMajorAxis),
    ("Sedna", Element::Eccentricity),
    ("2012 VP113", Element::SemiMajorAxis),
    ("2012 VP113", Element::Eccentricity),
    ("2004 VN112", Element::SemiMajorAxis),
    ("2010 GB174", Element::SemiMajorAxis),
    ("2010 GB174", Element::ArgPerihelion),
    ("2010 VZ98", Element::SemiMajorAxis),
];

#[test]
#[ignore = "hits the live JPL SBDB API"]
fn stable_kbo_table_matches_live_sbdb_modulo_documented_drift() {
    use p9_core::data::refresh::refresh_stable_kbos_from_sbdb;
    let client = SbdbClient::new().expect("HTTP client");
    let diffs = refresh_stable_kbos_from_sbdb(&client).expect("live SBDB refresh");
    let mut unexpected = 0;
    for d in &diffs {
        let known = STABLE_KBO_KNOWN_DRIFT
            .iter()
            .any(|(obj, el)| *obj == d.object && *el == d.element);
        eprintln!("{}: {d}", if known { "known drift" } else { "UNEXPECTED" });
        if !known {
            unexpected += 1;
        }
    }
    assert_eq!(
        unexpected, 0,
        "{unexpected} element(s) of stable_kbos drifted out of tolerance vs live SBDB          beyond the drift documented 2026-07-15 (transcription error or new          orbit-solution update; see stderr)"
    );
}
