//! End-to-end offline acceptance tests (design/02): real-data fixtures for
//! seven objects spanning the class table, through parse → merge → sync →
//! battery.

use p9_rubin_watch::battery;
use p9_rubin_watch::ledger::{Ledger, Tier};
use p9_rubin_watch::sources::{merge_census, parse_mpc_distant, parse_sbdb_query};
use p9_rubin_watch::sync::sync;

fn fixture(name: &str) -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name);
    std::fs::read_to_string(p).expect("fixture readable")
}

fn seeded() -> Ledger {
    let mpc = parse_mpc_distant(&fixture("mpc_distant.json")).expect("mpc parses");
    let sbdb = parse_sbdb_query(&fixture("sbdb_query.json")).expect("sbdb parses");
    assert_eq!(mpc.len(), 7);
    assert_eq!(sbdb.len(), 7);
    let merged = merge_census(sbdb, mpc);
    // The join must collapse both sources to seven objects.
    assert_eq!(merged.len(), 7, "cross-source alias join failed");
    let (ledger, delta) = sync(Ledger::empty("t0"), &merged, "2026-07-16T00:00:00Z");
    assert_eq!(delta.new.len(), 7);
    ledger.check_invariants().expect("invariants");
    ledger
}

#[test]
fn seed_reproduces_expected_classes_and_tiers() {
    let ledger = seeded();

    let get = |d: &str| &ledger.objects[ledger.find_by_desig(d).expect(d)];

    // Sedna: vetted ETNO, sedna_like, aliases from both sources.
    let sedna = get("Sedna");
    assert_eq!(sedna.tier, Tier::Vetted);
    assert!(sedna.class_tags.contains(&"etno_vetted_class".into()));
    assert!(sedna.class_tags.contains(&"sedna_like".into()));
    assert!(sedna.aliases.contains(&"2003 VB12".into()));

    // 2025 LS2: ETNO class but provisional (U = 9 despite 2 oppositions).
    let ls2 = get("2025 LS2");
    assert_eq!(ls2.tier, Tier::Provisional);
    assert_eq!(ls2.class_tags, vec!["etno_vetted_class".to_string()]);
    assert!((ls2.elements.a_au - 523.0).abs() < 15.0);
    assert!((ls2.elements.varpi_deg - 336.0).abs() < 2.0);
    assert_eq!(ls2.quality.n_oppositions, Some(2));

    // Drac 2008 KV42: high inclination + retrograde, a ~ 42 (census via the
    // i-branch).
    let drac = get("2008 KV42");
    assert!(drac.class_tags.contains(&"high_inclination".into()));
    assert!(drac.class_tags.contains(&"retrograde".into()));

    // 2013 SY99: q = 49.90 sits a hair BELOW the 50 AU gap edge — an ETNO,
    // not a gap-band member (documents the razor's-edge case; the synthetic
    // gap object in classify's unit tests covers the in-band path).
    let sy99 = get("2013 SY99");
    assert!(sy99.class_tags.contains(&"etno_vetted_class".into()));
    assert!(
        !sy99.class_tags.contains(&"gap_band".into()),
        "{:?}",
        sy99.class_tags
    );

    // Leleakuhonua 2015 TG387: ETNO + sedna_like (q ~ 65).
    let lele = get("2015 TG387");
    assert!(lele.class_tags.contains(&"etno_vetted_class".into()));
}

#[test]
fn second_sync_on_same_fixtures_is_a_noop() {
    let ledger = seeded();
    let mpc = parse_mpc_distant(&fixture("mpc_distant.json")).unwrap();
    let sbdb = parse_sbdb_query(&fixture("sbdb_query.json")).unwrap();
    let before = ledger.to_canonical_json();
    let (ledger2, delta) = sync(ledger, &merge_census(sbdb, mpc), "2026-07-17T00:00:00Z");
    assert!(delta.is_empty(), "{delta:?}");
    assert!(p9_rubin_watch::ledger::content_equal(
        &serde_json::from_str(&before).unwrap(),
        &ledger2
    ));
}

#[test]
fn battery_runs_on_the_seeded_ledger() {
    let ledger = seeded();
    let out = battery::run(&ledger, &[]);
    // Three vetted ETNOs in the fixture set (Sedna, GB174, Leleakuhonua) —
    // enough for stats; LS2 joins in the provisional column.
    let vetted = out.vetted.expect("vetted stats");
    let both = out.with_provisional.expect("both stats");
    assert_eq!(both.n, vetted.n + 1);
    // LS2 at varpi 336 deg sits off the cluster: including it lowers R-bar.
    assert!(both.r_bar < vetted.r_bar);
    // Leleakuhonua (q = 64.7) sits inside the 50-65 AU band and is surfaced
    // for curation; SY99 (q = 49.90) misses the low edge by 0.1 AU.
    assert_eq!(out.gap.new_gap_band, vec!["2015 TG387".to_string()]);
    // Fixture subset is far smaller than BROWN_2017 -> staleness flags.
    assert!(out.posterior_stale);
    // Frozen baseline reproduces the crate's published-scale numbers.
    assert!(out.baseline_frozen.n >= 10);
    assert!(out.baseline_frozen.r_bar > 0.4);
}

#[test]
fn report_renders_with_frontmatter() {
    let ledger = seeded();
    let mpc = parse_mpc_distant(&fixture("mpc_distant.json")).unwrap();
    let sbdb = parse_sbdb_query(&fixture("sbdb_query.json")).unwrap();
    // Re-sync a fresh ledger so the delta is the seed itself.
    let (ledger2, delta) = sync(Ledger::empty("t0"), &merge_census(sbdb, mpc), "t1");
    drop(ledger);
    let battery = battery::run(&ledger2, &[]);
    let md =
        p9_rubin_watch::report::render("2026-07-16T00:00:00Z", &ledger2, &delta, &battery, &[]);
    assert!(md.starts_with("---\n"));
    assert!(md.contains("schema_version: 1"));
    assert!(md.contains("## Battery"));
    assert!(md.contains("2025 LS2"));
    assert!(md.contains("Interpretation guardrails"));
    // Frontmatter battery blob is machine-parseable JSON.
    let line = md
        .lines()
        .find(|l| l.starts_with("battery: "))
        .expect("battery line");
    let parsed: serde_json::Value =
        serde_json::from_str(line.trim_start_matches("battery: ")).expect("battery json");
    assert!(parsed.get("gap").is_some());
}
