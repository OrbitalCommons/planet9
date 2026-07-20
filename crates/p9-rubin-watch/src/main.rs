//! `p9-rubin-watch` CLI. Phase 1: the `mpc-sync` subcommand
//! (rubin_watch/design/02).
//!
//! Exit codes: 0 = clean/no-op, 1 = delta emitted, 2 = source failure.

use p9_rubin_watch::{battery, clock, ledger::Ledger, report, sources, sync};

struct Args {
    dry_run: bool,
    skip_battery: bool,
    ledger_path: std::path::PathBuf,
    reports_dir: std::path::PathBuf,
    fixtures: Option<std::path::PathBuf>,
}

fn usage() -> ! {
    eprintln!(
        "usage: p9-rubin-watch mpc-sync [--dry-run] [--skip-battery] \
         [--ledger PATH] [--reports DIR] [--fixtures DIR]\n\
       p9-rubin-watch coverage [--map PATH]\n\
         \n\
         mpc-sync: --fixtures DIR reads DIR/mpc_distant.json and \
         DIR/sbdb_query.json instead of the network.\n\
         coverage: validate + summarize the LSST coverage map (schema gate)."
    );
    std::process::exit(2)
}

fn coverage_main(mut argv: impl Iterator<Item = String>) -> ! {
    let mut map_path = std::path::PathBuf::from("rubin_watch/lsst_coverage.json");
    while let Some(a) = argv.next() {
        match a.as_str() {
            "--map" => map_path = argv.next().unwrap_or_else(|| usage()).into(),
            _ => usage(),
        }
    }
    let map = match p9_core::data::lsst_coverage::CoverageMap::load(&map_path) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("coverage map INVALID: {e}");
            std::process::exit(2);
        }
    };
    let npix = p9_core::coords::healpix::npix(map.nside);
    let pix_area =
        4.0 * std::f64::consts::PI / npix as f64 * (180.0 / std::f64::consts::PI).powi(2);
    println!(
        "coverage map OK: nside {}, window {}..{}, {} visits, source {:?}",
        map.nside, map.first_night, map.last_night, map.n_reconstructed_visits, map.source
    );
    for (band, b) in &map.bands {
        let linkable = b.n_visits.iter().filter(|&&n| n >= 3).count();
        println!(
            "  {band}: {} covered pixels ({:.0} deg2), {} linkable (>=3 visits, {:.0} deg2), fiducial {:.1}",
            b.pixels.len(),
            b.pixels.len() as f64 * pix_area,
            linkable,
            linkable as f64 * pix_area,
            map.fiducial_depth[band]
        );
    }
    println!(
        "  flags: {} crowding pixels, {} template-epoch pixels",
        map.crowding_pixels.len(),
        map.template_epoch_pixels.len()
    );
    std::process::exit(0)
}

fn parse_args() -> Args {
    let mut argv = std::env::args().skip(1);
    match argv.next().as_deref() {
        Some("mpc-sync") => {}
        Some("coverage") => coverage_main(argv),
        _ => usage(),
    }
    let mut args = Args {
        dry_run: false,
        skip_battery: false,
        ledger_path: "rubin_watch/etno_ledger.json".into(),
        reports_dir: "reports/rubin-watch".into(),
        fixtures: None,
    };
    while let Some(a) = argv.next() {
        match a.as_str() {
            "--dry-run" => args.dry_run = true,
            "--skip-battery" => args.skip_battery = true,
            "--ledger" => args.ledger_path = argv.next().unwrap_or_else(|| usage()).into(),
            "--reports" => args.reports_dir = argv.next().unwrap_or_else(|| usage()).into(),
            "--fixtures" => args.fixtures = Some(argv.next().unwrap_or_else(|| usage()).into()),
            _ => usage(),
        }
    }
    args
}

fn main() {
    let args = parse_args();
    let (now_iso, today) = clock::now_utc();

    // --- census ------------------------------------------------------------
    let mut degraded: Vec<String> = Vec::new();
    let (mpc, mpc_stamp, sbdb) = match &args.fixtures {
        Some(dir) => {
            let read = |name: &str| {
                std::fs::read_to_string(dir.join(name))
                    .unwrap_or_else(|e| panic!("fixture {name}: {e}"))
            };
            (
                sources::parse_mpc_distant(&read("mpc_distant.json")).expect("mpc fixture"),
                Some("fixture".to_string()),
                sources::parse_sbdb_query(&read("sbdb_query.json")).expect("sbdb fixture"),
            )
        }
        None => {
            let (mpc, stamp) = match sources::fetch_mpc_distant() {
                Ok((m, s)) => (m, s),
                Err(e) => {
                    eprintln!("mpc source failed: {e}");
                    degraded.push("mpc".into());
                    (Vec::new(), None)
                }
            };
            let sbdb = match sources::fetch_sbdb_census() {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("sbdb source failed: {e}");
                    degraded.push("sbdb".into());
                    Vec::new()
                }
            };
            (mpc, stamp, sbdb)
        }
    };
    if mpc.is_empty() && sbdb.is_empty() {
        eprintln!("both census sources failed; ledger untouched");
        std::process::exit(2);
    }
    eprintln!(
        "census: {} sbdb + {} mpc records{}",
        sbdb.len(),
        mpc.len(),
        if degraded.is_empty() {
            String::new()
        } else {
            format!(" (degraded: {})", degraded.join(","))
        }
    );
    let merged = sources::merge_census(sbdb, mpc);

    // --- sync ----------------------------------------------------------------
    let prior = if args.ledger_path.exists() {
        match Ledger::load(&args.ledger_path) {
            Ok(l) => l,
            Err(e) => {
                eprintln!("ledger load failed: {e}");
                std::process::exit(2);
            }
        }
    } else {
        eprintln!("no ledger at {:?}; seeding fresh", args.ledger_path);
        Ledger::empty(&now_iso)
    };
    let seed_run = prior.objects.is_empty();
    let (mut ledger, delta) = sync::sync(prior, &merged, &now_iso);
    ledger.source_stamps.mpc_distant_last_modified = mpc_stamp;
    ledger.source_stamps.sbdb_query_fetched = Some(now_iso.clone());
    ledger.check_invariants().expect("post-sync invariants");

    eprintln!(
        "delta: {} new, {} promoted, {} demoted, {} reclassified, {} renamed, {} retired, {} discrepant, {} refreshed",
        delta.new.len(),
        delta.promoted.len(),
        delta.demoted.len(),
        delta.reclassified.len(),
        delta.renamed.len(),
        delta.retired.len(),
        delta.discrepant.len(),
        delta.refreshed
    );

    // --- battery + report -------------------------------------------------------
    let emit_report = !delta.is_empty() || seed_run;
    if emit_report && !args.skip_battery {
        let new_ids: Vec<String> = delta
            .new
            .iter()
            .filter_map(|d| ledger.find_by_desig(d))
            .map(|i| ledger.objects[i].id.clone())
            .collect();
        // On a seed run everything is "new"; a without-new column would be
        // empty and meaningless, so pass no exclusions.
        let battery = battery::run(&ledger, if seed_run { &[] } else { &new_ids });
        let md = report::render(&now_iso, &ledger, &delta, &battery, &degraded);
        if args.dry_run {
            println!("{md}");
        } else {
            std::fs::create_dir_all(&args.reports_dir).expect("mkdir reports");
            let path = args.reports_dir.join(format!("{today}-mpc.md"));
            std::fs::write(&path, md).expect("write report");
            eprintln!("report: {path:?}");
        }
    }

    if !args.dry_run {
        match ledger.save_if_changed(&args.ledger_path) {
            Ok(true) => eprintln!("ledger written: {:?}", args.ledger_path),
            Ok(false) => eprintln!("ledger unchanged"),
            Err(e) => {
                eprintln!("ledger write failed: {e}");
                std::process::exit(2);
            }
        }
    }

    std::process::exit(if delta.is_empty() { 0 } else { 1 });
}
