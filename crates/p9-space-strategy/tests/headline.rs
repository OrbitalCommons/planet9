//! The strategy's conclusions, held to account.

use std::sync::OnceLock;

use p9_space_strategy::report::{build_with, Report, BUDGETS_H};
use p9_space_strategy::zones::Zone;

fn report() -> &'static Report {
    static R: OnceLock<Report> = OnceLock::new();
    R.get_or_init(|| build_with(150_000))
}

#[test]
fn the_ground_searches_have_found_about_three_quarters() {
    // Brown, Holman & Batygin (2024) quote 78% for ZTF + DES + Pan-STARRS1.
    let found = report().policies[0].found_by_ground;
    assert!(
        (0.65..0.85).contains(&found),
        "found from the ground: {found}"
    );
}

#[test]
fn conceding_more_to_rubin_leaves_less() {
    let p = &report().policies;
    assert!(p[0].unique > p[1].unique && p[1].unique > p[2].unique);
    assert!(p[0].conceded_to_rubin == 0.0);
    for row in p {
        let total = row.found_by_ground + row.conceded_to_rubin + row.unique;
        assert!((total - 1.0).abs() < 1e-9, "{}: {total}", row.policy);
    }
}

#[test]
fn more_time_never_buys_less_and_never_more_than_is_there() {
    for row in &report().policies {
        assert_eq!(row.captured.len(), BUDGETS_H.len());
        for w in row.captured.windows(2) {
            assert!(w[1] >= w[0]);
        }
        assert!(*row.captured.last().unwrap() <= row.unique + 1e-12);
    }
}

#[test]
fn the_reference_campaign_stays_out_of_rubins_sky() {
    let r = report();
    assert!(r.tiles.iter().all(|t| t.zone != Zone::RubinOverlap.label()));
    let zones: Vec<&str> = r.zones.iter().map(|z| z.zone).collect();
    for z in [
        Zone::AnticentreCrossing,
        Zone::CentreCrossing,
        Zone::NorthOfRubin,
    ] {
        assert!(
            zones.contains(&z.label()),
            "{} missing from the plan",
            z.label()
        );
    }
    let hours: f64 = r.zones.iter().map(|z| z.hours).sum();
    let captured: f64 = r.zones.iter().map(|z| z.captured).sum();
    assert!((hours - r.reference_hours).abs() < 1e-6);
    assert!((captured - r.reference_captured).abs() < 1e-9);
    // A season or two of telescope time, not a decade.
    assert!(
        (500.0..3000.0).contains(&r.reference_hours),
        "{} h",
        r.reference_hours
    );
}

#[test]
fn area_is_dear_and_depth_is_cheap() {
    // Planet Nine is V ≈ 19–23 and this telescope reaches V ≈ 23 in two
    // minutes: the plan is wide and shallow.
    let r = report();
    let short = r.tiles.iter().filter(|t| t.integration_s <= 300.0).count();
    assert!(short as f64 > 0.9 * r.tiles.len() as f64);
    assert!(r.reference_area_deg2 > 800.0);
}

#[test]
fn the_plan_survives_a_different_orbit() {
    for row in &report().robustness {
        let kept = row.reference_plan / row.own_plan;
        assert!(
            kept > 0.6,
            "{}: the reference plan keeps only {:.0}% of the best plan",
            row.prior,
            100.0 * kept
        );
    }
}

#[test]
fn the_winter_season_carries_most_of_the_campaign() {
    // Aphelion is near ecliptic longitude 70°, at opposition in December.
    let m = &report().hours_by_month;
    let winter: f64 = [10, 11, 0, 1, 2].iter().map(|&k| m[k]).sum();
    let total: f64 = m.iter().sum();
    assert!(winter > 0.7 * total, "winter {winter} of {total} h");
}

#[test]
fn distant_objects_are_a_by_product_not_a_programme() {
    let r = report();
    assert!(r.tno_plan_new > 20.0 && r.tno_plan_new < 500.0);
    assert!(r.tno_plan_distant_new < 10.0);
    assert!(r.deep_program_distant_new < 10.0);
}
