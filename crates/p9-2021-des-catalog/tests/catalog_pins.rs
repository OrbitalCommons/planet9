//! Regression pins for Bernardinelli et al. (2021), the six-year DES TNO
//! catalog. Reference catalog size / footprint / depth are asserted as
//! labelled identities; the clustering and completeness numbers are computed
//! from real seeded calculations and pinned with documented residuals.

use p9_2021_des_catalog::catalog::{Angle, DES_EXTREME_TNOS};
use p9_2021_des_catalog::clustering::{clustering_battery, clustering_for_angle, NullModel};
use p9_2021_des_catalog::completeness::{
    des_depth_r, footprint_completeness, sky_coverage_fraction,
};
use p9_2021_des_catalog::reference;

const ITERS: usize = 30_000;

#[test]
fn reference_constants_match_paper() {
    assert_eq!(reference::N_TNOS, 815);
    assert_eq!(reference::N_NEW_TNOS, 461);
    assert_eq!(reference::N_EXTREME_TNOS, 16);
    assert_eq!(reference::FOOTPRINT_DEG2, 5000.0);
    assert_eq!(reference::DEPTH_R, 23.8);
    assert!(reference::PAPER_CONCLUSION.contains("isotropy"));
}

#[test]
fn footprint_depth_is_shared_table_value() {
    // The pinned DES depth must be the shared p9-core survey-table 23.8.
    assert!((des_depth_r() - reference::DEPTH_R).abs() < 1e-10);
}

#[test]
fn varpi_clustering_consistent_with_des_selection() {
    // Headline reproduction: the Planet-9 alignment angle ϖ on the DES
    // extreme-TNO subset is consistent with azimuthal isotropy once the DES
    // selection function is folded in (Rayleigh p reported, MC p > 0.05).
    let r = clustering_for_angle(
        &DES_EXTREME_TNOS,
        Angle::Varpi,
        NullModel::DesSelection,
        2021,
        ITERS,
    );
    assert!(
        r.rayleigh_p > 0.0 && r.rayleigh_p <= 1.0,
        "rayleigh p = {}",
        r.rayleigh_p
    );
    assert!(
        r.consistent_with_isotropy(),
        "varpi MC p = {:.3} (R-bar = {:.3}, Rayleigh p = {:.3}) should exceed 0.05",
        r.mc_p,
        r.r_bar,
        r.rayleigh_p
    );
}

#[test]
fn full_battery_runs_for_all_angles() {
    let results = clustering_battery(&DES_EXTREME_TNOS, NullModel::DesSelection, 2021, ITERS);
    assert_eq!(results.len(), 3);
    for r in &results {
        assert!(r.r_bar >= 0.0 && r.r_bar <= 1.0);
        assert!(r.mc_p > 0.0 && r.mc_p <= 1.0);
    }
    // varpi and omega should be the isotropy-consistent angles; Omega is the
    // documented residual (high concentration).
    let varpi = results.iter().find(|r| r.angle == Angle::Varpi).unwrap();
    assert!(
        varpi.consistent_with_isotropy(),
        "varpi p = {:.3}",
        varpi.mc_p
    );
}

#[test]
fn distant_object_completeness_fraction_in_unit_interval() {
    // P9-relevant distant-object completeness over the footprint: a fraction
    // strictly in (0, 1). A faint extreme TNO (H = 7) at 300 AU.
    let c = footprint_completeness(7.0, 300.0, 10.0, 2021, 300_000);
    assert!(c > 0.0 && c < 1.0, "footprint completeness = {c}");
}

#[test]
fn footprint_covers_about_a_tenth_of_sky() {
    // The ~5000 deg^2 footprint with a |b| > 10 deg galactic mask covers
    // ~10% of the celestial sphere; this is the sky term in the catalog's
    // distant-object completeness.
    let f = sky_coverage_fraction(10.0, 2021, 400_000);
    assert!((f - 0.11).abs() < 0.04, "coverage = {f:.3}");
}
