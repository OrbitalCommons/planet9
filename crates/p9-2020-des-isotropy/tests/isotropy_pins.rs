//! Regression pins for the DES isotropy reproduction and synthetic controls.
//!
//! The central reproduced result (Bernardinelli et al. 2020): the DES eTNO
//! orbital angles look clustered against a *flat* null, but folding an
//! isotropic population through the DES footprint selection shows the data are
//! consistent with isotropy. We pin this for the longitude of perihelion ϖ
//! (the Planet 9 alignment angle) and document the residual Ω clustering the
//! static-footprint model cannot fully account for.

use p9_2020_des_isotropy::analysis::{angle_p_values, run_battery, summarize};
use p9_2020_des_isotropy::des_sample::{case_sample, Angle, SampleCase};
use p9_2020_des_isotropy::selection::NullModel;
use p9_2020_des_isotropy::statistics::{
    circular_correlation, mc_correlation_p_value, IsotropySuite,
};
use p9_core::constants::TWO_PI;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

const SEED: u64 = 2020;
const MC: usize = 8000;

/// Headline reproduction: longitude of perihelion ϖ — the Planet 9 alignment
/// angle — is consistent with isotropy (every test p > 0.05) once the DES
/// footprint selection is accounted for, in the three larger Cases (n ≥ 4).
///
/// Case 4 (n = 3, a > 250 AU & q > 37 AU) is excluded: its three objects have
/// ϖ ≈ {7.7°, 0.2°, 349.4°}, a near-coincidence that yields R̄ ≈ 0.99, and no
/// selection model removes a 3-point pile-up — every test lands at p ≈ 0.03,
/// borderline (and right where the paper's own only-significant cells, the
/// a > 250 AU Ω tests, landed at p = 0.03). Treated as the small-n borderline.
#[test]
fn varpi_consistent_with_isotropy_under_selection() {
    for case in [SampleCase::Case1, SampleCase::Case2, SampleCase::Case3] {
        let sample = case_sample(case);
        let p = angle_p_values(&sample, Angle::Varpi, NullModel::DesSelection, SEED, MC);
        assert_eq!(
            p.n_consistent(0.05),
            4,
            "{case:?} ϖ: only {} of 4 tests consistent (p = {:?})",
            p.n_consistent(0.05),
            p.p_values()
        );
    }
    // Case 4 (n = 3): borderline, p ≈ 0.03 across the suite — documented, not
    // forced to isotropy.
    let c4 = angle_p_values(
        &case_sample(SampleCase::Case4),
        Angle::Varpi,
        NullModel::DesSelection,
        SEED,
        MC,
    );
    assert!(c4.r_bar > 0.95, "Case4 ϖ R̄ = {:.3}", c4.r_bar);
    for p in c4.p_values() {
        assert!(
            (0.01..0.05).contains(&p),
            "Case4 ϖ p = {p:.4} (expected borderline)"
        );
    }
}

/// The DES selection null raises the p-values relative to the flat null
/// (footprint-induced clustering is partly explained away): the count of
/// significant tests strictly decreases.
#[test]
fn selection_null_reduces_significance() {
    let flat = summarize(&run_battery(NullModel::Uniform, SEED, MC), 0.05);
    let sel = summarize(&run_battery(NullModel::DesSelection, SEED, MC), 0.05);
    assert!(
        sel.n_significant < flat.n_significant,
        "selection {} should be < flat {}",
        sel.n_significant,
        flat.n_significant
    );
    // The minimum p-value also rises (the most extreme cell is softened).
    assert!(sel.min_p >= flat.min_p);
}

/// Under the selection null, the great majority of the per-angle test results
/// for ω and ϖ (the apsidal angles) are consistent with isotropy. (Ω remains
/// the most clustered — see `omega_node_is_the_residual` — matching the paper,
/// where the only 2-of-12 significant Kuiper tests were both Ω at a > 250 AU.)
#[test]
fn apsidal_angles_mostly_consistent_under_selection() {
    let mut consistent = 0usize;
    let mut total = 0usize;
    for case in SampleCase::all() {
        let sample = case_sample(case);
        for angle in [Angle::ArgPeri, Angle::Varpi] {
            let p = angle_p_values(&sample, angle, NullModel::DesSelection, SEED, MC);
            consistent += p.n_consistent(0.05);
            total += 4;
        }
    }
    let frac = consistent as f64 / total as f64;
    assert!(
        frac > 0.5,
        "only {frac:.2} of ω/ϖ tests consistent ({consistent}/{total})"
    );
}

/// Documented residual: the ascending node Ω is the most concentrated angle
/// (R̄ ≈ 0.79–0.99) and our static-footprint selection cannot fully explain it
/// away — several Ω tests stay significant. This mirrors the paper's finding
/// that Ω (not ϖ) was the only angle to trip any test, but the paper's
/// per-exposure selection model removes more of the effect than our static
/// RA/Dec footprint does. We pin the *direction* of the residual (Ω is the
/// least-isotropic angle), not isotropy for Ω.
#[test]
fn omega_node_is_the_residual() {
    let sample = case_sample(SampleCase::Case1);
    let node = angle_p_values(&sample, Angle::Node, NullModel::DesSelection, SEED, MC);
    let varpi = angle_p_values(&sample, Angle::Varpi, NullModel::DesSelection, SEED, MC);
    // Ω has the higher R̄ and the smaller (more significant) Rayleigh p.
    assert!(
        node.r_bar > varpi.r_bar,
        "R̄: Ω {} ϖ {}",
        node.r_bar,
        varpi.r_bar
    );
    assert!(
        node.rayleigh_p < varpi.rayleigh_p,
        "Ω p {} should be < ϖ p {}",
        node.rayleigh_p,
        varpi.rayleigh_p
    );
}

/// Power check: a strongly clustered synthetic sample yields small p across
/// every statistic in the suite (the tests can detect clustering).
#[test]
fn clustered_control_yields_small_p() {
    let clustered: Vec<f64> = (0..7).map(|k| 1.0 + 0.05 * (k as f64 - 3.0)).collect();
    let suite = IsotropySuite::compute(&clustered, 123, MC);
    assert!(suite.r_bar > 0.9, "R̄ = {:.3}", suite.r_bar);
    for (name, p) in [
        ("rayleigh", suite.rayleigh_p),
        ("kuiper", suite.kuiper_p),
        ("watson", suite.watson_p),
        ("beran", suite.beran_p),
    ] {
        assert!(p < 0.05, "{name} p = {p:.4} on clustered control");
    }
}

/// Calibration: averaged over many genuinely-uniform samples the fraction of
/// analytic tests rejecting isotropy at 0.05 is small (not pathological).
#[test]
fn uniform_control_is_calibrated() {
    let trials = 600usize;
    let n = 7usize;
    let mut rng = StdRng::seed_from_u64(7);
    let mut rayleigh_rej = 0usize;
    let mut kuiper_rej = 0usize;
    for _ in 0..trials {
        let sample: Vec<f64> = (0..n).map(|_| rng.gen::<f64>() * TWO_PI).collect();
        if p9_core::analysis::circular::rayleigh_p_value(&sample) < 0.05 {
            rayleigh_rej += 1;
        }
        if p9_core::analysis::circular::kuiper_p_value(&sample) < 0.05 {
            kuiper_rej += 1;
        }
    }
    assert!(rayleigh_rej as f64 / (trials as f64) < 0.2);
    assert!(kuiper_rej as f64 / (trials as f64) < 0.2);
}

/// The MC Watson null is calibrated: a uniform sample's MC p-value is not
/// systematically small.
#[test]
fn mc_null_calibrated_on_uniform() {
    let trials = 300usize;
    let n = 7usize;
    let mut rng = StdRng::seed_from_u64(99);
    let mut watson_rej = 0usize;
    for t in 0..trials {
        let sample: Vec<f64> = (0..n).map(|_| rng.gen::<f64>() * TWO_PI).collect();
        let suite = IsotropySuite::compute(&sample, 1000 + t as u64, 1500);
        if suite.watson_p < 0.05 {
            watson_rej += 1;
        }
    }
    assert!(watson_rej as f64 / (trials as f64) < 0.2);
}

/// Two-angle correlation: ω vs Ω in the DES sample shows no strong joint
/// structure and the permutation p-value is not significant.
#[test]
fn des_two_angle_correlation_not_significant() {
    let sample = case_sample(SampleCase::Case1);
    let omega: Vec<f64> = sample.iter().map(|o| o.argp()).collect();
    let node: Vec<f64> = sample.iter().map(|o| o.node()).collect();
    let r = circular_correlation(&omega, &node);
    assert!(r.abs() < 0.95, "|r_T| = {:.3}", r.abs());
    let p = mc_correlation_p_value(&omega, &node, SEED, MC);
    assert!(
        p > 0.05,
        "ω–Ω correlation p = {p:.4} (unexpectedly significant)"
    );
}

/// Correlation power: a perfectly coupled angle pair has |r_T| ≈ 1 and a tiny
/// permutation p-value.
#[test]
fn correlation_detects_coupling() {
    let a: Vec<f64> = (0..7).map(|k| 0.3 + 0.7 * k as f64).collect();
    let b: Vec<f64> = a.iter().map(|x| (x + 0.1).rem_euclid(TWO_PI)).collect();
    let r = circular_correlation(&a, &b);
    assert!(r.abs() > 0.8, "coupled |r_T| = {:.3}", r.abs());
    let p = mc_correlation_p_value(&a, &b, 5, MC);
    assert!(p < 0.05, "coupled correlation p = {p:.4}");
}
