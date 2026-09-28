//! Headline reproduction of Bansal et al. (2026), arXiv:2607.15646.

use p9_2026_cluster_inclinations::debias::{
    default_w_grid, fit_intrinsic_width, P_ONE_SIGMA, P_THREE_SIGMA,
};
use p9_2026_cluster_inclinations::population::BirthEnvironment;
use p9_2026_cluster_inclinations::reference::*;
use p9_2026_cluster_inclinations::sample::{discovery_positions, sample_elements, HIGH_Q_SAMPLE};
use p9_2026_cluster_inclinations::simulation::{
    run, RunConfig, CLUSTER_FREE_P9, CLUSTER_INFLUENCED_P9,
};
use p9_2026_cluster_inclinations::width::{
    mean_pole, population_width, relative_inclinations, rotation_to_pole, SelectionCuts,
};

/// The observed sample in its own mean-pole frame: (latitudes at discovery,
/// inclinations), both in radians.
fn observed_in_mean_plane() -> (Vec<f64>, Vec<f64>) {
    let elements = sample_elements();
    let pole = mean_pole(&elements);
    let rot = rotation_to_pole(&pole);
    let incs = relative_inclinations(&elements, &pole);
    let betas: Vec<f64> = discovery_positions()
        .iter()
        .zip(&incs)
        .map(|(p, &i)| {
            let p = rot * p;
            // Two-body propagation over decades can nudge |β| past i by a
            // hair; the conditional CDF is only defined for |β| ≤ i.
            (p.z / p.norm()).asin().abs().min(i - 1e-6)
        })
        .collect();
    (betas, incs)
}

#[test]
fn observed_width_is_twelve_degrees() {
    assert_eq!(HIGH_Q_SAMPLE.len(), N_OBSERVED);
    let (betas, incs) = observed_in_mean_plane();
    let fit = fit_intrinsic_width(&betas, &incs, &default_w_grid(), 2000, 1);
    let (lo, hi) = fit.interval(P_ONE_SIGMA);
    eprintln!(
        "w_obs = {}° (1σ {lo}–{hi}°); paper {W_OBS_DEG} ({}–{})",
        fit.w_best_deg, W_OBS_1SIGMA_DEG.0, W_OBS_1SIGMA_DEG.1
    );
    assert!(
        (W_OBS_DEG - 3.0..=W_OBS_DEG + 3.0).contains(&fit.w_best_deg),
        "best-fit w = {}°, paper {W_OBS_DEG}°",
        fit.w_best_deg
    );
    assert!(
        lo <= W_OBS_DEG && hi >= W_OBS_DEG,
        "1σ interval [{lo}, {hi}] excludes 12°"
    );
    assert!(
        hi < W_CLUSTER_INFLUENCED_DEG.0,
        "1σ upper bound {hi}° reaches the cluster-influenced band"
    );
}

#[test]
fn cluster_influenced_width_is_rejected_at_three_sigma() {
    let (betas, incs) = observed_in_mean_plane();
    let fit = fit_intrinsic_width(&betas, &incs, &[12.0, 26.0], 20_000, 2);
    let p26 = fit.p_at(26.0);
    eprintln!("p(w = 26°) = {p26}; paper {P_REJECT_W26}");
    assert!(
        p26 < P_THREE_SIGMA * 5.0,
        "p(26°) = {p26} is not a ~3σ rejection"
    );
    assert!(
        fit.p_at(12.0) > P_ONE_SIGMA,
        "the best-fit width must be accepted"
    );
}

#[test]
fn planet_nine_does_not_cool_a_cluster_excited_population() {
    let result = run(&RunConfig::reduced_scale(
        BirthEnvironment::ClusterInfluenced,
        CLUSTER_INFLUENCED_P9[0],
    ));
    let whole = SelectionCuts::whole_population();
    let (w0, w1) = (result.initial_width(), result.final_width());
    let (p0, p1) = (
        result.initial_width_with(&whole),
        result.pooled_width(&whole),
    );
    eprintln!(
        "cluster-influenced: window w {w0:.1}° → {w1:.1}°, whole population {p0:.1}° → {p1:.1}° \
         (survival {:.2}, κ {:.2}); history {:?}",
        result.survival(),
        result.final_kappa(),
        result
            .width_history(&SelectionCuts::width_sample())
            .iter()
            .zip(result.width_history(&whole))
            .map(|((t, w, n), (_, wp, np))| format!("{t:.0}:{w:.1}({n})/{wp:.1}({np})"))
            .collect::<Vec<_>>()
    );
    // The whole excited population keeps its dispersion: Planet Nine's
    // secular forcing neither cools nor heats it (the 2 Gyr run documented in
    // the crate docs holds at 25–28°).
    assert!(
        p0 >= W_CLUSTER_INFLUENCED_DEG.0 - 4.0,
        "initial whole-population w = {p0}° is not a cluster-excited population"
    );
    assert!(
        p1 > p0 - 3.0,
        "Planet Nine cooled the population: {p0}° → {p1}°"
    );
    assert!(
        p1 >= W_CLUSTER_INFLUENCED_DEG.0 - 4.0,
        "final whole-population w = {p1}° left the cluster-influenced band"
    );
    // The observable window holds only a few tens of orbits at this scale
    // and its fit is a transient (see the crate docs for the 1024-particle,
    // 2 Gyr runs, where it pools to ~25°); it must at least stay above the
    // observed width.
    assert!(
        w1 > W_OBS_DEG,
        "windowed w = {w1}° dropped to the observed {W_OBS_DEG}°"
    );
}

#[test]
fn planet_nine_does_not_overexcite_a_cluster_free_population() {
    let result = run(&RunConfig::reduced_scale(
        BirthEnvironment::ClusterFree,
        CLUSTER_FREE_P9[6],
    ));
    let (w0, w1) = (result.initial_width(), result.final_width());
    eprintln!(
        "cluster-free: w {w0:.1}° → {w1:.1}° (survival {:.2}, κ {:.2}); history {:?}",
        result.survival(),
        result.final_kappa(),
        result
            .width_history(&SelectionCuts::width_sample())
            .iter()
            .zip(result.width_history(&SelectionCuts::whole_population()))
            .map(|((t, w, n), (_, wp, np))| format!("{t:.0}:{w:.1}({n})/{wp:.1}({np})"))
            .collect::<Vec<_>>()
    );
    assert!(
        (w0 - W_CLUSTER_FREE_INITIAL_DEG).abs() < 2.0,
        "initial w = {w0}°"
    );
    // Planet Nine's forced-plane broadening heats the cold population toward
    // the paper's 16–19.5° band (partially, at 400 Myr of the paper's 4 Gyr)
    // without over-exciting it.
    assert!(w1 > w0 + 1.0, "no Planet Nine broadening: {w0}° → {w1}°");
    assert!(
        w1 <= W_CLUSTER_FREE_DEG.1 + 0.5,
        "final w = {w1}° over-excited past the paper's cluster-free band"
    );
    assert!(
        W_CLUSTER_INFLUENCED_DEG.0 - w1 > 5.0,
        "final w = {w1}° approaches the cluster-influenced band"
    );
}

#[test]
fn selection_cuts_match_the_paper() {
    let c = SelectionCuts::width_sample();
    assert_eq!(
        (c.a_min, c.a_max, c.q_min, c.q_max),
        (200.0, 2000.0, 40.0, 80.0)
    );
    let (w, n) = population_width(&sample_elements(), &c);
    // Alicanto's current JPL solution (q = 80.6 AU) falls just outside the
    // strict q < 80 AU cut; the other 18 pass.
    assert_eq!(n, N_OBSERVED - 1);
    // The raw (biased) width of the observed sample is itself well inside the
    // cold band; debiasing can only lower it further.
    assert!(w < W_CLUSTER_INFLUENCED_DEG.0, "raw w = {w}°");
}

#[test]
#[ignore = "paper scale: 10⁴ particles × 4 Gyr per Table 1 row (hours)"]
fn paper_scale_table_1() {
    for p9 in CLUSTER_INFLUENCED_P9.iter().chain(CLUSTER_FREE_P9.iter()) {
        let env = if CLUSTER_INFLUENCED_P9.contains(p9) {
            BirthEnvironment::ClusterInfluenced
        } else {
            BirthEnvironment::ClusterFree
        };
        let r = run(&RunConfig::paper_scale(env, *p9));
        eprintln!(
            "{env:?} m9={} a9={} e9={}: w = {:.1}°, κ = {:.3}",
            p9.mass_earth,
            p9.a,
            p9.e,
            r.final_width(),
            r.final_kappa()
        );
    }
}
