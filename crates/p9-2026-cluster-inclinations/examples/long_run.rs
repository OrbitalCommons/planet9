//! Long integration driver: `long_run <cluster|free> <n> <t_myr> [q_max_au] [seed]`.
//! Prints the windowed and whole-population width history plus the pooled
//! Table 1 statistics over the final quarter of the run.

use p9_2026_cluster_inclinations::population::BirthEnvironment;
use p9_2026_cluster_inclinations::simulation::{
    run, RunConfig, CLUSTER_FREE_P9, CLUSTER_INFLUENCED_P9,
};
use p9_2026_cluster_inclinations::width::SelectionCuts;
use p9_core::constants::YEAR_DAYS;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let env = match args.get(1).map(String::as_str) {
        Some("cluster") => BirthEnvironment::ClusterInfluenced,
        Some("free") => BirthEnvironment::ClusterFree,
        _ => panic!("usage: long_run <cluster|free> <n> <t_myr> [q_max_au] [seed]"),
    };
    let n: usize = args[2].parse().unwrap();
    let t_myr: f64 = args[3].parse().unwrap();
    let p9 = match env {
        BirthEnvironment::ClusterInfluenced => CLUSTER_INFLUENCED_P9[0],
        BirthEnvironment::ClusterFree => CLUSTER_FREE_P9[6],
    };
    let mut cfg = RunConfig::reduced_scale(env, p9);
    cfg.n_particles = n;
    cfg.t_days = t_myr * 1.0e6 * YEAR_DAYS;
    cfg.analysis_from_days = 0.75 * cfg.t_days;
    if let Some(q) = args.get(4) {
        cfg.q_max = q.parse().unwrap();
    }
    if let Some(s) = args.get(5) {
        cfg.seed = s.parse().unwrap();
    }
    let r = run(&cfg);
    for ((t, w, n), (_, wp, np)) in r
        .width_history(&SelectionCuts::width_sample())
        .iter()
        .zip(r.width_history(&SelectionCuts::whole_population()))
    {
        println!("t = {t:6.0} Myr  window w = {w:5.1} deg (n = {n:4})  whole w = {wp:5.1} deg (n = {np:4})");
    }
    println!(
        "pooled: w = {:.1} deg, kappa = {:.3}, survival = {:.3}",
        r.final_width(),
        r.final_kappa(),
        r.survival()
    );
}
