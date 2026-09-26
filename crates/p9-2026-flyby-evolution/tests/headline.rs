//! Headline reproduction of Pfalzner, Wagner & Bischoff (2026),
//! arXiv:2609.03575.

use p9_2026_flyby_evolution::colours::colour_dichotomy;
use p9_2026_flyby_evolution::evolution::{
    evolve, in_group, near_neptune, EvolutionConfig, EvolutionResult,
};
use p9_2026_flyby_evolution::flyby::{run_flyby, FlybyConfig, FlybyOutcome};
use p9_2026_flyby_evolution::groups::{Census, DynamicalGroup, GROUPS};
use p9_2026_flyby_evolution::reference::*;
use p9_core::constants::RAD2DEG;
use std::sync::OnceLock;

fn flyby() -> &'static FlybyOutcome {
    static OUT: OnceLock<FlybyOutcome> = OnceLock::new();
    OUT.get_or_init(|| run_flyby(&FlybyConfig::model_a1(4096)))
}

fn evolved() -> &'static EvolutionResult {
    static OUT: OnceLock<EvolutionResult> = OnceLock::new();
    OUT.get_or_init(|| evolve(&flyby().bound(), &EvolutionConfig::reduced_scale()))
}

fn report(census: &Census, label: &str) {
    let f = census.fractions();
    eprintln!("{label}: N = {} orbits with p <= 100 AU", census.n);
    for (k, g) in GROUPS.iter().enumerate() {
        eprintln!(
            "  {:<11} {:.3}  (paper t=0 {:.3})",
            g.label(),
            f[k],
            TABLE1_FRACTION_T0[k]
        );
    }
}

#[test]
fn flyby_produces_every_table_1_family() {
    let out = flyby();
    let census = out.census();
    report(&census, "post-flyby");
    eprintln!(
        "unbound from the Sun: {:.3} (paper captures {CAPTURED_BY_PERTURBER})",
        out.unbound_fraction()
    );
    for g in GROUPS {
        assert!(census.weight(g) > 0.0, "{} is empty", g.label());
    }
    // The encounter strips a few tenths of the (1/r-weighted) disc: the
    // paper's 8.3% capture plus ejections.
    assert!(
        (0.1..0.35).contains(&out.unbound_fraction()),
        "unbound fraction {}",
        out.unbound_fraction()
    );
    // The two belts dominate the census just after the flyby, as in Table 1
    // (0.185 + 0.231 ≈ 42%); Sedna-like and injected tracers follow at the
    // ten-percent level.
    let belts = census.fraction(DynamicalGroup::ColdKb) + census.fraction(DynamicalGroup::HotKb);
    assert!((0.3..0.55).contains(&belts), "belt fraction {belts}");
    assert!((0.05..0.2).contains(&census.fraction(DynamicalGroup::SednaLike)));
    assert!((0.1..0.3).contains(&census.fraction(DynamicalGroup::Inner)));
    // Summed absolute deviation from the Table 1 row. The reproduction
    // under-produces the inclined and retrograde families (~4% and ~1%
    // against 13% and 3.5%), which in the paper's own parameter study depend
    // on the unstated disc extent and resolution; everything else is within
    // ~25%.
    let dev: f64 = GROUPS
        .iter()
        .enumerate()
        .map(|(k, g)| (census.fraction(*g) - TABLE1_FRACTION_T0[k]).abs())
        .sum();
    assert!(
        dev < 0.45,
        "summed deviation from Table 1 (t = 0) is {dev:.3}"
    );
}

#[test]
fn the_inner_disc_stays_cold_and_the_outer_disc_is_scattered() {
    // Pfalzner et al. (2024): the disc is left undisturbed out to
    // r_d ≈ 0.28 M_p^{-0.32} r_peri ≈ 33 AU.
    let out = flyby();
    let inner: Vec<_> = out
        .bound()
        .into_iter()
        .filter(|(r0, _)| *r0 < 33.0)
        .collect();
    let cold = inner
        .iter()
        .filter(|(_, e)| e.e < 0.15 && e.i * RAD2DEG < 10.0)
        .count();
    assert!(
        cold as f64 > 0.8 * inner.len() as f64,
        "{cold}/{} inner tracers stayed cold",
        inner.len()
    );
    let outer: Vec<_> = out
        .bound()
        .into_iter()
        .filter(|(r0, _)| *r0 > 100.0)
        .collect();
    let excited = outer.iter().filter(|(_, e)| e.e > 0.3).count();
    assert!(
        excited as f64 > 0.8 * outer.len() as f64,
        "{excited}/{} outer tracers excited",
        outer.len()
    );
}

#[test]
fn neptune_clears_the_near_neptune_low_eccentricity_tracers_first() {
    let ev = evolved();
    let at10 = ev.snapshot_at(10.0);
    let fin = ev.final_snapshot();
    let loss10 = ev.region_loss(at10, near_neptune);
    let loss_final = ev.region_loss(fin, near_neptune);
    let low_e = ev.region_loss(fin, |e| near_neptune(e) && e.e < 0.2);
    let high_e = ev.region_loss(fin, |e| near_neptune(e) && e.e > 0.4);
    let low_i = ev.region_loss(fin, |e| near_neptune(e) && e.i * RAD2DEG < 10.0);
    let high_i = ev.region_loss(fin, |e| near_neptune(e) && e.i * RAD2DEG > 20.0);
    eprintln!(
        "30<q<40 loss: {loss10:.3} at 10 Myr (paper {NEAR_NEPTUNE_LOSS_10MYR}), {loss_final:.3} at {:.0} Myr; \
         e<0.2 {low_e:.3} vs e>0.4 {high_e:.3}; i<10 {low_i:.3} vs i>20 {high_i:.3}",
        fin.t_myr()
    );
    assert!(
        (0.05..0.35).contains(&loss10),
        "10 Myr loss {loss10} far from the paper's {NEAR_NEPTUNE_LOSS_10MYR}"
    );
    assert!(loss_final >= loss10, "losses must accumulate");
    assert!(
        loss_final < NEAR_NEPTUNE_LOSS_4_5GYR,
        "20 Myr cannot exceed the 4.5 Gyr loss"
    );
    assert!(
        low_e > high_e,
        "low-e tracers should be depleted faster: {low_e} vs {high_e}"
    );
    assert!(
        low_i > high_i,
        "low-i tracers should be depleted faster: {low_i} vs {high_i}"
    );
}

#[test]
fn detached_and_sedna_like_populations_are_invariant() {
    let ev = evolved();
    let fin = ev.final_snapshot();
    for g in GROUPS {
        eprintln!(
            "{:<11} retention {:.3} (paper 4.5 Gyr {:.3})",
            g.label(),
            ev.retention(fin, g),
            TABLE1_RETENTION_T4_5GYR[g.index()]
        );
    }
    for g in [
        DynamicalGroup::Detached,
        DynamicalGroup::SednaLike,
        DynamicalGroup::Retrograde,
    ] {
        let r = ev.retention(fin, g);
        assert!(r > 0.9, "{} retention {r}", g.label());
        assert!(
            ev.ejected_fraction(fin, in_group(g)) < 0.05,
            "{} ejected",
            g.label()
        );
    }
    // The cold belt is the family Neptune erodes fastest.
    let cold = ev.retention(fin, DynamicalGroup::ColdKb);
    let hot = ev.retention(fin, DynamicalGroup::HotKb);
    assert!(cold < 1.0 && cold < ev.retention(fin, DynamicalGroup::Detached));
    assert!(cold <= hot + 0.05, "cold {cold} vs hot {hot}");
    // Injected tracers are on their way out (the paper's 99% is at
    // 4.5 Gyr; here the inner giants are only a ring, so the loss runs
    // slower than with Uranus and Saturn scattering directly).
    let inner = ev.ejected_fraction(fin, in_group(DynamicalGroup::Inner));
    eprintln!(
        "inner group ejected by {:.0} Myr: {inner:.3} (paper 4.5 Gyr {INJECTED_EJECTED})",
        fin.t_myr()
    );
    assert!(inner > 0.15, "inner ejection {inner}");
    assert!(
        inner > ev.ejected_fraction(fin, in_group(DynamicalGroup::HotKb)),
        "the planetary region must empty faster than the belt"
    );
}

#[test]
fn very_red_bodies_are_scarce_at_high_inclination_and_eccentricity() {
    let out = flyby();
    let bound = out.bound();
    let (lo_i, hi_i, lo_e, hi_e) = colour_dichotomy(bound.iter().map(|(r, e)| (*r, e)));
    eprintln!("post-flyby very-red fractions: i<21 {lo_i:?} i>21 {hi_i:?}; e<0.42 {lo_e:?} e>0.42 {hi_e:?}");
    assert!(lo_i.unwrap() > hi_i.unwrap());
    assert!(lo_e.unwrap() > hi_e.unwrap());

    let ev = evolved();
    let fin = ev.final_snapshot();
    let survivors = ev.survivors(fin);
    let (lo_i, hi_i, lo_e, hi_e) = colour_dichotomy(survivors.iter().map(|(r, e)| (*r, e)));
    eprintln!(
        "evolved very-red fractions: i<21 {lo_i:?} i>21 {hi_i:?}; e<0.42 {lo_e:?} e>0.42 {hi_e:?}"
    );
    assert!(lo_i.unwrap() > hi_i.unwrap());
    assert!(lo_e.unwrap() > hi_e.unwrap());
}

#[test]
#[ignore = "paper scale: every bound tracer for 4.56 Gyr (hours)"]
fn paper_scale_table_1() {
    let out = run_flyby(&FlybyConfig::model_a1(50_000));
    let ev = evolve(&out.bound(), &EvolutionConfig::paper_scale());
    for t in [0.0, 100.0, 1000.0, 4560.0] {
        report(&ev.census(ev.snapshot_at(t)), &format!("t = {t} Myr"));
    }
    let fin = ev.final_snapshot();
    for g in GROUPS {
        eprintln!(
            "{:<11} retention {:.3} (paper {:.3})",
            g.label(),
            ev.retention(fin, g),
            TABLE1_RETENTION_T4_5GYR[g.index()]
        );
    }
    eprintln!(
        "30<q<40 loss {:.3} (paper {NEAR_NEPTUNE_LOSS_4_5GYR})",
        ev.region_loss(fin, near_neptune)
    );
}
