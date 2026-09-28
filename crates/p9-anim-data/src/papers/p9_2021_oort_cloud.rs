//! Film export for `p9-2021-oort-cloud`: the numbers its scene and ledger entry draw.

use p9_2021_oort_cloud::injection_simulation::{
    InjectionConfig, evolve_population, f_varpi_anti_aligned,
};
use p9_2021_oort_cloud::oort_cloud::{
    OortCloudConfig, generate_ioc_population, generate_scattered_disk,
};
use p9_core::analysis::circular::wrap_to_pi;
use p9_core::constants::GYR_DAYS;
use rand::SeedableRng;
use serde_json::{Value, json};

/// Seeds pooled for the film (the crate's seed-averaged regression set).
const SEEDS: [u64; 5] = [42, 43, 44, 45, 46];
/// Particles per seed, and the secular time-compression factor, as in the
/// crate's reduced-scale regression test.
const N_IOC: usize = 40;
const N_SCATTERED: usize = 48;
const MASS_BOOST: f64 = 300.0;

/// Published confinement fractions (Batygin & Brown 2021).
const PUBLISHED_F_IOC: f64 = 0.67;
const PUBLISHED_F_SCATTERED: f64 = 0.88;

pub fn export() -> Value {
    let config = InjectionConfig {
        ioc_config: OortCloudConfig {
            n_particles: N_IOC,
            a_min: 800.0,
            a_max: 2500.0,
            q_min: 60.0,
            q_max: 300.0,
            ..OortCloudConfig::nominal()
        },
        n_scattered: N_SCATTERED,
        mass_boost: MASS_BOOST,
        ..InjectionConfig::nominal()
    };
    let p9_varpi = config.p9.omega + config.p9.omega_big;
    let duration = config.duration_gyr * GYR_DAYS / config.mass_boost;
    let q_cut = config.q_injection_threshold;

    let mut particles = Vec::new();
    let mut injected_dvarpi = Vec::new();
    let mut scattered_dvarpi = Vec::new();
    let mut control_dvarpi = Vec::new();
    let mut n_injectable = 0usize;
    let mut n_control_crossed = 0usize;

    // One seed: the inner Oort cloud with Planet Nine, the same cloud without
    // it, and a scattered disk with it (shorter run, finer step, as in the
    // crate's `simulate_injection`).
    let one_seed = |seed: u64| {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let ioc = generate_ioc_population(&config.ioc_config, &mut rng);
        let disk = generate_scattered_disk(config.n_scattered, &mut rng);
        let (boost, dt) = (config.mass_boost, config.dt_days);
        (
            evolve_population(&ioc, Some(&config.p9), boost, duration, dt),
            evolve_population(&ioc, None, boost, duration, dt),
            evolve_population(&disk, Some(&config.p9), boost, duration / 10.0, 7.0e4),
        )
    };
    let runs: Vec<_> = std::thread::scope(|s| {
        let handles: Vec<_> = SEEDS
            .iter()
            .map(|&seed| s.spawn(move || one_seed(seed)))
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    for (with_p9, control, disk_run) in &runs {
        for (o, c) in with_p9.iter().zip(control) {
            let q0 = o.initial.a * (1.0 - o.initial.e);
            let injectable = q0 > q_cut && o.initial.a > config.a_dkb_min;
            let injected = injectable && o.min_q < q_cut;
            n_injectable += injectable as usize;
            n_control_crossed += (injectable && c.min_q < q_cut) as usize;
            let dvarpi = o
                .final_elements
                .as_ref()
                .map(|f| wrap_to_pi(f.longitude_of_perihelion() - p9_varpi));
            if let (true, Some(d)) = (injected, dvarpi) {
                injected_dvarpi.push(d);
            }
            particles.push(json!({
                "a_au": o.initial.a,
                "q0_au": q0,
                "q_min_au": o.min_q,
                "q_min_control_au": c.min_q,
                "injectable": injectable,
                "injected": injected,
                "survived": o.final_elements.is_some(),
            }));
        }
        control_dvarpi.extend(
            control
                .iter()
                .filter_map(|o| o.final_elements.as_ref())
                .map(|f| wrap_to_pi(f.longitude_of_perihelion() - p9_varpi)),
        );
        scattered_dvarpi.extend(
            disk_run
                .iter()
                .filter_map(|o| o.final_elements.as_ref())
                .map(|f| wrap_to_pi(f.longitude_of_perihelion() - p9_varpi)),
        );
    }

    let n_injected = particles
        .iter()
        .filter(|p| p["injected"] == json!(true))
        .count();
    let degrees = |v: &[f64]| v.iter().map(|x| x.to_degrees()).collect::<Vec<_>>();

    json!({
        "p9": {
            "mass_earth": config.p9.mass_earth,
            "a_au": config.p9.a,
            "e": config.p9.e,
            "i_deg": config.p9.i.to_degrees(),
        },
        "q_injection_au": q_cut,
        "equivalent_gyr": config.duration_gyr,
        "mass_boost": config.mass_boost,
        "n_seeds": SEEDS.len(),
        "particles": particles,
        "n_injectable": n_injectable,
        "n_injected": n_injected,
        "n_control_crossed": n_control_crossed,
        "injection_fraction": n_injected as f64 / n_injectable.max(1) as f64,
        "f_varpi_ioc": f_varpi_anti_aligned(&injected_dvarpi),
        "f_varpi_scattered": f_varpi_anti_aligned(&scattered_dvarpi),
        "f_varpi_control": f_varpi_anti_aligned(&control_dvarpi),
        "dvarpi_ioc_deg": degrees(&injected_dvarpi),
        "dvarpi_scattered_deg": degrees(&scattered_dvarpi),
        "dvarpi_control_deg": degrees(&control_dvarpi),
        "published": {"f_varpi_ioc": PUBLISHED_F_IOC, "f_varpi_scattered": PUBLISHED_F_SCATTERED},
    })
}
