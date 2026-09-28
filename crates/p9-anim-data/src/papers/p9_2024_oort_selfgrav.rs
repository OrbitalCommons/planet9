//! Film export for `p9-2024-oort-selfgrav`: the numbers its scene and ledger entry draw.

use p9_2024_oort_selfgrav::hamiltonian::{HamiltonianParams, j2_fraction};
use p9_2024_oort_selfgrav::miyamoto_nagai::MiyamotoNagaiParams;
use p9_2024_oort_selfgrav::vzlk::{
    VzlkConfig, evolutionary_timescale_gyr, integrate_vzlk, minimum_perihelion,
};
use p9_core::constants::{DEG2RAD, EARTH_MASS_SOLAR, GYR_DAYS};
use serde_json::{Value, json};

/// Age of the solar system (Gyr).
const AGE_GYR: f64 = 4.5;
/// Steps per integrated trajectory.
const N_STEPS: usize = 600;
/// Inclination of the orbit the timescale curves are drawn for.
const I_CURVE_DEG: f64 = 30.0;
/// Perihelion of the orbit the timescale curves are drawn for (AU).
const Q_CURVE_AU: f64 = 100.0;

fn params(n_quadrature: usize, m_ioc_earth: f64) -> HamiltonianParams {
    HamiltonianParams {
        n_quadrature,
        mn_params: MiyamotoNagaiParams {
            m_ioc_solar: m_ioc_earth * EARTH_MASS_SOLAR,
            ..MiyamotoNagaiParams::default_paper()
        },
        ..HamiltonianParams::default_paper()
    }
}

/// One secular trajectory at the paper's (a, J_z), followed for `span` of its
/// own estimated periods.
fn trajectory(
    cfg: &VzlkConfig,
    q0: f64,
    omega0_deg: f64,
    span: f64,
    p: &HamiltonianParams,
) -> Value {
    let e0 = 1.0 - q0 / cfg.a;
    let eta0 = (1.0 - e0 * e0).sqrt();
    let i0 = (cfg.j_z / eta0).clamp(-1.0, 1.0).acos();
    let omega0 = omega0_deg * DEG2RAD;
    let tau = evolutionary_timescale_gyr(cfg.a, e0, i0, omega0, p);
    let dt = span * tau * GYR_DAYS / N_STEPS as f64;
    let track = integrate_vzlk(cfg.a, cfg.j_z, e0, omega0, p, dt, N_STEPS);
    let q_at_age = track
        .iter()
        .take_while(|s| s.t_days / GYR_DAYS <= AGE_GYR)
        .last()
        .map(|s| s.q);
    json!({
        "q0_au": q0,
        "omega0_deg": omega0_deg,
        "i0_deg": i0 / DEG2RAD,
        "timescale_gyr": tau,
        "t_gyr": track.iter().map(|s| s.t_days / GYR_DAYS).collect::<Vec<_>>(),
        "omega_deg": track.iter().map(|s| s.omega / DEG2RAD).collect::<Vec<_>>(),
        "q_au": track.iter().map(|s| s.q).collect::<Vec<_>>(),
        "i_deg": track.iter().map(|s| s.i / DEG2RAD).collect::<Vec<_>>(),
        "q_min_au": track.iter().map(|s| s.q).fold(f64::INFINITY, f64::min),
        "q_max_au": track.iter().map(|s| s.q).fold(0.0, f64::max),
        "q_at_age_au": q_at_age,
    })
}

pub fn export() -> Value {
    let nominal = params(64, 3.0);
    let cfg = VzlkConfig::default_paper();

    let tracks: Vec<Value> = [
        (100.0, 90.0),
        (150.0, 90.0),
        (220.0, 90.0),
        (300.0, 90.0),
        (400.0, 90.0),
        (150.0, 0.0),
        (300.0, 0.0),
    ]
    .iter()
    .map(|&(q0, w0)| trajectory(&cfg, q0, w0, 2.0, &nominal))
    .collect();
    let headline = trajectory(&cfg, 300.0, 90.0, 2.0, &nominal);

    // Timescale against semi-major axis for three cloud masses.
    let a_grid: Vec<f64> = (0..=25).map(|k| 500.0 + 100.0 * k as f64).collect();
    let masses = [1.0, 3.0, 10.0];
    let curves: Vec<Value> = masses
        .iter()
        .map(|&m| {
            let p = params(64, m);
            let tau: Vec<f64> = a_grid
                .iter()
                .map(|&a| {
                    evolutionary_timescale_gyr(
                        a,
                        1.0 - Q_CURVE_AU / a,
                        I_CURVE_DEG * DEG2RAD,
                        0.0,
                        &p,
                    )
                })
                .collect();
            json!({"m_ioc_earth": m, "timescale_gyr": tau})
        })
        .collect();

    // Planetary share of the Hamiltonian against the paper's quoted values.
    let full = HamiltonianParams::default_paper();
    let share: Vec<Value> = [(75.0, 0.02), (50.0, 0.07), (30.0, 0.25)]
        .iter()
        .map(|&(q, published)| {
            let f = j2_fraction(1000.0, 1.0 - q / 1000.0, 30.0 * DEG2RAD, 0.0, &full);
            json!({"q_au": q, "computed": f, "published": published})
        })
        .collect();

    let timescale = evolutionary_timescale_gyr(1000.0, 0.9, 30.0 * DEG2RAD, 0.0, &nominal);

    json!({
        "cloud": {
            "m_ioc_earth": 3.0,
            "a_tilde_au": nominal.mn_params.a_tilde,
            "b_tilde_au": nominal.mn_params.b_tilde,
        },
        "a_au": cfg.a,
        "j_z": cfg.j_z,
        "age_gyr": AGE_GYR,
        "timescale_gyr": timescale,
        "timescale_over_age": timescale / AGE_GYR,
        "tracks": tracks,
        "headline_track": headline,
        "q_min_headline_au": minimum_perihelion(cfg.a, cfg.j_z, 1.0 - 300.0 / cfg.a, 90.0 * DEG2RAD, &nominal),
        "timescale_curves": {
            "a_au": a_grid,
            "q_au": Q_CURVE_AU,
            "i_deg": I_CURVE_DEG,
            "curves": curves,
        },
        "j2_share": share,
        "j2_share_q75": share[0]["computed"],
    })
}
