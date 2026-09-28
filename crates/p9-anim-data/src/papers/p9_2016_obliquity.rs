//! Film export for `p9-2016-obliquity`: the numbers its scene and ledger entry draw.

use p9_2016_obliquity::parameter_survey::find_required_inclination;
use p9_2016_obliquity::secular_hamiltonian::{SecularParams, SpinOrbitState, integrate_obliquity};
use p9_core::constants::{EARTH_MASS_SOLAR, GYR_DAYS, YEAR_DAYS};
use p9_core::initial_conditions::giant_planets::{
    giant_planet_angular_momentum, p9_angular_momentum,
};
use serde_json::{Value, json};
use std::f64::consts::PI;
use std::thread;

/// Observed solar obliquity relative to the invariable plane (degrees).
const OBSERVED_OBLIQUITY_DEG: f64 = 6.0;

/// Age of the Solar System (Gyr).
const AGE_GYR: f64 = 4.5;

/// Planet Nine masses surveyed (Earth masses), Bailey et al. Fig. 2.
const MASSES_EARTH: [f64; 3] = [10.0, 15.0, 20.0];

/// Survey grid: semi-major axes (AU) and eccentricities.
const A_GRID_AU: [f64; 6] = [400.0, 500.0, 600.0, 700.0, 800.0, 900.0];
const E_GRID: [f64; 7] = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];

/// Perihelion window the paper surveys (AU).
const PERIHELION_RANGE_AU: (f64, f64) = (150.0, 350.0);

/// How closely the solved inclination must reproduce the observed tilt (degrees).
const SOLVE_TOLERANCE_DEG: f64 = 0.02;

/// Illustrative perturber of Bailey et al. Fig. 4.
const SHOWCASE: (f64, f64, f64) = (15.0, 500.0, 0.5);

/// Batygin & Brown (2016) nominal perturber.
const NOMINAL: (f64, f64, f64, f64) = (10.0, 700.0, 0.6, 30.0);

/// Inclinations whose 4.5 Gyr histories are drawn (degrees).
const HISTORY_INCLINATIONS_DEG: [f64; 3] = [10.0, 20.0, 30.0];

/// Published: inclinations needed across the surveyed grid (degrees), and the
/// tilt the nominal Batygin & Brown perturber reaches (degrees).
const PUBLISHED_REQUIRED_I_DEG: (f64, f64) = (15.0, 30.0);
const PUBLISHED_NOMINAL_TILT_DEG: (f64, f64) = (3.0, 5.0);

fn history(mass_earth: f64, a9: f64, e9: f64, i9_deg: f64) -> Value {
    let m9_solar = mass_earth * EARTH_MASS_SOLAR;
    let initial = SpinOrbitState::from_inclinations(
        i9_deg.to_radians(),
        PI,
        giant_planet_angular_momentum(),
        p9_angular_momentum(m9_solar, a9, e9),
    );
    let params = SecularParams {
        m9_solar,
        a9,
        e9,
        t_total: AGE_GYR * GYR_DAYS,
        dt: 5e4 * YEAR_DAYS,
    };
    let snaps = integrate_obliquity(initial, &params, 0.05 * GYR_DAYS);
    let last = snaps.last().expect("the integration returns snapshots");
    json!({
        "mass_earth": mass_earth,
        "a_au": a9,
        "e": e9,
        "i9_deg": i9_deg,
        "t_gyr": snaps.iter().map(|s| s.t / GYR_DAYS).collect::<Vec<_>>(),
        "obliquity_deg": snaps.iter().map(|s| s.obliquity.to_degrees()).collect::<Vec<_>>(),
        "final_obliquity_deg": last.obliquity.to_degrees(),
    })
}

/// One survey cell: the inclination that yields the observed tilt, if any.
struct Cell {
    mass_earth: f64,
    a_au: f64,
    e: f64,
    required_i_deg: Option<f64>,
}

fn survey() -> Vec<Cell> {
    let mut grid = Vec::new();
    for &mass_earth in &MASSES_EARTH {
        for &a_au in &A_GRID_AU {
            for &e in &E_GRID {
                let q = a_au * (1.0 - e);
                if q > PERIHELION_RANGE_AU.0 && q < PERIHELION_RANGE_AU.1 {
                    grid.push((mass_earth, a_au, e));
                }
            }
        }
    }
    thread::scope(|scope| {
        let handles: Vec<_> = grid
            .iter()
            .map(|&(mass_earth, a_au, e)| {
                scope.spawn(move || Cell {
                    mass_earth,
                    a_au,
                    e,
                    required_i_deg: find_required_inclination(
                        mass_earth,
                        a_au,
                        e,
                        OBSERVED_OBLIQUITY_DEG,
                        SOLVE_TOLERANCE_DEG,
                    )
                    .map(|(i, _)| i.to_degrees()),
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("survey cell"))
            .collect()
    })
}

pub fn export() -> Value {
    let cells = survey();
    let solved: Vec<f64> = cells.iter().filter_map(|c| c.required_i_deg).collect();
    let mut sorted = solved.clone();
    sorted.sort_by(f64::total_cmp);

    let surveys: Vec<Value> = MASSES_EARTH
        .iter()
        .map(|&m| {
            let rows: Vec<Value> = cells
                .iter()
                .filter(|c| c.mass_earth == m)
                .map(|c| {
                    json!({
                        "a_au": c.a_au,
                        "e": c.e,
                        "perihelion_au": c.a_au * (1.0 - c.e),
                        "required_i_deg": c.required_i_deg,
                    })
                })
                .collect();
            json!({"mass_earth": m, "cells": rows})
        })
        .collect();

    let showcase: Vec<Value> = HISTORY_INCLINATIONS_DEG
        .iter()
        .map(|&i| history(SHOWCASE.0, SHOWCASE.1, SHOWCASE.2, i))
        .collect();
    let showcase_required = cells
        .iter()
        .find(|c| c.mass_earth == SHOWCASE.0 && c.a_au == SHOWCASE.1 && c.e == SHOWCASE.2)
        .and_then(|c| c.required_i_deg);
    let nominal = history(NOMINAL.0, NOMINAL.1, NOMINAL.2, NOMINAL.3);

    json!({
        "observed_obliquity_deg": OBSERVED_OBLIQUITY_DEG,
        "age_gyr": AGE_GYR,
        "a_grid_au": A_GRID_AU,
        "e_grid": E_GRID,
        "surveys": surveys,
        "cells_total": cells.len(),
        "cells_solved": solved.len(),
        "required_i_min_deg": sorted.first(),
        "required_i_max_deg": sorted.last(),
        "required_i_median_deg": sorted.get(sorted.len() / 2),
        "showcase": {
            "mass_earth": SHOWCASE.0,
            "a_au": SHOWCASE.1,
            "e": SHOWCASE.2,
            "required_i_deg": showcase_required,
            "histories": showcase,
        },
        "nominal": nominal,
        "nominal_obliquity_deg": nominal["final_obliquity_deg"],
        "published": {
            "required_i_deg": PUBLISHED_REQUIRED_I_DEG,
            "nominal_tilt_deg": PUBLISHED_NOMINAL_TILT_DEG,
        },
    })
}
