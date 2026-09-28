//! Film export for `p9-2018-chaotic-dynamics`: the numbers its scene and ledger entry draw.

use p9_2018_chaotic_dynamics::chaos::{
    ChaosMap, chirikov_overlap_parameter_general, critical_eccentricity, is_chaotic,
    j_one_location_typed, overlap_zone_width_typed,
};
use p9_2018_chaotic_dynamics::reference::{FIDUCIAL_A9_AU, FIDUCIAL_E9, FIDUCIAL_M9_EARTH};
use p9_2018_chaotic_dynamics::width::libration_width_typed;
use p9_core::analysis::resonance::resonance_semi_major_axis;
use p9_core::constants::{A_NEPTUNE_AU, MASS_NEPTUNE_SOLAR};
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::units::{Length, au};
use serde_json::{Value, json};

/// Planet Nine masses swept for the chaos boundary (Earth masses).
const MASSES: [f64; 3] = [5.0, 10.0, 20.0];

/// Semi-major-axis window of the chaos map (AU).
const A_LO: f64 = 150.0;
const A_HI: f64 = 495.0;

fn in_au(l: Length) -> f64 {
    (l / au(1.0)).value
}

/// Neptune's resonance-overlap parameter at (a, e): the crate's general
/// Chirikov chain form fed Neptune's orbit and mass.
fn neptune_k(a: f64, e: f64) -> f64 {
    chirikov_overlap_parameter_general(a, a * (1.0 - e), A_NEPTUNE_AU, MASS_NEPTUNE_SOLAR)
}

pub fn export() -> Value {
    let (a9, m9, e9) = (FIDUCIAL_A9_AU, FIDUCIAL_M9_EARTH, FIDUCIAL_E9);

    // Eccentricity grid shared by the resonance widths and the overlap zone.
    let e_grid: Vec<f64> = (0..=19).map(|k| 0.05 * k as f64).collect();

    // The strong interior j:1 resonances and their libration widths against e.
    let j_one: Vec<Value> = (2..=5i64)
        .map(|j| {
            let widths: Vec<f64> = e_grid
                .iter()
                .map(|&e| in_au(libration_width_typed(j, e, m9, a9, e9)))
                .collect();
            json!({"j": j, "a": in_au(j_one_location_typed(j, a9)), "width": widths})
        })
        .collect();

    // First-order N:(N+1) resonances, which crowd together toward the planet.
    let first_order: Vec<Value> = (1..=40u32)
        .map(|n| json!({"n": n, "a": resonance_semi_major_axis(n, n + 1, a9)}))
        .collect();

    // Inward extent of the overlapped (chaotic) zone next to Planet Nine.
    let zone: Vec<f64> = e_grid
        .iter()
        .map(|&e| in_au(overlap_zone_width_typed(e, a9, m9)))
        .collect();

    // K = 1 boundary of the Planet Nine chaotic zone for each mass.
    let a_line: Vec<f64> = (0..=138).map(|k| A_LO + 2.5 * k as f64).collect();
    let boundaries: Vec<Value> = MASSES
        .iter()
        .map(|&m| {
            let e_crit: Vec<Option<f64>> = a_line
                .iter()
                .map(|&a| critical_eccentricity(a, a9, m, e9))
                .collect();
            let map = ChaosMap::build(A_LO, A_HI, 69, 0.0, 0.95, 38, a9, m, e9);
            json!({
                "mass_earth": m,
                "e_crit": e_crit,
                "circular_zone_au": in_au(overlap_zone_width_typed(0.0, a9, m)),
                "chaotic_fraction": map.chaotic_fraction(),
            })
        })
        .collect();

    // Neptune's own overlap boundary: the eccentricity above which the
    // perihelion dips into Neptune's resonance web (K_N = 1).
    let neptune_e: Vec<Option<f64>> = a_line
        .iter()
        .map(|&a| {
            let grid: Vec<f64> = (0..=950).map(|k| 0.001 * k as f64).collect();
            grid.into_iter().find(|&e| neptune_k(a, e) > 1.0)
        })
        .collect();

    // The fiducial map, each cell tagged 0 regular, 1 Planet Nine overlap,
    // 2 Neptune overlap.
    let map = ChaosMap::build(A_LO, A_HI, 69, 0.0, 0.95, 38, a9, m9, e9);
    let cells: Vec<Vec<u8>> = map
        .a_vals
        .iter()
        .zip(&map.k)
        .map(|(&a, col)| {
            map.e_vals
                .iter()
                .zip(col)
                .map(|(&e, &k)| {
                    if k > 1.0 {
                        1
                    } else if neptune_k(a, e) > 1.0 {
                        2
                    } else {
                        0
                    }
                })
                .collect()
        })
        .collect();
    let n_cells = (map.a_vals.len() * map.e_vals.len()) as f64;
    let regular_fraction = cells.iter().flatten().filter(|&&c| c == 0).count() as f64 / n_cells;

    // The observed distant objects placed on the map.
    let etnos: Vec<Value> = BROWN_2017_SAMPLE
        .iter()
        .filter(|o| o.a < a9)
        .map(|o| {
            json!({
                "name": o.name,
                "a": o.a,
                "e": o.e,
                "chaotic_p9": is_chaotic(o.a, o.e, a9, m9, e9),
                "chaotic_neptune": neptune_k(o.a, o.e) > 1.0,
            })
        })
        .collect();
    let n_regular = etnos
        .iter()
        .filter(|o| o["chaotic_p9"] == json!(false) && o["chaotic_neptune"] == json!(false))
        .count();
    let n_etno = etnos.len();

    json!({
        "a9": a9,
        "m9_earth": m9,
        "e9": e9,
        "e_grid": e_grid,
        "j_one": j_one,
        "first_order": first_order,
        "overlap_zone_au": zone,
        "zone_e085_au": in_au(overlap_zone_width_typed(0.85, a9, m9)),
        "a_line": a_line,
        "boundaries": boundaries,
        "neptune_e": neptune_e,
        "map": {"a": map.a_vals, "e": map.e_vals, "cells": cells},
        "chaotic_fraction": map.chaotic_fraction(),
        "regular_fraction": regular_fraction,
        "etnos": etnos,
        "n_etno": n_etno,
        "n_etno_regular": n_regular,
        "etno_regular_fraction": n_regular as f64 / n_etno as f64,
    })
}
