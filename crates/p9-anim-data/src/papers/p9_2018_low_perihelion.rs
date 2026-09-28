//! Film export for `p9-2018-low-perihelion`: the numbers its scene and ledger entry draw.
//!
//! The crate measures how tightly a single eccentric Planet Nine confines the
//! apse of a test ETNO (the half-width of the Δϖ libration island, and the
//! depth of the apsidal well) as P9's perihelion q9 = a9(1 − e9) is lowered at
//! fixed a9. The scene sweeps q9 from the canonical 280 AU down to 70 AU at
//! a9 = 700 AU, and repeats the sweep at the paper's other planet size,
//! a9 = 1500 AU, over the paper's q9 range.

use p9_2018_low_perihelion::confinement::{Confinement, FavoredApse, analyze};
use p9_2018_low_perihelion::reference::{
    BB_A9_AU, BB_E9, BB_Q9_AU, E9_SWEEP_HI, E9_SWEEP_LO, LOW_Q9_MIN_AU, P9_MASS_EARTH,
    TEST_ETNO_A_AU, etno_median_a,
};
use p9_2018_low_perihelion::sweep::{lowest_confining_q9, sweep_q9};
use p9_core::data::ephemeris_constraint::brown_batygin_orbit;
use p9_core::data::etno::BROWN_2017_SAMPLE;
use serde_json::{Value, json};

/// Sweep points along each q9 track.
const N_SWEEP: usize = 13;
/// The paper's wider planet (Cáceres & Gomes run a9 = 700 and 1500 AU).
const WIDE_A9_AU: f64 = 1500.0;
/// q9 range of the paper's models (60–300 AU).
const PAPER_Q9_LO_AU: f64 = 60.0;
const PAPER_Q9_HI_AU: f64 = 300.0;

fn row(e9: f64, q9: f64, c: &Confinement) -> Value {
    json!({
        "e9": e9,
        "q9_au": q9,
        "half_width_deg": c.libration_half_width.to_degrees(),
        "well_depth": c.well_depth,
        "forced_e": c.forced_eccentricity,
        "aligned": c.favored_apse == FavoredApse::Aligned,
        "confined": c.is_confined(),
    })
}

pub fn export() -> Value {
    // Headline sweep: a9 = 700 AU, e9 0.6 → 0.9, q9 280 → 70 AU.
    let rows = sweep_q9(TEST_ETNO_A_AU, BB_A9_AU, P9_MASS_EARTH, N_SWEEP);
    let track_700: Vec<Value> = rows
        .iter()
        .map(|r| row(r.e9, r.q9_au, &r.confinement))
        .collect();
    let first = rows.first().expect("sweep has rows").confinement;
    let last = rows.last().expect("sweep has rows").confinement;

    // Same measurement for the paper's a9 = 1500 AU planet over q9 = 60–300 AU.
    let track_1500: Vec<Value> = (0..N_SWEEP)
        .map(|k| {
            let q9 = PAPER_Q9_HI_AU
                - (PAPER_Q9_HI_AU - PAPER_Q9_LO_AU) * k as f64 / (N_SWEEP - 1) as f64;
            let e9 = 1.0 - q9 / WIDE_A9_AU;
            row(
                e9,
                q9,
                &analyze(TEST_ETNO_A_AU, WIDE_A9_AU, e9, P9_MASS_EARTH),
            )
        })
        .collect();

    let bb = brown_batygin_orbit();
    let etnos: Vec<Value> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| {
            json!({
                "name": o.name,
                "a": o.a,
                "e": o.e,
                "q": o.perihelion(),
                "varpi_deg": o.longitude_of_perihelion().to_degrees(),
            })
        })
        .collect();

    json!({
        "a9_au": BB_A9_AU,
        "e9_canonical": BB_E9,
        "q9_canonical_au": BB_Q9_AU,
        "e9_range": [E9_SWEEP_LO, E9_SWEEP_HI],
        "q9_floor_au": LOW_Q9_MIN_AU,
        "p9_varpi_deg": (bb.omega + bb.omega_big).to_degrees().rem_euclid(360.0),
        "mass_earth": P9_MASS_EARTH,
        "test_a_au": TEST_ETNO_A_AU,
        "etno_median_a_au": etno_median_a(),
        "wide_a9_au": WIDE_A9_AU,
        "track_700": track_700,
        "track_1500": track_1500,
        "half_width_canonical_deg": first.libration_half_width.to_degrees(),
        "half_width_low_deg": last.libration_half_width.to_degrees(),
        "narrowing": 1.0 - last.libration_half_width / first.libration_half_width,
        "lowest_confining_q9_au": lowest_confining_q9(&rows),
        "etnos": etnos,
    })
}
