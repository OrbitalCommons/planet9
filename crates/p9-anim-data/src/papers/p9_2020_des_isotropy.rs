//! Film export for `p9-2020-des-isotropy`: the numbers its scene and ledger entry draw.

use std::cell::RefCell;

use p9_2020_des_isotropy::analysis::{
    BatteryCell, PAPER_KUIPER_SIGNIFICANT_TESTS, PAPER_KUIPER_TOTAL_TESTS, run_battery,
};
use p9_2020_des_isotropy::des_sample::{Angle, DES_ETNOS, SampleCase, case_sample};
use p9_2020_des_isotropy::selection::{
    NullModel, ecliptic_to_equatorial, perihelion_direction, selection_aware_p,
};
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::analysis::surveys::{
    DES_FOOTPRINT_BANDS, des_footprint_contains, des_footprint_solid_angle_deg2,
};
use p9_core::constants::{DEG2RAD, RAD2DEG};
use serde_json::{Value, json};

const SEED: u64 = 2020;
const MC_ITERS: usize = 8000;
const ALPHA: f64 = 0.05;
/// Width of the longitude-of-perihelion bins of the null histogram (degrees).
const VARPI_BIN_DEG: f64 = 30.0;

fn cells_json(cells: &[BatteryCell]) -> Vec<Value> {
    cells
        .iter()
        .map(|c| {
            json!({
                "case": c.case,
                "angle": c.angle,
                "n": c.p.n,
                "r_bar": c.p.r_bar,
                "kuiper_p": c.p.kuiper_p,
                "rayleigh_p": c.p.rayleigh_p,
                "watson_p": c.p.watson_p,
                "beran_p": c.p.beran_p,
            })
        })
        .collect()
}

fn kuiper_significant(cells: &[BatteryCell]) -> usize {
    cells.iter().filter(|c| c.p.kuiper_p < ALPHA).count()
}

/// Every longitude of perihelion the crate's null draws for `case`, binned.
/// The statistic closure is handed each synthetic sample in turn, so it doubles
/// as a tap on the null population.
fn null_varpi_histogram(case: SampleCase, null: NullModel) -> Vec<f64> {
    let n_bins = (360.0 / VARPI_BIN_DEG) as usize;
    let counts = RefCell::new(vec![0usize; n_bins]);
    let calls = RefCell::new(0usize);
    selection_aware_p(
        &case_sample(case),
        Angle::Varpi,
        null,
        SEED,
        MC_ITERS,
        |angles| {
            let mut k = calls.borrow_mut();
            if *k > 0 {
                let mut c = counts.borrow_mut();
                for a in angles {
                    let bin = ((a * RAD2DEG).rem_euclid(360.0) / VARPI_BIN_DEG) as usize;
                    c[bin.min(n_bins - 1)] += 1;
                }
            }
            *k += 1;
            0.0
        },
    );
    let counts = counts.into_inner();
    let total: usize = counts.iter().sum();
    counts
        .iter()
        .map(|&c| c as f64 / total as f64 / VARPI_BIN_DEG)
        .collect()
}

pub fn export() -> Value {
    let objects: Vec<Value> = DES_ETNOS
        .iter()
        .map(|o| {
            let (lambda, beta) =
                perihelion_direction(o.i_deg * DEG2RAD, Angle::ArgPeri.of(o), Angle::Node.of(o));
            let (ra, dec) = ecliptic_to_equatorial(lambda, beta);
            let (ra_deg, dec_deg) = (ra * RAD2DEG, dec * RAD2DEG);
            json!({
                "name": o.name,
                "a": o.a,
                "e": o.e,
                "q": o.a * (1.0 - o.e),
                "i_deg": o.i_deg,
                "node_deg": Angle::Node.of(o) * RAD2DEG,
                "argp_deg": Angle::ArgPeri.of(o) * RAD2DEG,
                "varpi_deg": Angle::Varpi.of(o) * RAD2DEG,
                "y4_discovery": o.y4_discovery,
                "ra_deg": ra_deg,
                "dec_deg": dec_deg,
                "in_footprint": des_footprint_contains(ra_deg, dec_deg),
            })
        })
        .collect();

    let footprint: Vec<Value> = DES_FOOTPRINT_BANDS
        .iter()
        .map(|b| {
            json!({
                "ra_start_deg": b.ra_start_deg,
                "ra_end_deg": b.ra_end_deg,
                "dec_min_deg": b.dec_min_deg,
                "dec_max_deg": b.dec_max_deg,
            })
        })
        .collect();

    let flat = run_battery(NullModel::Uniform, SEED, MC_ITERS);
    let sel = run_battery(NullModel::DesSelection, SEED, MC_ITERS);

    let cases: Vec<Value> = SampleCase::all()
        .iter()
        .map(|c| {
            let (a_min, q_min) = c.cuts();
            json!({
                "case": format!("{c:?}"),
                "a_min": a_min,
                "q_min": q_min,
                "n": case_sample(*c).len(),
                "names": case_sample(*c).iter().map(|o| o.name).collect::<Vec<_>>(),
            })
        })
        .collect();

    // The Planet Nine angle in the largest sample (Case 1, a > 150, q > 30).
    let varpi_case1 = |cells: &[BatteryCell]| {
        cells
            .iter()
            .find(|c| c.case == "Case1" && c.angle == Angle::Varpi.label())
            .map(|c| c.p.kuiper_p)
            .unwrap()
    };
    let p_flat = varpi_case1(&flat);
    let p_sel = varpi_case1(&sel);
    // The paper's verdict rests on the strongest of its twelve tests, read
    // with the number of tests in mind: the chance that the smallest of
    // `n` independent p-values is at least this small.
    let min_p_flat = flat.iter().map(|c| c.p.kuiper_p).fold(1.0_f64, f64::min);
    let min_p_sel = sel.iter().map(|c| c.p.kuiper_p).fold(1.0_f64, f64::min);
    let after_trials = |p: f64| 1.0 - (1.0 - p).powi(sel.len() as i32);

    let bin_centres: Vec<f64> = (0..(360.0 / VARPI_BIN_DEG) as usize)
        .map(|k| (k as f64 + 0.5) * VARPI_BIN_DEG)
        .collect();

    json!({
        "objects": objects,
        "footprint": footprint,
        "footprint_deg2": des_footprint_solid_angle_deg2(),
        "cases": cases,
        "battery_flat": cells_json(&flat),
        "battery_selection": cells_json(&sel),
        "alpha": ALPHA,
        "n_tests": sel.len(),
        "n_significant_flat": kuiper_significant(&flat),
        "n_significant_selection": kuiper_significant(&sel),
        "min_p_selection": min_p_sel,
        "paper_n_tests": PAPER_KUIPER_TOTAL_TESTS,
        "paper_n_significant": PAPER_KUIPER_SIGNIFICANT_TESTS,
        "p_varpi_flat": p_flat,
        "p_varpi_selection": p_sel,
        "min_p_flat": min_p_flat,
        "p_after_trials": after_trials(min_p_sel),
        "sigma_flat": p_value_to_sigma(after_trials(min_p_flat)),
        "sigma": p_value_to_sigma(after_trials(min_p_sel)),
        "n_sample": case_sample(SampleCase::Case1).len(),
        "null_varpi": {
            "bin_deg": VARPI_BIN_DEG,
            "centre_deg": bin_centres,
            "flat": null_varpi_histogram(SampleCase::Case1, NullModel::Uniform),
            "selection": null_varpi_histogram(SampleCase::Case1, NullModel::DesSelection),
        },
    })
}
