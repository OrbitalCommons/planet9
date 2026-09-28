//! Film export for `p9-2024-primordial-alignment`: the numbers its scene and ledger entry draw.

use p9_2024_primordial_alignment::convergence::{
    PUBLISHED_CONVERGENCE_GYR, Precessor, observed_sednoids, primordial_convergence_epoch,
    resultant_at_lookback, scan_convergence,
};
use p9_2024_primordial_alignment::precession::apsidal_precession_rate;
use p9_2026_cluster_inclinations::sample::HIGH_Q_SAMPLE;
use p9_core::analysis::circular::circular_mean;
use p9_core::constants::TWO_PI;
use serde_json::{Value, json};

/// The paper's three sednoids, as named in the workspace's JPL-epoch table.
const SEDNOIDS: [&str; 3] = [
    "Sedna",
    "Alicanto (2012 VP113)",
    "Leleakuhonua (2015 TG387)",
];
/// Look-back window and resolution of the scan.
const MAX_GYR: f64 = 5.0;
const N_SCAN: usize = 2500;
/// Samples of each plotted curve.
const N_PLOT: usize = 500;
/// Published primordial longitude of perihelion (deg), Huang & Gladman abstract.
const PUBLISHED_VARPI_DEG: f64 = 200.0;
/// Mean resultant length counted as a tight alignment.
const TIGHT_R_BAR: f64 = 0.95;

fn wrap_deg(rad: f64) -> f64 {
    rad.to_degrees().rem_euclid(360.0)
}

fn curve(set: &[Precessor]) -> Vec<f64> {
    (0..=N_PLOT)
        .map(|k| resultant_at_lookback(set, MAX_GYR * 1e9 * k as f64 / N_PLOT as f64))
        .collect()
}

pub fn export() -> Value {
    let sednoids: Vec<Precessor> = SEDNOIDS
        .iter()
        .map(|&name| {
            let o = HIGH_Q_SAMPLE
                .iter()
                .find(|o| o.name == name)
                .expect("sednoid present in the high-q sample");
            let el = o.elements();
            Precessor {
                name: o.name,
                varpi0: el.longitude_of_perihelion(),
                rate_rad_per_yr: apsidal_precession_rate(el.a, el.e, el.i),
            }
        })
        .collect();

    let tau_gyr: Vec<f64> = (0..=N_PLOT)
        .map(|k| MAX_GYR * k as f64 / N_PLOT as f64)
        .collect();
    let objects: Vec<Value> = SEDNOIDS
        .iter()
        .zip(&sednoids)
        .map(|(&name, p)| {
            let o = HIGH_Q_SAMPLE.iter().find(|o| o.name == name).unwrap();
            let track: Vec<f64> = tau_gyr
                .iter()
                .map(|&t| wrap_deg(p.varpi0 - p.rate_rad_per_yr * t * 1e9))
                .collect();
            json!({
                "name": name,
                "a_au": o.a,
                "e": o.e,
                "q_au": o.perihelion(),
                "i_deg": o.i_deg,
                "varpi_now_deg": wrap_deg(p.varpi0),
                "period_gyr": TWO_PI / p.rate_rad_per_yr / 1e9,
                "turns_in_age": p.rate_rad_per_yr * PUBLISHED_CONVERGENCE_GYR * 1e9 / TWO_PI,
                "varpi_at_birth_deg":
                    wrap_deg(p.varpi0 - p.rate_rad_per_yr * PUBLISHED_CONVERGENCE_GYR * 1e9),
                "varpi_track_deg": track,
            })
        })
        .collect();

    let scan = scan_convergence(&sednoids, MAX_GYR, N_SCAN);
    let r_at_birth = resultant_at_lookback(&sednoids, PUBLISHED_CONVERGENCE_GYR * 1e9);

    let varpi_at = |t_gyr: f64| -> f64 {
        let v: Vec<f64> = sednoids
            .iter()
            .map(|p| p.varpi0 - p.rate_rad_per_yr * t_gyr * 1e9)
            .collect();
        wrap_deg(circular_mean(&v).unwrap_or(0.0))
    };

    // Fraction of the look-back window in which the three apsides are tightly
    // aligned: how special the alignment at the best epoch is.
    let r_curve = curve(&sednoids);
    let tight = r_curve.iter().filter(|&&r| r >= TIGHT_R_BAR).count() as f64 / r_curve.len() as f64;

    let brown = observed_sednoids();
    let brown_scan = scan_convergence(&brown, MAX_GYR, N_SCAN);
    let mechanism = primordial_convergence_epoch(PUBLISHED_CONVERGENCE_GYR, N_SCAN);

    json!({
        "published_epoch_gyr": PUBLISHED_CONVERGENCE_GYR,
        "published_varpi_deg": PUBLISHED_VARPI_DEG,
        "objects": objects,
        "tau_gyr": tau_gyr,
        "r_bar": r_curve,
        "r_bar_now": scan.r_bar_present,
        "r_bar_at_birth": r_at_birth,
        "best_epoch_gyr": scan.t_star_gyr,
        "best_r_bar": scan.r_bar_max,
        "best_varpi_deg": varpi_at(scan.t_star_gyr),
        "tight_r_bar": TIGHT_R_BAR,
        "fraction_tight": tight,
        "brown2017": {
            "n": brown.len(),
            "r_bar": curve(&brown),
            "r_bar_now": brown_scan.r_bar_present,
            "best_epoch_gyr": brown_scan.t_star_gyr,
            "best_r_bar": brown_scan.r_bar_max,
        },
        "mechanism": {
            "recovered_epoch_gyr": mechanism.t_star_gyr,
            "r_bar_max": mechanism.r_bar_max,
            "r_bar_today": mechanism.r_bar_present,
        },
    })
}
