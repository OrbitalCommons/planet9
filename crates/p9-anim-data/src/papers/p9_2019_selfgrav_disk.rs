//! Film export for `p9-2019-selfgrav-disk`: the numbers its scene and ledger entry draw.

use p9_2016_inclination_instability::suppression::giant_planet_precession_rate_typed;
use p9_2019_selfgrav_disk::disk::reference::{
    M_DISK_EARTH_FIDUCIAL, PRECESSION_PERIOD_HI_YR, PRECESSION_PERIOD_LO_YR,
};
use p9_2019_selfgrav_disk::disk::{DiskProfile, surface_density};
use p9_2019_selfgrav_disk::secular::{ApsidalMode, apsidal_mode, solve};
use p9_core::constants::{EARTH_MASS_SOLAR, RAD2DEG, TWO_PI, YEAR_DAYS};
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::units::{days, radians};
use serde_json::{Value, json};

/// Test-particle semi-major axis for the headline numbers (AU).
const A_TEST_AU: f64 = 250.0;
/// Eccentricity of the test orbit used for the giant-planet comparison.
const E_TEST: f64 = 0.7;
/// Age of the Solar System (yr).
const AGE_YR: f64 = 4.5e9;
/// Each `solve` sums 200 rings of numerically integrated Laplace coefficients
/// (about 10 ms), so the curves are sampled sparsely.
const N_RADII: usize = 21;
/// A lighter disc drawn for comparison (Earth masses).
const LIGHT_DISK_EARTH: f64 = 1.0;
/// Rings drawn in the scene (of the profile's 200).
const N_DRAWN_RINGS: usize = 16;
/// Starting (eccentricity, apse relative to the disc's) of the phase-plane
/// trajectories.
const TRAJECTORY_STARTS: [(f64, f64); 7] = [
    (0.10, 0.0),
    (0.20, 0.0),
    (0.30, 0.0),
    (0.45, 0.0),
    (0.60, 0.0),
    (0.75, 0.0),
    (0.90, 0.0),
];

fn log_grid(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|k| lo * (hi / lo).powf(k as f64 / (n - 1) as f64))
        .collect()
}

/// Apsidal precession period (yr) the four giant planets force on (a, e).
fn planet_period_yr(a: f64, e: f64) -> f64 {
    let rate_per_day = (giant_planet_precession_rate_typed(a, e) * days(1.0) / radians(1.0)).value;
    TWO_PI / rate_per_day / YEAR_DAYS
}

pub fn export() -> Value {
    let disk = DiskProfile::fiducial(M_DISK_EARTH_FIDUCIAL);
    let sol = solve(A_TEST_AU, &disk);

    let rings = disk.rings();
    let step = rings.len() / N_DRAWN_RINGS;
    let drawn: Vec<Value> = rings
        .iter()
        .skip(step / 2)
        .step_by(step)
        .map(|r| {
            json!({
                "a_au": r.a,
                "e": r.e,
                "varpi_deg": r.varpi * RAD2DEG,
                "sigma_earth_per_au2": surface_density(&disk, r.a) / EARTH_MASS_SOLAR,
            })
        })
        .collect();

    // Phase-plane trajectories z(t) = z_f + (z_0 - z_f) exp(iAt), in the frame
    // of the disc apse: k = e cos(dvarpi), h = e sin(dvarpi).
    let phi_f = sol.forced_delta_varpi();
    let (zf_k, zf_h) = (sol.e_forced * phi_f.cos(), sol.e_forced * phi_f.sin());
    let trajectories: Vec<Value> = TRAJECTORY_STARTS
        .iter()
        .map(|&(e0, dphi0)| {
            let (dk, dh) = (e0 * dphi0.cos() - zf_k, e0 * dphi0.sin() - zf_h);
            let (k, h): (Vec<f64>, Vec<f64>) = (0..=96)
                .map(|s| {
                    let (sin, cos) = (TWO_PI * s as f64 / 96.0).sin_cos();
                    (zf_k + dk * cos - dh * sin, zf_h + dk * sin + dh * cos)
                })
                .unzip();
            let librates = apsidal_mode(A_TEST_AU, &disk, e0, dphi0) == ApsidalMode::Libration;
            json!({"e0": e0, "dvarpi0_deg": dphi0 * RAD2DEG, "k": k, "h": h, "librates": librates})
        })
        .collect();

    let radii = log_grid(60.0, 700.0, N_RADII);
    let light = DiskProfile::fiducial(LIGHT_DISK_EARTH);
    let light_period: Vec<f64> = radii
        .iter()
        .map(|&a| solve(a, &light).precession_period_yr())
        .collect();

    // Where the disc starts to turn an orbit faster than the planets do.
    let (mut lo, mut hi) = (60.0_f64, 700.0_f64);
    for _ in 0..18 {
        let mid = (lo * hi).sqrt();
        if solve(mid, &disk).precession_period_yr() < planet_period_yr(mid, E_TEST) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let crossover_au = (lo * hi).sqrt();

    let (disk_period, e_forced): (Vec<f64>, Vec<f64>) = radii
        .iter()
        .map(|&a| {
            let s = solve(a, &disk);
            (s.precession_period_yr(), s.e_forced)
        })
        .unzip();
    let planets_period: Vec<f64> = radii.iter().map(|&a| planet_period_yr(a, E_TEST)).collect();

    json!({
        "disk": {
            "mass_earth": disk.mass_earth(),
            "a_in_au": disk.a_in,
            "a_out_au": disk.a_out,
            "e": disk.e_disk,
            "sigma_index": disk.sigma_index,
            "softening_frac": disk.softening_frac,
            "rings": drawn,
        },
        "a_test_au": A_TEST_AU,
        "e_test": E_TEST,
        "age_yr": AGE_YR,
        "period_yr": sol.precession_period_yr(),
        "period_myr": sol.precession_period_yr() / 1e6,
        "precession_prograde": sol.a_coeff > 0.0,
        "e_forced": sol.e_forced,
        // The libration circle about the forced point that just touches e = 0.
        "e_libration_max": 2.0 * sol.e_forced,
        "etno_eccentricities": BROWN_2017_SAMPLE.iter().map(|o| o.e).collect::<Vec<_>>(),
        "etno_e_min": BROWN_2017_SAMPLE.iter().map(|o| o.e).fold(f64::INFINITY, f64::min),
        "n_etno_librating": BROWN_2017_SAMPLE
            .iter()
            .filter(|o| o.e < 2.0 * sol.e_forced)
            .count(),
        "forced_dvarpi_deg": phi_f * RAD2DEG,
        "paper_forced_dvarpi_deg": 180.0,
        "paper_period_range_yr": [PRECESSION_PERIOD_LO_YR, PRECESSION_PERIOD_HI_YR],
        "planets_period_yr": planet_period_yr(A_TEST_AU, E_TEST),
        "trajectories": trajectories,
        "crossover_au": crossover_au,
        "light_disk_earth": LIGHT_DISK_EARTH,
        "light_disk_period_yr": light_period,
        "radii_au": radii,
        "disk_period_yr": disk_period,
        "planets_period_vs_a_yr": planets_period,
        "e_forced_vs_a": e_forced,
    })
}
