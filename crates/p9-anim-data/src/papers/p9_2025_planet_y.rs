//! Film export for `p9-2025-planet-y`: the numbers its scene and ledger entry draw.

use p9_2025_planet_y::laplace_plane::published::{
    A_Y_HIGH_AU, A_Y_LOW_AU, I_Y_DEG, M_Y_HIGH_EARTH, M_Y_LOW_EARTH,
};
use p9_2025_planet_y::laplace_plane::{
    PlanetY, forced_inclination_deg, inner_torque_coeff, outer_torque_coeff,
};
use p9_2025_planet_y::mean_plane::mean_plane_tilt_deg;
use serde_json::{Value, json};

/// Semi-major-axis bins of the paper's mean-plane measurement (AU) and whether
/// a warp was reported there.
const BINS: [(f64, f64, bool); 3] = [
    (50.0, 80.0, false),
    (80.0, 200.0, true),
    (200.0, 400.0, false),
];
/// Published warp of the 80-200 AU and 80-400 AU bins: inclination and node
/// (deg), and the confidence of each detection.
const PUBLISHED_WARP_DEG: f64 = 15.0;
const PUBLISHED_NODE_DEG: f64 = 120.0;
const PUBLISHED_CONFIDENCE: [(&str, f64); 2] = [("80-200 AU", 0.96), ("80-400 AU", 0.98)];
/// Non-resonant objects behind the measurement: 50-400 AU, and 80-400 AU.
const PUBLISHED_SAMPLE: (usize, usize) = (154, 46);
/// Inclination of the illustrated planets (deg): the paper prefers >= 10 deg
/// and the measured warp itself is ~15 deg.
const I_PLANET_DEG: f64 = 15.0;
/// Planets whose forced-plane profile is drawn: (mass in Earth masses, a in AU).
const CASES: [(f64, f64); 3] = [
    (M_Y_LOW_EARTH, 150.0),
    (0.3, 150.0),
    (M_Y_HIGH_EARTH, 150.0),
];
/// Profile grid (AU).
const A_MIN_AU: f64 = 40.0;
const A_MAX_AU: f64 = 400.0;
const N_PROFILE: usize = 181;
/// Grid of the (mass, semi-major axis) map.
const N_MASS: usize = 13;
const N_AXIS: usize = 11;
/// Samples per bin for the bin-averaged tilt.
const N_BIN: usize = 24;

/// Mean forced tilt (deg) over a bin, uniform in semi-major axis.
fn bin_tilt(lo: f64, hi: f64, planet: &PlanetY) -> f64 {
    (0..N_BIN)
        .map(|k| mean_plane_tilt_deg(lo + (hi - lo) * (k as f64 + 0.5) / N_BIN as f64, planet))
        .sum::<f64>()
        / N_BIN as f64
}

/// Semi-major axis (AU) inside the planet's orbit where its torque equals the
/// giant planets', by bisection.
fn crossover_au(planet: &PlanetY) -> f64 {
    let (mut lo, mut hi) = (30.0, planet.a_au);
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if outer_torque_coeff(mid, planet) > inner_torque_coeff(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

pub fn export() -> Value {
    let a_grid: Vec<f64> = (0..N_PROFILE)
        .map(|k| A_MIN_AU + (A_MAX_AU - A_MIN_AU) * k as f64 / (N_PROFILE - 1) as f64)
        .collect();
    let profiles: Vec<Value> = CASES
        .iter()
        .map(|&(m, a)| {
            let planet = PlanetY::new(m, a, I_PLANET_DEG);
            let tilt: Vec<f64> = a_grid
                .iter()
                .map(|&x| forced_inclination_deg(x, &planet))
                .collect();
            let bins: Vec<f64> = BINS
                .iter()
                .map(|&(lo, hi, _)| bin_tilt(lo, hi, &planet))
                .collect();
            json!({
                "mass_earth": m,
                "a_au": a,
                "i_deg": I_PLANET_DEG,
                "tilt_deg": tilt,
                "bin_tilt_deg": bins,
                "crossover_au": crossover_au(&planet),
            })
        })
        .collect();

    // The (mass, semi-major axis) map: tilt forced on the quiet 50-80 AU bin
    // and on the warped 80-200 AU bin.
    let masses: Vec<f64> = (0..N_MASS)
        .map(|k| 0.01 * (400.0_f64).powf(k as f64 / (N_MASS - 1) as f64))
        .collect();
    let axes: Vec<f64> = (0..N_AXIS)
        .map(|k| 80.0 + 220.0 * k as f64 / (N_AXIS - 1) as f64)
        .collect();
    let mut inner_bin = Vec::new();
    let mut warp_bin = Vec::new();
    for &a in &axes {
        for &m in &masses {
            let planet = PlanetY::new(m, a, I_PLANET_DEG);
            inner_bin.push(bin_tilt(BINS[0].0, BINS[0].1, &planet));
            warp_bin.push(bin_tilt(BINS[1].0, BINS[1].1, &planet));
        }
    }

    let earth = PlanetY::new(M_Y_HIGH_EARTH, 150.0, I_PLANET_DEG);
    let mercury = PlanetY::new(M_Y_LOW_EARTH, 150.0, I_PLANET_DEG);

    json!({
        "bins": BINS.iter().map(|&(lo, hi, warped)| json!({"lo_au": lo, "hi_au": hi, "warped": warped})).collect::<Vec<_>>(),
        "published": {
            "warp_deg": PUBLISHED_WARP_DEG,
            "node_deg": PUBLISHED_NODE_DEG,
            "confidence": PUBLISHED_CONFIDENCE.iter().map(|&(b, c)| json!({"bin": b, "confidence": c})).collect::<Vec<_>>(),
            "n_objects": PUBLISHED_SAMPLE.0,
            "n_objects_warp_bin": PUBLISHED_SAMPLE.1,
            "mass_range_earth": [M_Y_LOW_EARTH, M_Y_HIGH_EARTH],
            "a_range_au": [A_Y_LOW_AU, A_Y_HIGH_AU],
            "i_min_deg": I_Y_DEG,
        },
        "i_planet_deg": I_PLANET_DEG,
        "a_grid_au": a_grid,
        "profiles": profiles,
        "warp_bin_tilt_earth_deg": bin_tilt(BINS[1].0, BINS[1].1, &earth),
        "inner_bin_tilt_earth_deg": bin_tilt(BINS[0].0, BINS[0].1, &earth),
        "warp_bin_tilt_mercury_deg": bin_tilt(BINS[1].0, BINS[1].1, &mercury),
        "inner_bin_tilt_mercury_deg": bin_tilt(BINS[0].0, BINS[0].1, &mercury),
        "map": {
            "mass_earth": masses,
            "a_au": axes,
            "inner_bin_tilt_deg": inner_bin,
            "warp_bin_tilt_deg": warp_bin,
        },
    })
}
