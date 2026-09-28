//! Film export for `p9-2022-uranus-tilt`: the numbers its scene and ledger entry draw.

use p9_2022_uranus_tilt::cassini::{cassini_state_2, critical_ratio};
use p9_2022_uranus_tilt::resonance_capture::{
    SweepConfig, final_obliquity, fixed_band_sweep, max_obliquity, run_sweep, uranus_sweep,
};
use p9_2022_uranus_tilt::spin_axis::{alpha_typed, precession_period_myr};
use p9_2022_uranus_tilt::{ALPHA_PRESENT_ARCSEC_PER_YR, URANUS_OBLIQUITY_DEG};
use p9_core::units::{arcseconds, julian_year};
use serde_json::{Value, json};
use std::thread;

/// Spin-axis precession constant of the paper's showcase run (arcsec/yr),
/// Lu & Laughlin (2022) Fig. 9, and the peak obliquity that run reached.
const SHOWCASE_ALPHA_ARCSEC: f64 = 5.305;
const PUBLISHED_SHOWCASE_PEAK_DEG: f64 = 103.49;

/// Range of peak obliquities across the paper's stable runs (degrees).
const PUBLISHED_PEAK_RANGE_DEG: (f64, f64) = (75.9, 105.6);

/// Range of precession constants the paper's spin integrations sample (arcsec/yr).
const PUBLISHED_ALPHA_RANGE_ARCSEC: (f64, f64) = (0.005, 6.0);

/// Snapshot spacing of the exported histories (yr).
const SNAPSHOT_YR: f64 = 1.0e6;

fn sweep(cfg: &SweepConfig) -> Value {
    let snaps = run_sweep(cfg, SNAPSHOT_YR);
    json!({
        "t_myr": snaps.iter().map(|s| s.t / 1.0e6).collect::<Vec<_>>(),
        "obliquity_deg": snaps.iter().map(|s| s.obliquity.to_degrees()).collect::<Vec<_>>(),
        "ratio": snaps.iter().map(|s| s.ratio).collect::<Vec<_>>(),
        "final_deg": final_obliquity(&snaps).to_degrees(),
        "peak_deg": max_obliquity(&snaps).to_degrees(),
    })
}

pub fn export() -> Value {
    let alpha_today = (alpha_typed() / (arcseconds(1.0) / julian_year())).value;

    // The paper's showcase precession constant, swept through the resonance.
    let showcase_cfg = uranus_sweep(SHOWCASE_ALPHA_ARCSEC);
    let inclination = showcase_cfg.inclination;
    let showcase = sweep(&showcase_cfg);

    // The equilibrium the axis follows: Cassini state 2 against alpha/|g|.
    let ratio: Vec<f64> = (0..=120)
        .map(|k| 10f64.powf(-1.0 + 2.4 * k as f64 / 120.0))
        .collect();
    let state2: Vec<f64> = ratio
        .iter()
        .map(|&r| cassini_state_2(r, inclination).to_degrees())
        .collect();

    // The same Planet Nine-driven band of nodal frequencies, with only the
    // precession constant varied.
    let alphas: Vec<f64> = (0..=24)
        .map(|k| 10f64.powf(-2.0 + 2.8 * k as f64 / 24.0))
        .chain([alpha_today, SHOWCASE_ALPHA_ARCSEC])
        .collect();
    let mut band: Vec<(f64, f64, f64)> = thread::scope(|scope| {
        let handles: Vec<_> = alphas
            .iter()
            .map(|&a| {
                scope.spawn(move || {
                    let snaps = run_sweep(&fixed_band_sweep(a), SNAPSHOT_YR);
                    (
                        a,
                        final_obliquity(&snaps).to_degrees(),
                        max_obliquity(&snaps).to_degrees(),
                    )
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("sweep"))
            .collect()
    });
    band.sort_by(|x, y| x.0.total_cmp(&y.0));
    let at = |alpha: f64| {
        band.iter()
            .find(|b| b.0 == alpha)
            .map(|b| b.2)
            .expect("alpha was swept")
    };
    let band_cfg = fixed_band_sweep(alpha_today);
    let arcsec_per_rad = 3600.0 * 180.0 / std::f64::consts::PI;

    json!({
        "alpha_today_arcsec_yr": alpha_today,
        "precession_period_today_myr": precession_period_myr(0.0),
        "uranus_obliquity_deg": URANUS_OBLIQUITY_DEG,
        "inclination_deg": inclination.to_degrees(),
        "critical_ratio": critical_ratio(inclination),
        "showcase_alpha_arcsec_yr": SHOWCASE_ALPHA_ARCSEC,
        "alpha_enhancement": SHOWCASE_ALPHA_ARCSEC / alpha_today,
        "showcase": showcase,
        "showcase_peak_deg": showcase["peak_deg"],
        "cassini_state_2": {"ratio": ratio, "obliquity_deg": state2},
        "band": {
            "g_initial_arcsec_yr": band_cfg.g_initial.abs() * arcsec_per_rad,
            "g_final_arcsec_yr": band_cfg.g_final.abs() * arcsec_per_rad,
            "alpha_arcsec_yr": band.iter().map(|b| b.0).collect::<Vec<_>>(),
            "final_deg": band.iter().map(|b| b.1).collect::<Vec<_>>(),
            "peak_deg": band.iter().map(|b| b.2).collect::<Vec<_>>(),
            "peak_today_deg": at(alpha_today),
            "peak_showcase_deg": at(SHOWCASE_ALPHA_ARCSEC),
        },
        "published": {
            "alpha_today_arcsec_yr": ALPHA_PRESENT_ARCSEC_PER_YR,
            "showcase_peak_deg": PUBLISHED_SHOWCASE_PEAK_DEG,
            "peak_range_deg": PUBLISHED_PEAK_RANGE_DEG,
            "alpha_range_arcsec_yr": PUBLISHED_ALPHA_RANGE_ARCSEC,
        },
    })
}
