//! Film export for `p9-2019-tess-yield`: the numbers its scene and ledger entry draw.

use p9_2019_tess_yield::tess::{
    TESS_FRAMES_PER_SECTOR, TESS_STACKED_DEPTH_IC, TESS_STACKED_DEPTH_SIGMA, TessStack,
};
use p9_2019_tess_yield::yield_grid::{ReferenceBox, p9_apparent_magnitude};
use p9_2020_tess_shiftstack::tess::{PIXEL_SCALE_ARCSEC, SECTOR_DAYS};
use p9_core::analysis::stacking::orbit_metric::apparent_sky_rate_at_opposition;
use serde_json::{Value, json};

/// Smallest displacement over a sector that separates a mover from a star
/// (pixels; Payne, Holman & Pál 2019).
const MIN_DISPLACEMENT_PX: f64 = 5.0;

/// Published distance limit set by that displacement (AU).
const PUBLISHED_MOTION_LIMIT_AU: f64 = 900.0;

fn displacement_px(distance_au: f64) -> f64 {
    apparent_sky_rate_at_opposition(distance_au) * SECTOR_DAYS / PIXEL_SCALE_ARCSEC
}

pub fn export() -> Value {
    let stack = TessStack::default();
    let bx = ReferenceBox::default();

    // Depth against the number of frames stacked, on a log grid.
    let log_frames: Vec<f64> = (0..=64)
        .map(|k| TESS_FRAMES_PER_SECTOR.log10() * k as f64 / 64.0)
        .collect();
    let depth: Vec<f64> = log_frames
        .iter()
        .map(|&lg| {
            TessStack {
                frames: 10f64.powf(lg),
                ..TessStack::default()
            }
            .stacked_depth()
        })
        .collect();

    let distances: Vec<f64> = (0..=80).map(|k| 200.0 + 10.0 * k as f64).collect();
    let v_at = |mass: f64| -> Vec<f64> {
        distances
            .iter()
            .map(|&d| p9_apparent_magnitude(mass, d, bx.albedo))
            .collect()
    };
    let shifted = |dm: f64| TessStack {
        single_frame_depth: stack.single_frame_depth + dm,
        ..TessStack::default()
    };
    let fraction_vs_distance = |s: &TessStack| -> Vec<f64> {
        distances
            .iter()
            .map(|&d| bx.detectable_fraction_at_distance(s, d))
            .collect()
    };

    // Distance at which a sector's drift falls to the minimum displacement.
    let (mut lo, mut hi) = (100.0_f64, 5000.0_f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if displacement_px(mid) > MIN_DISPLACEMENT_PX {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    json!({
        "single_frame_depth": stack.single_frame_depth,
        "frames_per_sector": stack.frames,
        "stacked_depth": stack.stacked_depth(),
        "published_depth": TESS_STACKED_DEPTH_IC,
        "published_depth_sigma": TESS_STACKED_DEPTH_SIGMA,
        "depth_vs_frames": {"log_frames": log_frames, "depth": depth},
        "box": {
            "mass_lo": bx.mass_earth.0,
            "mass_hi": bx.mass_earth.1,
            "distance_lo": bx.distance_au.0,
            "distance_hi": bx.distance_au.1,
            "albedo": bx.albedo,
        },
        "distance_au": distances,
        "v_mass_lo": v_at(bx.mass_earth.0),
        "v_mass_hi": v_at(bx.mass_earth.1),
        "detectable_fraction": bx.detectable_fraction(&stack),
        "detectable_fraction_shallow": bx.detectable_fraction(&shifted(-TESS_STACKED_DEPTH_SIGMA)),
        "detectable_fraction_deep": bx.detectable_fraction(&shifted(TESS_STACKED_DEPTH_SIGMA)),
        "fraction_vs_distance": fraction_vs_distance(&stack),
        "fraction_vs_distance_shallow": fraction_vs_distance(&shifted(-TESS_STACKED_DEPTH_SIGMA)),
        "fraction_vs_distance_deep": fraction_vs_distance(&shifted(TESS_STACKED_DEPTH_SIGMA)),
        "displacement_px": distances.iter().map(|&d| displacement_px(d)).collect::<Vec<f64>>(),
        "pixel_scale_arcsec": PIXEL_SCALE_ARCSEC,
        "sector_days": SECTOR_DAYS,
        "min_displacement_px": MIN_DISPLACEMENT_PX,
        "motion_limit_au": 0.5 * (lo + hi),
        "published_motion_limit_au": PUBLISHED_MOTION_LIMIT_AU,
    })
}
