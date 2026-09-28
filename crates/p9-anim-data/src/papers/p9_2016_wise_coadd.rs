//! Film export for `p9-2016-wise-coadd`: the numbers its scene and ledger entry draw.

use p9_2016_holman_payne::published::{
    PREFERRED_DEC_DEG, PREFERRED_RA_DEG, PREFERRED_SKY_HALF_EXTENT_DEG,
};
use p9_2016_wise_coadd::band::{W1, magnitude};
use p9_2016_wise_coadd::coadd::{
    PUBLISHED_COADD_W1_DEPTH, SINGLE_FRAME_W1_DEPTH, coadd_depth, frames_for_gain,
};
use p9_2016_wise_coadd::detectability::max_detectable_distance;
use p9_2018_wise_search::thermal_model::{
    P9Thermal, W1_ALBEDO_BRIGHT, W1_ALBEDO_DARK, W1_ALBEDO_DEFAULT,
};
use p9_core::analysis::thermal::max_detectable_distance as crossing_distance;
use serde_json::{Value, json};

/// Mass of the planet whose brightness is traced (Earth masses).
const MASS_EARTH: f64 = 10.0;

/// The most W1-luminous Fortney et al. (2016) model atmosphere as quoted by
/// Meisner et al. (2016): W1 = 16.1 at 622 AU. Self-luminous, so it fades as
/// the inverse square of distance.
const LUMINOUS_W1: f64 = 16.1;
const LUMINOUS_AT_AU: f64 = 622.0;

/// Published reach with single exposures and with the coadds (AU).
const PUBLISHED_SINGLE_REACH_AU: f64 = 430.0;
const PUBLISHED_COADD_REACH_AU: f64 = 800.0;

/// Area actually searched (square degrees; Meisner et al. 2016, Fig. 1).
const PUBLISHED_AREA_DEG2: f64 = 1840.0;

fn luminous_w1(distance_au: f64) -> f64 {
    LUMINOUS_W1 + 5.0 * (distance_au / LUMINOUS_AT_AU).log10()
}

pub fn export() -> Value {
    let frames: Vec<u32> = (1..=40).collect();
    let depth: Vec<f64> = frames.iter().map(|&n| coadd_depth(W1, n)).collect();
    let frames_per_coadd = frames_for_gain(PUBLISHED_COADD_W1_DEPTH - SINGLE_FRAME_W1_DEPTH);
    let coadd = coadd_depth(W1, frames_per_coadd.round() as u32);

    let distances: Vec<f64> = (0..=90).map(|k| 100.0 + 10.0 * k as f64).collect();
    let reflected = |albedo: f64| -> Vec<f64> {
        distances
            .iter()
            .map(|&d| magnitude(&P9Thermal::new(MASS_EARTH, d).with_albedo(albedo), W1))
            .collect()
    };
    let luminous: Vec<f64> = distances.iter().map(|&d| luminous_w1(d)).collect();

    let half = PREFERRED_SKY_HALF_EXTENT_DEG;
    let (dec_lo, dec_hi) = (PREFERRED_DEC_DEG - half, PREFERRED_DEC_DEG + half);
    let box_area =
        2.0 * half * (dec_hi.to_radians().sin() - dec_lo.to_radians().sin()).to_degrees();

    json!({
        "mass_earth": MASS_EARTH,
        "single_depth": SINGLE_FRAME_W1_DEPTH,
        "coadd_depth": coadd,
        "published_coadd_depth": PUBLISHED_COADD_W1_DEPTH,
        "frames_per_coadd": frames_per_coadd,
        "depth_vs_frames": {"frames": frames, "w1_depth": depth},
        "distance_au": distances,
        "w1_reflected": {
            "dark": reflected(W1_ALBEDO_DARK),
            "default": reflected(W1_ALBEDO_DEFAULT),
            "bright": reflected(W1_ALBEDO_BRIGHT),
            "albedo_dark": W1_ALBEDO_DARK,
            "albedo_default": W1_ALBEDO_DEFAULT,
            "albedo_bright": W1_ALBEDO_BRIGHT,
        },
        "w1_luminous": luminous,
        "luminous_anchor": {"w1": LUMINOUS_W1, "distance_au": LUMINOUS_AT_AU},
        "reach_reflected_single_au": max_detectable_distance(W1, SINGLE_FRAME_W1_DEPTH, MASS_EARTH),
        "reach_reflected_coadd_au": max_detectable_distance(W1, coadd, MASS_EARTH),
        "reach_luminous_single_au":
            crossing_distance(50.0, 5000.0, SINGLE_FRAME_W1_DEPTH, luminous_w1),
        "reach_luminous_coadd_au": crossing_distance(50.0, 5000.0, coadd, luminous_w1),
        "published_single_reach_au": PUBLISHED_SINGLE_REACH_AU,
        "published_coadd_reach_au": PUBLISHED_COADD_REACH_AU,
        "region": {
            "ra_lo": PREFERRED_RA_DEG - half,
            "ra_hi": PREFERRED_RA_DEG + half,
            "dec_lo": dec_lo,
            "dec_hi": dec_hi,
            "area_deg2": box_area,
            "published_area_deg2": PUBLISHED_AREA_DEG2,
        },
    })
}
