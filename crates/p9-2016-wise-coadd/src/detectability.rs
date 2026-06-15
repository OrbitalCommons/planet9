//! Maximum heliocentric distance to which a coadd detects a given-mass P9, per
//! band.
//!
//! As in the sibling crate, band magnitude faints monotonically with distance,
//! so for any mass bright enough at the near edge there is a single crossing of
//! the (coadd) depth, found by bisection. The 2016 paper's deep coadd reaches
//! larger distances than single frames precisely because the depth is deeper by
//! the √N gain (`coadd`). The W2-vs-W1 contrast — W2 reaching a cold body
//! farther out — comes from the band flux (`band`): on the Wien tail W2 is far
//! brighter for a cold planet.

use p9_2018_wise_search::thermal_model::P9Thermal;

use crate::band::{band_magnitude, WiseBand};

/// Whether a P9 of `mass_earth` at `distance_au` is brighter than `depth` in
/// `band` (modelled band magnitude at or below the depth).
pub fn is_detectable(band: WiseBand, depth: f64, mass_earth: f64, distance_au: f64) -> bool {
    band_magnitude(mass_earth, distance_au, band) <= depth
}

/// Maximum heliocentric distance (AU) at which a P9 of `mass_earth` reaches
/// `depth` in `band`, by bisection on the monotone magnitude–distance relation.
/// `None` if undetectable even at the near edge of the search range.
pub fn max_detectable_distance(band: WiseBand, depth: f64, mass_earth: f64) -> Option<f64> {
    let d_min = 50.0_f64;
    let d_max = 50_000.0_f64;

    if !is_detectable(band, depth, mass_earth, d_min) {
        return None;
    }
    if is_detectable(band, depth, mass_earth, d_max) {
        return Some(d_max);
    }

    let mut lo = d_min;
    let mut hi = d_max;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if is_detectable(band, depth, mass_earth, mid) {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-4 {
            break;
        }
    }
    Some(0.5 * (lo + hi))
}

/// Detectability of a fully specified planet against a band depth (uses the
/// planet's own albedo / internal temperature).
pub fn planet_is_detectable(band: WiseBand, depth: f64, p9: &P9Thermal) -> bool {
    crate::band::magnitude(p9, band) <= depth
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::band::{W1, W2};
    use crate::coadd::{coadd_depth, single_frame_depth};

    #[test]
    fn max_distance_increases_with_mass() {
        let depth = coadd_depth(W1, 100);
        let d_heavy = max_detectable_distance(W1, depth, 15.0).expect("heavy detectable");
        let d_light = max_detectable_distance(W1, depth, 8.0).expect("light detectable");
        assert!(
            d_heavy > d_light,
            "heavier planet reaches farther: {d_heavy:.0} vs {d_light:.0} AU"
        );
    }

    #[test]
    fn deeper_coadd_reaches_farther_than_single_frame() {
        // The whole point of the 2016 coadd: a deeper limit reaches a larger
        // maximum distance than the single-frame depth, at fixed mass/band.
        let mass = 12.0;
        let single = max_detectable_distance(W1, single_frame_depth(W1), mass)
            .expect("detectable single-frame");
        let deep =
            max_detectable_distance(W1, coadd_depth(W1, 100), mass).expect("detectable coadd");
        assert!(
            deep > single,
            "coadd reach {deep:.0} AU must exceed single-frame {single:.0} AU"
        );
    }

    #[test]
    fn w2_reaches_a_cold_p9_farther_than_w1() {
        // Headline of the paper: for a cold body W2 is more favorable. At the
        // SAME numerical depth, W2 reaches the planet at larger distance than
        // W1 because the cold-body band flux is higher at 4.6 µm.
        let mass = 10.0;
        let depth = 18.0; // common reference depth for an apples-to-apples band test
        let d_w1 = max_detectable_distance(W1, depth, mass).expect("W1 detectable");
        let d_w2 = max_detectable_distance(W2, depth, mass).expect("W2 detectable");
        assert!(
            d_w2 > d_w1,
            "W2 reach {d_w2:.0} AU should exceed W1 reach {d_w1:.0} AU for a cold P9"
        );
    }

    #[test]
    fn crossover_magnitude_equals_depth() {
        let depth = coadd_depth(W1, 100);
        let d = max_detectable_distance(W1, depth, 12.0).expect("detectable");
        let m = band_magnitude(12.0, d, W1);
        assert!((m - depth).abs() < 1e-2, "W1({d:.1} AU) = {m:.3}");
    }
}
