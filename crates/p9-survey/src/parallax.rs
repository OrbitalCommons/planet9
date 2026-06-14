//! Observer-baseline parallax of a distant body.
//!
//! Parallax is just baseline / distance. For a body at hundreds of AU the only
//! baseline that matters is Earth's *heliocentric* orbit (2 AU peak-to-peak),
//! which gives Planet Nine a huge annual parallax (~arcminutes). A low-Earth
//! orbit adds a second, far shorter baseline — the orbit *diameter*,
//! ~14,000 km ≈ 9.4e-5 AU — so its per-orbit parallax is ~2e4× smaller than the
//! annual one and lands well below an arcsecond (sub-pixel for the JBT/SPENCER
//! 0.130″/px scale). Being in LEO helps by getting above the atmosphere
//! (depth, continuous all-sky access), not by parallax.
//!
//! This module quantifies that so the comparison is explicit, not assumed.

/// Astronomical unit in kilometers (IAU 2012).
pub const AU_KM: f64 = 1.495_978_707e8;
/// Mean Earth radius used elsewhere in the workspace.
pub const EARTH_RADIUS_KM: f64 = p9_core::constants::EARTH_RADIUS_KM;
/// Arcseconds per radian (= PC_AU, the parsec definition).
pub const ARCSEC_PER_RAD: f64 = p9_core::constants::PC_AU;
/// The LEO altitude requested for this comparison (km).
pub const LEO_ALTITUDE_KM: f64 = 650.0;

/// Full peak-to-peak parallax (arcsec) of a body at `distance_au` seen across a
/// baseline of `baseline_km`.
pub fn parallax_arcsec(baseline_km: f64, distance_au: f64) -> f64 {
    (baseline_km / AU_KM) / distance_au * ARCSEC_PER_RAD
}

/// LEO orbit diameter (km) at the given altitude — the maximum instantaneous
/// parallax baseline a single satellite sweeps each orbit.
pub fn leo_baseline_km(altitude_km: f64) -> f64 {
    2.0 * (EARTH_RADIUS_KM + altitude_km)
}

/// Per-orbit (peak-to-peak) parallax of a body at `distance_au`, observed from
/// a circular LEO at `altitude_km`.
pub fn leo_parallax_arcsec(altitude_km: f64, distance_au: f64) -> f64 {
    parallax_arcsec(leo_baseline_km(altitude_km), distance_au)
}

/// Annual (peak-to-peak) parallax from Earth's 2-AU heliocentric baseline.
pub fn annual_parallax_arcsec(distance_au: f64) -> f64 {
    parallax_arcsec(2.0 * AU_KM, distance_au)
}

/// Distance (AU) at which the LEO per-orbit parallax equals `arcsec` — e.g.
/// pass one detector pixel to find where the wobble becomes a 1-pixel effect.
pub fn leo_parallax_crossover_au(altitude_km: f64, arcsec: f64) -> f64 {
    (leo_baseline_km(altitude_km) / AU_KM) * ARCSEC_PER_RAD / arcsec
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn leo_baseline_is_orbit_diameter() {
        // 650 km altitude -> radius 7028.1 km -> diameter 14056.2 km.
        assert_relative_eq!(leo_baseline_km(650.0), 14_056.2, epsilon = 0.1);
    }

    #[test]
    fn leo_parallax_is_tens_of_milliarcsec() {
        // At a representative 460 AU the per-orbit swing is ~0.042".
        let p = leo_parallax_arcsec(650.0, 460.0);
        assert_relative_eq!(p, 0.0421, epsilon = 0.002);
        // Sub-pixel for the JBT/SPENCER 0.130"/px scale.
        assert!(p < 0.130);
    }

    #[test]
    fn annual_parallax_dwarfs_leo() {
        // Annual parallax at 460 AU is ~897" (15'); ratio ~ 2 AU / 14,056 km.
        let ann = annual_parallax_arcsec(460.0);
        let leo = leo_parallax_arcsec(650.0, 460.0);
        assert_relative_eq!(ann, 897.0, epsilon = 5.0);
        let ratio = ann / leo;
        assert!(
            ratio > 2.0e4,
            "annual/LEO ratio = {ratio:.0}, expected ~2e4"
        );
    }

    #[test]
    fn crossover_is_far_inside_planet_nine() {
        // LEO wobble reaches one JBT pixel (0.130") only inside ~149 AU —
        // i.e. never for P9 at 300-800 AU, but yes for nearer scattered-disk
        // objects.
        let d = leo_parallax_crossover_au(650.0, 0.130);
        assert_relative_eq!(d, 149.0, epsilon = 3.0);
        assert!(d < 300.0);
    }
}
