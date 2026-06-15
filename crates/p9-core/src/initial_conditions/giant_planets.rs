//! Coplanar giant-planet "wire" set and angular-momentum helpers.
//!
//! A flat secular representation of the four giant planets as
//! `(mass [M_sun], semi-major axis [AU])` pairs, distinct from the full
//! J2000 state vectors in [`super::planets`]: the secular obliquity/precession
//! theories (Bailey+ 2016, Gomes+ 2016, Lai 2016) treat the planets as rigid
//! coplanar wires and only need mass and semi-major axis.
//!
//! Masses come from p9-core's [`MASS_JUPITER_SOLAR`] … [`MASS_NEPTUNE_SOLAR`]
//! constants (the same DE440 planetary-system masses used by
//! [`super::planets`]); Neptune's semi-major axis is [`A_NEPTUNE_AU`]. The
//! Jupiter/Saturn/Uranus semi-major axes are the secular-theory values used by
//! the obliquity reproductions and are canonical here for those crates.

use crate::constants::*;

/// The four giant planets as (mass in M_sun, semi-major axis in AU) — the
/// rigid coplanar "wire" set whose couplings are summed *per planet*
/// (Bailey+ 2016). Collapsing them to a single angular-momentum-conserving
/// ring underestimates Σmᵢ/aᵢ³-type couplings by ~40% and Σmᵢaᵢ² ones
/// by ~45%.
pub const GIANT_PLANETS: [(f64, f64); 4] = [
    (MASS_JUPITER_SOLAR, 5.2026),
    (MASS_SATURN_SOLAR, 9.5549),
    (MASS_URANUS_SOLAR, 19.2184),
    (MASS_NEPTUNE_SOLAR, A_NEPTUNE_AU),
];

/// Combined angular momentum of the four giant planets.
///
/// L_planets = Σ m_i * sqrt(G * M_sun * a_i) for each giant planet.
/// Returns in M_sun * AU² / day.
pub fn giant_planet_angular_momentum() -> f64 {
    GIANT_PLANETS
        .iter()
        .map(|&(m, a)| m * (G_AU3_MSUN_DAY2 * a).sqrt())
        .sum()
}

/// Planet Nine angular momentum magnitude.
///
/// L₉ = m₉ * sqrt(G * M_sun * a₉) * sqrt(1 - e₉²)
pub fn p9_angular_momentum(mass_solar: f64, a: f64, e: f64) -> f64 {
    mass_solar * (G_AU3_MSUN_DAY2 * a).sqrt() * (1.0 - e * e).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_giant_planet_angular_momentum() {
        let l = giant_planet_angular_momentum();
        // Should be positive and dominated by Jupiter
        assert!(l > 0.0);

        let l_jup = MASS_JUPITER_SOLAR * (G_AU3_MSUN_DAY2 * 5.2026_f64).sqrt();
        // Jupiter should be > 50% of total
        assert!(l_jup / l > 0.5);
    }

    #[test]
    fn test_p9_angular_momentum() {
        let m9 = 10.0 * EARTH_MASS_SOLAR;
        let l9 = p9_angular_momentum(m9, 700.0, 0.6);
        let l_planets = giant_planet_angular_momentum();

        // P9's angular momentum is significant but smaller than giant planets
        // At a=700 AU, L_9/L_gp ~ 0.2 (large a compensates small mass)
        assert!(l9 < l_planets);
        assert!(l9 > 0.0);
    }
}
