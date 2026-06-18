//! The "canonical" solar system planets (Mercury → Neptune) that Lai (2016)
//! sums over for the total orbital angular momentum `L`, the P9-forced
//! precession rate `Ω_L`, and the solar spin precession `Ω★ps`.
//!
//! The four giants reuse
//! [`p9_core::initial_conditions::giant_planets::GIANT_PLANETS`]; the four
//! terrestrials are added here (they contribute < 0.1% of `L` and `Ω_L` but
//! matter at the percent level for `Ω★ps ∝ Σ mⱼ/aⱼ³`, which is weighted toward
//! the inner planets relative to the angular-momentum sums).

use p9_core::constants::*;
use p9_core::initial_conditions::giant_planets::GIANT_PLANETS;
use p9_core::units::{days, radians, AngularVelocity};

/// Terrestrial planet (mass in M_sun, semi-major axis in AU). Masses from
/// JPL DE440 GM ratios; semi-major axes are mean osculating values.
pub const TERRESTRIAL_PLANETS: [(f64, f64); 4] = [
    (GM_MERCURY / GM_SUN, 0.387_098),
    (GM_VENUS / GM_SUN, 0.723_332),
    // Earth alone (the Moon adds ~1.2% but is negligible for these sums).
    (EARTH_MASS_SOLAR, 1.000_000),
    (GM_MARS / GM_SUN, 1.523_679),
];

/// All eight canonical planets, Mercury → Neptune, as (mass M_sun, a AU).
pub fn canonical_planets() -> Vec<(f64, f64)> {
    TERRESTRIAL_PLANETS
        .iter()
        .copied()
        .chain(GIANT_PLANETS.iter().copied())
        .collect()
}

/// Jupiter's semi-major axis (AU) — the angular-momentum normalization unit
/// in Lai (2016).
pub const A_JUPITER_AU: f64 = 5.2026;

/// Mean motion `n = √(GM_sun / a³)` in rad/day for a circular orbit at `a`.
pub fn mean_motion(a_au: f64) -> f64 {
    (G_AU3_MSUN_DAY2 * (a_au * a_au * a_au).recip()).sqrt()
}

/// [`mean_motion`] as a typed [`AngularVelocity`].
pub fn mean_motion_typed(a_au: f64) -> AngularVelocity {
    (radians(mean_motion(a_au)) / days(1.0)).into()
}

/// Orbital angular momentum of a circular planet, `L = m √(GM_sun a)`,
/// in M_sun·AU²/day.
pub fn orbital_angular_momentum(mass_solar: f64, a_au: f64) -> f64 {
    mass_solar * (G_AU3_MSUN_DAY2 * a_au).sqrt()
}

/// Jupiter's orbital angular momentum `L_J` (the normalization in Lai 2016).
pub fn l_jupiter() -> f64 {
    orbital_angular_momentum(MASS_JUPITER_SOLAR, A_JUPITER_AU)
}

/// Total canonical-planet orbital angular momentum `L = Σ_j Lⱼ`.
/// Lai (2016) quotes `L = 1.624 L_J`.
pub fn total_orbital_angular_momentum() -> f64 {
    canonical_planets()
        .iter()
        .map(|&(m, a)| orbital_angular_momentum(m, a))
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use uom::si::angular_velocity::radian_per_second;

    #[test]
    fn typed_mean_motion_matches_f64() {
        assert_relative_eq!(
            mean_motion_typed(A_JUPITER_AU).get::<radian_per_second>(),
            mean_motion(A_JUPITER_AU) / 86_400.0,
            epsilon = 1e-30
        );
    }

    #[test]
    fn total_l_matches_lai_1p624_lj() {
        // Lai (2016) Eq. (4)/(5): L = 1.624 L_J.
        let ratio = total_orbital_angular_momentum() / l_jupiter();
        assert!(
            (ratio - 1.624).abs() < 0.01,
            "L/L_J = {ratio:.4} (Lai: 1.624)"
        );
    }

    #[test]
    fn giants_dominate_angular_momentum() {
        let l_total = total_orbital_angular_momentum();
        let l_terr: f64 = TERRESTRIAL_PLANETS
            .iter()
            .map(|&(m, a)| orbital_angular_momentum(m, a))
            .sum();
        // Terrestrials are ~0.16% of the canonical L (giants dominate).
        let share = l_terr / l_total;
        assert!(share < 3e-3, "terrestrial L share {share:.3e}");
    }

    #[test]
    fn jupiter_mean_motion_period_about_11_86_yr() {
        let n = mean_motion(A_JUPITER_AU);
        let period_yr = (2.0 * std::f64::consts::PI / n) / YEAR_DAYS;
        assert!(
            (period_yr - 11.86).abs() < 0.05,
            "Jupiter period {period_yr:.3} yr"
        );
    }
}
