/// Gravitational constant in AU^3 / (M_sun * day^2)
/// Derived from GM_SUN = 1.32712440041279419e+20 m^3/s^2
/// with AU = 1.495978707e+11 m and day = 86400 s
pub const G_AU3_MSUN_DAY2: f64 = 2.959122082855911e-4;

/// GM values in AU^3/day^2 (heliocentric gravitational parameters)
/// Source: JPL DE440 / Park et al. (2021)
///
/// Note: starfield's `GM_SUN` is in km^3/s^2 and from a slightly older
/// determination (1.32712440042e20 vs DE440's 1.32712440041279419e20 m^3/s^2,
/// a 5.4e-12 relative difference); p9-core keeps the DE440 value in its
/// AU^3/day^2 convention. The agreement is pinned by a test below.
pub const GM_SUN: f64 = 2.959122082855911e-4;
pub const GM_MERCURY: f64 = 4.912_486_6e-11;
pub const GM_VENUS: f64 = 7.243_452_5e-10;
/// GM of the Earth-Moon *system* (barycenter). Note this is NOT the same body
/// as `EARTH_MASS_SOLAR` below, which is the Earth alone (the Moon adds ~1.2%).
pub const GM_EARTH_MOON: f64 = 8.997_011_4e-10;
pub const GM_MARS: f64 = 9.549_535_2e-11;
pub const GM_JUPITER: f64 = 2.825_345_818e-7;
pub const GM_SATURN: f64 = 8.459_715_6e-8;
pub const GM_URANUS: f64 = 1.292_024_9e-8;
pub const GM_NEPTUNE: f64 = 1.524_358_9e-8;
pub const GM_PLUTO: f64 = 2.175_1e-12;

/// Mass ratios relative to the Sun
pub const MASS_JUPITER_SOLAR: f64 = 9.547_919_384e-4;
pub const MASS_SATURN_SOLAR: f64 = 2.858_859_81e-4;
pub const MASS_URANUS_SOLAR: f64 = 4.366_244_0e-5;
pub const MASS_NEPTUNE_SOLAR: f64 = 5.151_389_0e-5;

/// Earth mass in solar masses (Earth only, excluding the Moon; see GM_EARTH_MOON)
pub const EARTH_MASS_SOLAR: f64 = 3.003_489e-6;

/// Earth equatorial radius in km (IAU 2015 nominal)
pub const EARTH_RADIUS_KM: f64 = 6_378.1;

/// Equatorial radii in AU
pub const RADIUS_JUPITER_AU: f64 = 71_492.0 / 1.495_978_707e8;
pub const RADIUS_SATURN_AU: f64 = 60_268.0 / 1.495_978_707e8;
pub const RADIUS_URANUS_AU: f64 = 25_559.0 / 1.495_978_707e8;
pub const RADIUS_NEPTUNE_AU: f64 = 24_764.0 / 1.495_978_707e8;
pub const RADIUS_SUN_AU: f64 = 696_000.0 / 1.495_978_707e8;

/// J2 oblateness coefficients (dimensionless)
pub const J2_SUN: f64 = 2.0e-7;
pub const J2_JUPITER: f64 = 1.4736e-2;
pub const J2_SATURN: f64 = 1.6298e-2;
pub const J2_URANUS: f64 = 3.343e-3;
pub const J2_NEPTUNE: f64 = 3.411e-3;

/// J4 oblateness coefficients (dimensionless)
pub const J4_JUPITER: f64 = -5.87e-4;
pub const J4_SATURN: f64 = -9.35e-4;
pub const J4_URANUS: f64 = -3.21e-5;
pub const J4_NEPTUNE: f64 = -3.48e-5;

/// Unit conversions. AU_M, AU_KM, and DAY_S are re-exported from starfield
/// (IAU 2012 AU = 149_597_870_700 m exactly) so the two crates cannot drift.
pub use starfield::constants::{AU_KM, AU_M, DAY_S};
pub const YEAR_DAYS: f64 = 365.25;
pub const GYR_DAYS: f64 = 365.25e9;

/// J2000.0 epoch (Julian Date, TT); re-exported from starfield.
pub use starfield::constants::J2000;

/// Mathematical constants; re-exported from starfield where they overlap
/// (TAU is starfield's name for 2π).
pub use starfield::constants::TAU as TWO_PI;
pub use starfield::constants::{DEG2RAD, RAD2DEG};

/// Neptune's semi-major axis in AU (single source for the resonance/Hansen
/// machinery; several paper crates previously re-defined this inconsistently)
pub const A_NEPTUNE_AU: f64 = 30.07;

/// Neptune's orbital period in days (~164.8 years)
pub const NEPTUNE_PERIOD_DAYS: f64 = 60_190.03;

/// Jupiter's orbital period in days (~11.86 years)
pub const JUPITER_PERIOD_DAYS: f64 = 4_332.59;

/// Galactic tidal parameters
/// Local disk density in M_sun/AU^3 (0.1 M_sun/pc^3 converted)
pub const RHO_MW_MSUN_AU3: f64 = 0.1 / {
    let pc_km = 3.085_677_6e13_f64;
    let ratio = pc_km / AU_KM;
    ratio * ratio * ratio
};

/// Parsec in AU
pub const PC_AU: f64 = 206_264.806;

/// km/s to AU/day
pub const KMS_TO_AUDAY: f64 = DAY_S / AU_KM;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gm_sun_au3_day2_matches_starfield_km3_s2() {
        // starfield's GM_SUN (km^3/s^2) converted to p9-core's AU^3/day^2
        // convention. The underlying solar GM determinations differ by
        // 5.4e-12 (DE440 vs starfield's older value), so agreement is pinned
        // at 1e-11 relative rather than full f64 precision.
        let starfield_au3_day2 =
            starfield::constants::GM_SUN * DAY_S * DAY_S / (AU_KM * AU_KM * AU_KM);
        let rel = ((starfield_au3_day2 - GM_SUN) / GM_SUN).abs();
        assert!(rel < 1e-11, "relative difference = {rel:e}");
    }

    #[test]
    fn test_g_au3_msun_day2_consistent_with_gm_sun() {
        // G in AU^3/(M_sun day^2) is numerically GM_sun in AU^3/day^2.
        assert_eq!(G_AU3_MSUN_DAY2, GM_SUN);
    }

    #[test]
    fn test_au_conversions_self_consistent() {
        // Re-exported starfield values: AU_M and AU_KM describe the same
        // IAU 2012 astronomical unit.
        assert_eq!(AU_M, 149_597_870_700.0);
        assert_eq!(AU_KM * 1_000.0, AU_M);
    }
}
