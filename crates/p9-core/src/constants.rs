/// Gravitational constant in AU^3 / (M_sun * day^2)
/// Derived from GM_SUN = 1.32712440041279419e+20 m^3/s^2
/// with AU = 1.495978707e+11 m and day = 86400 s
pub const G_AU3_MSUN_DAY2: f64 = 2.959122082855911e-4;

/// GM values in AU^3/day^2 (heliocentric gravitational parameters)
/// Source: JPL DE440 / Park et al. (2021)
pub const GM_SUN: f64 = 2.959122082855911e-4;
pub const GM_MERCURY: f64 = 4.912_486_6e-11;
pub const GM_VENUS: f64 = 7.243_452_5e-10;
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

/// Earth mass in solar masses
pub const EARTH_MASS_SOLAR: f64 = 3.003_489e-6;

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

/// Unit conversions
pub const AU_M: f64 = 1.495_978_707e11;
pub const AU_KM: f64 = 1.495_978_707e8;
pub const DAY_S: f64 = 86_400.0;
pub const YEAR_DAYS: f64 = 365.25;
pub const GYR_DAYS: f64 = 365.25e9;

/// J2000.0 epoch (Julian Date, TT)
pub const J2000: f64 = 2_451_545.0;

/// Mathematical constants
pub const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
pub const DEG2RAD: f64 = std::f64::consts::PI / 180.0;
pub const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;

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
