//! Reproduction of Gomes, Deienno & Morbidelli (2016)
//! "The inclination of the planetary system relative to the solar equator
//!  may be explained by the presence of Planet 9" (arXiv:1607.05111).
//!
//! This is an *independent* calculation of the Planet-Nine solar-obliquity
//! mechanism, distinct in method from Bailey, Batygin & Brown (2016) (the
//! sibling crate `p9-2016-obliquity`) and from Lai (2016).
//!
//! ## The Gomes et al. mechanism
//!
//! A distant, inclined Planet Nine exerts a secular torque on the orbital
//! plane of the four giant planets. Because Planet Nine carries a comparable
//! (or larger) orbital angular momentum, the *total* angular momentum vector
//! of the system is tilted away from the original giant-planet plane. The
//! giant-planet angular-momentum vector L_p then *precesses* about the
//! conserved total angular momentum L_tot on a cone whose half-angle is the
//! mutual inclination i_p between L_p and L_tot.
//!
//! The Sun's spin axis is dynamically decoupled on these timescales: the
//! solar quadrupole couples the spin axis to the planetary plane only very
//! weakly (the spin precession period far exceeds the planet-plane precession
//! period), so to leading order the spin axis stays where it started — pointing
//! along the *original* planetary plane normal. As the planetary plane
//! precesses out from under it, the angle between the spin axis and the
//! instantaneous planetary plane (the observed "solar obliquity") grows.
//!
//! Over one full precession cycle the obliquity oscillates between 0 and
//! 2·i_p. The solar-system age is a fraction of the precession period, so the
//! achieved obliquity at 4.5 Gyr is 2·i_p·sin²(π·t/T_prec).
//!
//! The forced inclination i_p, the precession period T_prec, and hence the
//! obliquity all depend on Planet Nine's mass, semimajor axis, and
//! inclination. This crate computes those real quantities (reusing the giant
//! planet model and constants from `p9-core` and the sibling crate's solar
//! model) and pins them against the published ~6° figure.

pub mod precession_model;

/// Published reference constants from Gomes, Deienno & Morbidelli (2016).
pub mod reference {
    use p9_core::constants::DEG2RAD;
    use p9_core::units::{degrees, Angle};

    /// Observed solar obliquity: the tilt of the Sun's equator relative to the
    /// invariable (giant-planet) plane, ≈ 6° (degrees).
    pub const OBSERVED_SOLAR_OBLIQUITY_DEG: f64 = 6.0;

    /// Observed solar obliquity in radians.
    pub fn observed_solar_obliquity_rad() -> f64 {
        OBSERVED_SOLAR_OBLIQUITY_DEG * DEG2RAD
    }

    /// Observed solar obliquity as a typed [`Angle`].
    pub fn observed_solar_obliquity() -> Angle {
        degrees(OBSERVED_SOLAR_OBLIQUITY_DEG)
    }

    /// Nominal Planet Nine mass adopted by Gomes et al. (Earth masses).
    pub const P9_MASS_EARTH: f64 = 10.0;

    /// Nominal Planet Nine semimajor axis (AU).
    pub const P9_SEMIMAJOR_AXIS_AU: f64 = 700.0;

    /// Nominal Planet Nine eccentricity.
    pub const P9_ECCENTRICITY: f64 = 0.6;

    /// Nominal Planet Nine inclination relative to the original planetary
    /// plane (degrees).
    pub const P9_INCLINATION_DEG: f64 = 30.0;

    /// Age of the Solar System used as the integration horizon (Gyr).
    pub const SOLAR_SYSTEM_AGE_GYR: f64 = 4.5;
}
