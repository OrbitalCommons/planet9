//! Reproduction of Fienga et al. (2016), "Constraints on the location of a
//! possible 9th planet derived from the Cassini data" (arXiv:1602.06116).
//!
//! # Headline result
//!
//! A 10 Earth-mass Planet Nine on the Brown & Batygin (2016) orbit perturbs
//! the heliocentric Sun–Saturn vector. Folded through the INPOP planetary
//! ephemeris and confronted with the Cassini Earth–Saturn range residuals,
//! the data REDUCE the residuals for a P9 located near true anomaly
//! `ν ≈ 108°` along the B&B orbit (a heliocentric distance of ~600 AU at that
//! phase) and WORSEN them (exclude P9) over large complementary ranges of ν.
//! The Cassini ranging thus localizes a preferred zone and excludes others.
//!
//! # What this crate computes (real, no hard-coded answers)
//!
//! 1. [`geometry`]: the P9 heliocentric position `r_p9(ν)` as a function of its
//!    own true anomaly along the B&B orbit, using `p9_core`'s element→Cartesian
//!    machinery. For `e ~ 0.6` the heliocentric distance `r(ν)` swings from the
//!    perihelion `a(1-e)` to the aphelion `a(1+e)`.
//! 2. [`perturbation`]: the DIFFERENTIAL (tidal) acceleration P9 exerts on
//!    Saturn relative to the Sun,
//!    `a_diff = GM_p9 · [ (r_p9 - r_sat)/|r_p9 - r_sat|³ - r_p9/|r_p9|³ ]`,
//!    and its quasi-static conversion to an Earth–Saturn RANGE perturbation
//!    amplitude across `ν ∈ [0, 360°)`.
//! 3. The ν that MINIMIZES the range perturbation (the best-fit zone), shown to
//!    sit near `ν ≈ 100–120°`, with antipodal/complementary ν giving large
//!    (excluded) perturbations.
//!
//! # Approximations (this is an analytic proxy, not the full INPOP fit)
//!
//! Fienga et al. re-fit the entire INPOP planetary ephemeris (all planet
//! orbits, Cassini tracking arcs, relativity) for each assumed P9 position and
//! read off the change in the post-fit Earth–Saturn range residuals. We instead
//! use a quasi-static analytic proxy: the P9 differential acceleration on
//! Saturn is projected onto the Earth–Saturn line of sight and integrated over a
//! characteristic Cassini observing baseline (`CASSINI_BASELINE_YEARS`) to give
//! the induced range displacement `δρ ≈ ½ |a_los| T²`. This reproduces the
//! SHAPE of the constraint (a single preferred ν zone, with a clear minimum
//! near 108° and large excluded lobes) and the SCALINGS (linear in `M_p9`,
//! steeply falling with P9 heliocentric distance), which is what the tests pin.
//! It does not reproduce the absolute residual amplitudes of the full fit, and
//! a real ephemeris re-fit partially absorbs a constant-acceleration term into
//! the planets' osculating elements (we do not). See the test module and
//! `REPRODUCTION_NOTES.md` for the honest residual discussion.

pub mod geometry;
pub mod perturbation;

use p9_core::constants::DEG2RAD;
use p9_core::types::P9Params;

/// Planet Nine on the Brown & Batygin (2016) orbit, exactly as tabulated by
/// Fienga et al. (2016), Table 1: `a = 700 AU`, `e = 0.6`, `i = 30°`,
/// `ω = 150°`, `Ω = 113°` (IERS ecliptic). The mean anomaly is irrelevant here
/// — the analysis sweeps the true anomaly explicitly. Note `Ω = 113°` differs
/// from `P9Params::nominal_2016()` (`Ω = 100°`), so this paper-specific orbit
/// is built directly.
pub fn brown_batygin_orbit() -> P9Params {
    P9Params {
        mass_earth: published::P9_MASS_EARTH,
        a: 700.0,
        e: 0.6,
        i: 30.0 * DEG2RAD,
        omega: 150.0 * DEG2RAD,
        omega_big: 113.0 * DEG2RAD,
        mean_anomaly: 0.0,
    }
}

/// Published reference values from Fienga et al. (2016), kept as labelled
/// constants (NOT used to derive any computed quantity).
pub mod published {
    /// Planet Nine mass assumed by the paper, in Earth masses.
    pub const P9_MASS_EARTH: f64 = 10.0;

    /// Most-probable true anomaly of P9 that REDUCES the Cassini Earth–Saturn
    /// range residuals (degrees): `v = 117.8°₋₁₀⁺¹¹` (Fienga et al. 2016, §4).
    pub const PREFERRED_TRUE_ANOMALY_DEG: f64 = 117.8;

    /// Favored interval for the true anomaly (degrees): `v ∈ [108°, 129°]`
    /// (Fienga et al. 2016, Conclusions / Fig 6 green zone).
    pub const FAVORED_INTERVAL_DEG: (f64, f64) = (108.0, 129.0);

    /// Approximate P9 heliocentric distance at the preferred true anomaly (AU).
    pub const PREFERRED_HELIO_DISTANCE_AU: f64 = 600.0;

    /// Tolerance (deg) within which our analytic-proxy favored anomaly should
    /// land of the published preferred true anomaly.
    pub const TRUE_ANOMALY_TOLERANCE_DEG: f64 = 30.0;
}

pub use geometry::{p9_position_at_true_anomaly, P9Geometry};
pub use perturbation::{
    differential_acceleration, favored_true_anomaly, prefit_amplitude,
    range_perturbation_amplitude, range_perturbation_curve, RangePerturbationCurve,
};
