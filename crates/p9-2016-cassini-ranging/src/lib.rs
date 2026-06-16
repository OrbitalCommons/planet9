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
//! 1. P9 heliocentric position `r_p9(ν)` as a function of its
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

pub mod perturbation;

/// Published reference values from Fienga et al. (2016), kept as labelled
/// constants (NOT used to derive any computed quantity). The favored-ν zone and
/// the Brown & Batygin reference orbit live in
/// [`p9_core::data::ephemeris_constraint`].
pub mod published {
    use p9_core::units::{au, Length};

    /// Approximate P9 heliocentric distance at the preferred true anomaly (AU).
    pub const PREFERRED_HELIO_DISTANCE_AU: f64 = 600.0;

    /// [`PREFERRED_HELIO_DISTANCE_AU`] as a typed [`Length`].
    pub fn preferred_helio_distance() -> Length {
        au(PREFERRED_HELIO_DISTANCE_AU)
    }
}

pub use perturbation::{
    differential_acceleration, favored_true_anomaly, prefit_amplitude,
    range_perturbation_amplitude, range_perturbation_curve, RangePerturbationCurve,
};
