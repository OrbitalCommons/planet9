//! Reproduction of Pfalzner, Wagner & Bischoff (2026), "Trans-Neptunian
//! Object dynamics even better explained by a stellar flyby after 4.5 Gyr of
//! evolution" (arXiv:2609.03575).
//!
//! # The claim
//!
//! Pfalzner et al. (2024) found that a single stellar flyby — 0.8 M☉,
//! parabolic, periastron 110 AU, inclination 70°, argument of periastron 80°
//! — turns a thin disc into a population resembling today's TNOs
//! immediately after the encounter: cold and hot belts, detached and
//! Sedna-like objects, high-inclination and retrograde orbits. Such an
//! encounter is most likely while the Sun was young, so this paper asks what
//! 4.56 Gyr of subsequent evolution under the giant planets does to that
//! population. Answer: the fit *improves*. The 30 < q < 35 AU excess is
//! cleared (63% of the 30–40 AU tracers are lost, 13% within the first
//! 10 Myr), resonant populations appear, the cold belt loses ~80% and the
//! hot belt ~40% of its members, 7–8% of tracers are injected into the
//! planetary region and 99% of those are ejected, while the detached,
//! Sedna-like and retrograde populations are essentially invariant
//! (Table 1, [`reference`]). The colour dichotomy of the 2025 companion
//! paper survives as well.
//!
//! # What this crate computes
//!
//! * [`flyby`] — the encounter itself: exact parabolic perturber motion,
//!   Dormand–Prince integration of the tracers in the Sun–perturber
//!   barycentric frame, post-flyby heliocentric elements.
//! * [`groups`] — the Table 1 classification and census.
//! * [`evolution`] — the long-term integration with p9-core's hybrid
//!   integrator (Neptune direct, Jupiter–Uranus as a J2/J4 ring), with
//!   per-group retention and regional loss statistics.
//! * [`colours`] — the radial colour gradient and the very-red deficit at
//!   high inclination and eccentricity.
//!
//! The default tests run 4096 tracers through the flyby and 512 of the
//! bound ones for 20 Myr; the 4.56 Gyr paper-scale run is behind
//! `#[ignore]`.

pub mod colours;
pub mod evolution;
pub mod flyby;
pub mod groups;

/// Published values, as labelled targets.
pub mod reference {
    /// Table 1, N_x/N at t = 0 (just after the flyby), Table order:
    /// cold KB, hot KB, detached, Sedna-like, inclined, retrograde, inner.
    pub const TABLE1_FRACTION_T0: [f64; 7] = [0.185, 0.231, 0.093, 0.117, 0.131, 0.035, 0.178];
    /// Table 1, N_x/N at 0.1 Gyr.
    pub const TABLE1_FRACTION_T0_1GYR: [f64; 7] = [0.078, 0.313, 0.122, 0.137, 0.136, 0.034, 0.042];
    /// Table 1, N_x/N at 1 Gyr.
    pub const TABLE1_FRACTION_T1GYR: [f64; 7] = [0.071, 0.290, 0.122, 0.138, 0.136, 0.034, 0.006];
    /// Table 1, N_x/N at 4.5 Gyr.
    pub const TABLE1_FRACTION_T4_5GYR: [f64; 7] = [0.071, 0.237, 0.185, 0.148, 0.019, 0.055, 0.001];
    /// Table 1, N_x/N_init at 4.5 Gyr: retention of each group.
    pub const TABLE1_RETENTION_T4_5GYR: [f64; 7] =
        [0.209, 0.598, 1.078, 1.045, 0.789, 0.833, 0.004];
    /// Loss of the 30 < q < 40 AU tracers within 10 Myr (Section 3.1).
    pub const NEAR_NEPTUNE_LOSS_10MYR: f64 = 0.13;
    /// Loss of the 30 < q < 40 AU tracers by 4.5 Gyr.
    pub const NEAR_NEPTUNE_LOSS_4_5GYR: f64 = 0.63;
    /// Losses in 30 < q < 40 AU by eccentricity: e < 0.2 vs e > 0.4.
    pub const LOW_E_LOSS: f64 = 0.80;
    pub const HIGH_E_LOSS: f64 = 0.47;
    /// Losses in 30 < q < 40 AU by inclination: low-i vs high-i.
    pub const LOW_I_LOSS: f64 = 0.51;
    pub const HIGH_I_LOSS: f64 = 0.21;
    /// Fraction of tracers injected into the planetary region, and the
    /// fraction of those subsequently ejected (Section 3.2.5).
    pub const INJECTED_FRACTION: (f64, f64) = (0.07, 0.08);
    pub const INJECTED_EJECTED: f64 = 0.991;
    /// Disc mass captured by the perturber, model A1 (Pfalzner et al. 2024).
    pub const CAPTURED_BY_PERTURBER: f64 = 0.083;
    /// Loss rate in the first 10–100 Myr relative to the present one.
    pub const EARLY_LOSS_RATE_RATIO: f64 = 100.0;
}
