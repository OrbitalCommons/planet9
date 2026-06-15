//! Reproduction of Oldroyd & Trujillo (2021), "The Perihelion Gap: An
//! Observational Signature of Distant Solar System Architecture"
//! (arXiv:2107.06296).
//!
//! Headline: the perihelion (q) distribution of the distant minor planets
//! shows a relative DEFICIT — a *gap* — at intermediate perihelia, roughly
//! q ≈ 50–65 AU at high eccentricity (e ≳ 0.65). The gap separates two
//! distinct distant populations:
//!
//!   - the Neptune-coupled SCATTERING extreme TNOs (ETNOs) at low perihelion
//!     (q ≳ 40 AU, a ≳ 150 AU), whose perihelia are pumped up the outer edge
//!     of the Neptune-scattering region; and
//!   - the DETACHED inner-Oort-cloud (IOC) / Sednoid population at high
//!     perihelion (q ≳ 65 AU, a ≳ 250 AU), dynamically decoupled from Neptune.
//!
//! Oldroyd & Trujillo argue the gap is very unlikely to arise from a single
//! continuous distribution and is instead the boundary between these two
//! sculpting mechanisms — a feature a distant planet could help carve.
//!
//! What this crate computes, all from real (seeded) calculations reusing
//! `p9_core`, with no hard-coded answers:
//!
//! - `sample`: the observed distant-TNO perihelion sample. The vetted
//!   `p9_core::data::etno::BROWN_2017_SAMPLE` (10 objects, a ≳ 230 AU) plus a
//!   documented extended distant-TNO table covering the lower-a ETNOs and the
//!   high-q Sednoids needed to see both flanks of the gap.
//! - `distribution`: histogramming of perihelia and the dip statistic — the
//!   ratio of the count density inside the gap to the density of the two
//!   flanking populations.
//! - `synthetic`: a seeded synthetic two-component (scattering + detached)
//!   population that reproduces a gap, contrasted with a single continuous
//!   control population that does not.
//! - `boundaries`: relates the gap edges to physics from p9-core: the high-q
//!   scattering edge to the Neptune resonance-overlap boundary
//!   (`p9_core::analysis::resonance::critical_perihelion`) and the low-q
//!   detached edge to the Sednoid onset.
//!
//! Documented model reductions (see `REPRODUCTION_NOTES.md`): the gap is
//! partly observational — the published feature folds in survey detection
//! biases against faint, distant, high-q objects that this reduced model does
//! not invert. We reproduce the *gap statistic* (a real deficit at q ≈ 50–65
//! AU in the vetted sample, and in a seeded two-population draw) and its
//! location, not the full debiased completeness argument.

pub mod boundaries;
pub mod distribution;
pub mod sample;
pub mod synthetic;

/// Published reference values from Oldroyd & Trujillo (2021), kept as labelled
/// constants for cross-checks (never asserted as if computed).
pub mod published {
    /// Low-perihelion edge of the gap (AU). Below this, the Neptune-coupled
    /// scattering ETNO population dominates. Abstract: the gap runs "between
    /// roughly 50 and 65 au".
    pub const GAP_Q_LOW_AU: f64 = 50.0;

    /// High-perihelion edge of the gap (AU). Above this, the detached
    /// inner-Oort-cloud (IOC) / Sednoid population dominates (q ≳ 65 AU).
    pub const GAP_Q_HIGH_AU: f64 = 65.0;

    /// Gap centre (AU), the midpoint of the published range — the labelled
    /// location the computed deficit is checked against.
    pub const GAP_Q_CENTER_AU: f64 = 0.5 * (GAP_Q_LOW_AU + GAP_Q_HIGH_AU);

    /// Eccentricity threshold at which the gap is clearest (e ≳ 0.65): the
    /// distant high-e orbits whose perihelia map this q range.
    pub const GAP_ECCENTRICITY_FLOOR: f64 = 0.65;

    /// Semi-major-axis floor of the ETNO sample (a ≳ 150 AU).
    pub const ETNO_A_FLOOR_AU: f64 = 150.0;

    /// Perihelion floor of the ETNO sample (q ≳ 40 AU).
    pub const ETNO_Q_FLOOR_AU: f64 = 40.0;

    /// Perihelion floor of the inner-Oort-cloud / detached population
    /// (q ≳ 65 AU).
    pub const IOC_Q_FLOOR_AU: f64 = 65.0;

    /// Predicted relative abundance of objects *inside* the gap compared with
    /// the nearby IOC population (65 ≲ q ≲ 100 AU): "only ∼20% as numerous".
    pub const GAP_RELATIVE_ABUNDANCE: f64 = 0.20;
}
