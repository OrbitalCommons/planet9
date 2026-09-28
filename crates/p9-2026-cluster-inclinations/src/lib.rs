//! Reproduction of Bansal, Brunton, Batygin, Pichierri & Nesvorný (2026),
//! "Distant TNO Inclinations as a Constraint on Primordial Cluster
//! Perturbations in the Presence of Planet Nine" (arXiv:2607.15646).
//!
//! # The question
//!
//! The distant (q > 40 AU, a > 200 AU) trans-Neptunian population has a
//! modest inclination dispersion. Hu et al. (2025) read that as evidence
//! against strong stellar perturbations in the Sun's birth cluster, while
//! Nesvorný et al. (2023) need strong cluster perturbations to explain the
//! radial extent of the scattered disk. Could Planet Nine reconcile the two
//! by dynamically *cooling* a cluster-excited population over 4 Gyr into
//! something that looks cold today?
//!
//! # The paper's answer
//!
//! No. Integrating cluster-excited and quiescent populations under Neptune,
//! a solar J2 for the inner giants, Planet Nine (m₉ ∈ {5, 7.07, 10} M⊕,
//! i₉ = 20°, q₉ ≈ 250–300 AU), the Galactic tide and passing stars for
//! 4 Gyr:
//!
//! * the observed 19 high-q TNOs have an intrinsic (Brown 2001 debiased)
//!   width `w_obs = 12° (+6/−5)` — [`REFERENCE`];
//! * cluster-influenced populations stay at `w ≈ 26–27.5°`, rejected at
//!   ~3σ (p = 0.0014) by the observed sample;
//! * cluster-free populations initialised at σ = 15° end at `w ≈ 16–19.5°`,
//!   within 1.5σ of the observations — Planet Nine neither cools nor
//!   appreciably heats the distant population;
//! * cluster-free runs reproduce strong perihelion clustering (von Mises
//!   κ > 1), cluster-influenced runs only weak clustering (κ < 1).
//!
//! # What this crate computes
//!
//! * [`sample`] — the 19-object observed sample (JPL SBDB elements) with
//!   discovery latitudes recomputed from the orbits.
//! * [`debias`] — the Brown (2001) latitude-conditional CDF, Kuiper-test scan
//!   over `w`, and Monte-Carlo calibrated confidence intervals; the
//!   `w_obs ≈ 12°` headline and the 3σ rejection of `w = 26°` are computed.
//! * [`population`] — cluster-influenced and cluster-free initial
//!   conditions (the latter is p9-core's scattered-disk generator).
//! * [`stars`] — impulse-approximation passing stars.
//! * [`simulation`] — the WHM integration (Planet Nine, averaged giant
//!   planets, Galactic tide, passing stars), parallel over particle chunks,
//!   at reduced (100 Myr) and paper (4 Gyr) scale.
//! * [`width`] / [`clustering`] — the paper's `w` and `κ` estimators.
//!
//! The reduced-scale integration is shorter than the paper's by 10×; what it
//! tests is the paper's mechanism — Planet Nine's secular forcing broadens
//! a cold population toward the 16–19.5° band and leaves an isotropically
//! excited one at w ≳ 26°, never cooling it — rather than the exact Table 1
//! values, which the `#[ignore]`d paper-scale test targets.
//!
//! A longer check with the `long_run` example (1024 cluster-influenced
//! particles, m₉ = 5 M⊕, a₉ = 367 AU, e₉ = 0.2, 2 Gyr, run 2026-09-26): the
//! windowed width wanders between 19° and 30° as orbits cycle through the
//! q ∈ (40, 80) AU window (a few dozen at a time), pooling to 25.6° over the
//! final 500 Myr, the whole-population width holds at 24–29°, 80% of the
//! particles survive, and κ = 0.78 — the paper's cluster-influenced row
//! (w ≈ 26–27°, κ < 1) at reduced scale. The same run with perihelia drawn
//! up to 300 AU (the excited population extends far beyond the window)
//! pools to 24.8° with κ = 0.26 and 87% survival. Note that the paper's own
//! `sin i · exp(−i²/2w²)` form is a poor description of a flat excited
//! distribution (its Figure 2), so the windowed fit is noisy at small n:
//! the default 192-particle test sees only 20–50 orbits in the window and
//! its windowed width drifts down to ~17° over 400 Myr while the whole
//! population holds at 25–34°.

pub mod clustering;
pub mod debias;
pub mod population;
pub mod sample;
pub mod simulation;
pub mod stars;
pub mod width;

/// Published values from the paper, kept as labelled targets for the tests.
pub mod reference {
    /// Intrinsic width of the observed high-q sample (deg), Section 3.1.
    pub const W_OBS_DEG: f64 = 12.0;
    /// Its 1σ interval (deg): 12 (+6/−5).
    pub const W_OBS_1SIGMA_DEG: (f64, f64) = (7.0, 18.0);
    /// Cluster-influenced simulation widths span this band (deg), Table 1.
    pub const W_CLUSTER_INFLUENCED_DEG: (f64, f64) = (26.0, 27.5);
    /// Cluster-free simulation widths span this band (deg), Table 1.
    pub const W_CLUSTER_FREE_DEG: (f64, f64) = (16.0, 19.5);
    /// Primordial half-Gaussian σ of the cluster-free runs (deg). In the
    /// paper's `sin i · exp(−i²/2w²)` measure this population starts at
    /// `w ≈ 10°` (see [`crate::population`]).
    pub const SIGMA_CLUSTER_FREE_INITIAL_DEG: f64 = 15.0;
    /// Initial fitted width of a σ = 15° half-Gaussian population (deg).
    pub const W_CLUSTER_FREE_INITIAL_DEG: f64 = 10.0;
    /// Monte-Carlo p-value at which w = 26° is rejected (~3σ).
    pub const P_REJECT_W26: f64 = 0.0014;
    /// Number of observed high-q TNOs.
    pub const N_OBSERVED: usize = 19;
    /// Table 1 κ: cluster-influenced runs stay below, cluster-free above.
    pub const KAPPA_THRESHOLD: f64 = 1.0;
}
