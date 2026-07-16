//! Reproduction of the 2016 resonance-*prediction* papers:
//!
//!   - Malhotra, Volk & Wang (2016), "Corralling a Distant Planet with
//!     Extreme Resonant Kuiper Belt Objects" (arXiv:1603.02196), and
//!   - Millholland & Laughlin (2016), "Constraints on Planet Nine's Orbit and
//!     Sky Position within a Framework of Mean-Motion Resonances"
//!     (arXiv:1612.07774).
//!
//! HEADLINE. If the longest-period clustered ETNOs are each trapped in a
//! small-integer (N:1 or N:2) mean-motion resonance with Planet Nine, then
//! their observed orbital periods are commensurate with P9's period and
//! therefore *predict* it. An object in the exterior p:q resonance sits at
//!
//!   a_ETNO = a9 · (q/p)^(2/3)        ⇔        P9 / P_ETNO = p / q,
//!
//! i.e. the interior ETNO completes p orbits while P9 completes q (p > q).
//! Scanning a9 and asking which value places the whole distant set closest
//! to low-order resonances localizes Planet Nine.
//!
//! Malhotra et al. report a9 ≈ 665 AU (period ≈ 17,117 yr), with the four
//! longest-period ETNOs (Sedna, 2004 VN112, 2010 GB174, 2012 VP113) near
//! period ratios of small integers (3:2, 5:2, 3:1, 4:1 — the N:1 / N:2
//! family the DEFAULT_Q_MAX = 2 scan admits). Millholland &
//! Laughlin independently localize P9 near a9 ≈ 654 AU.
//!
//! HONEST CAVEAT. These resonance predictions are NOT confirmed. The sibling
//! crate `p9-2018-resonance` reproduces Bailey, Brown & Batygin (2018), which
//! shows that the space of high-order resonances is so dense that *any* a9
//! can be made "commensurate" with a handful of objects; the probability that
//! all of the observed ETNOs genuinely sit in N:1 or N:2 resonances is < 5%.
//! We therefore report the best-fit a9 together with its residuals and do not
//! treat the prediction as a measurement.
//!
//! This crate reuses `p9_core::analysis::resonance::resonance_semi_major_axis`
//! for the resonance-location arithmetic and the vetted
//! `p9_core::data::etno::BROWN_2017_SAMPLE` for the observational sample.

pub mod prediction;

pub use prediction::{
    best_fit_a9, longest_period_etnos, malhotra_constraining_set, nearest_resonance,
    orbital_period_years, period_ratio, resonance_residual, scan_a9, score_a9, A9Score,
    ResonanceMatch, DEFAULT_P_MAX, DEFAULT_Q_MAX, MALHOTRA_2016_A9_AU, MALHOTRA_2016_PERIOD_YR,
    MILLHOLLAND_2016_A9_AU,
};
