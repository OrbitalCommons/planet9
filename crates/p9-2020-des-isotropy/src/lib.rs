//! Reproduction of Bernardinelli et al. (2020),
//! "Testing the isotropy of the Dark Energy Survey's extreme trans-Neptunian
//! objects" (arXiv:2003.08901).
//!
//! Headline: applying a battery of isotropy test statistics to the DES Y4
//! eTNO detections, the DES sample on its own is *consistent with an isotropic
//! (uniform) distribution* of orbital angles — DES alone shows no
//! statistically significant clustering in longitude of perihelion (ϖ),
//! argument of perihelion (ω), or longitude of ascending node (Ω). The paper
//! ran Kuiper's test on 3 angles × 4 sample definitions (12 tests); only 2
//! reached p = 0.03, both for Ω at a > 250 AU.
//!
//! This crate:
//!   * encodes the DES eTNO table from the paper (`des_sample`),
//!   * implements a suite of ≥5 isotropy statistics, reusing the
//!     `p9_core::analysis::circular` primitives (Rayleigh, Kuiper, R̄) and
//!     adding Watson U², Beran/Ajne and a two-angle circular correlation
//!     (`statistics`),
//!   * runs the full Case × angle battery with seeded Monte Carlo nulls and
//!     summarizes consistency with isotropy (`analysis`).
//!
//! Two nulls are provided (`selection::NullModel`). Against a *flat* isotropic
//! null the encoded DES angles look clustered (Ω has R̄ ≈ 0.8) — this is the
//! "clustering vs isotropy" framing the paper reacts against. The paper's
//! actual null is an isotropic population folded through the DES selection
//! function; we reproduce that with a documented static-footprint acceptance
//! (`selection::des_footprint_weight`). Under the selection null the longitude
//! of perihelion ϖ — the Planet 9 alignment angle — is consistent with
//! isotropy in every Case, reproducing the paper's headline.
//!
//! Residual: the static RA/Dec footprint cannot fully explain away the very
//! high Ω concentration (R̄ up to 0.99 in the n = 3–4 sub-samples); several Ω
//! tests stay significant. This matches the paper's finding that Ω was the
//! only angle to trip any test (2 of 12 Kuiper tests, both Ω at a > 250 AU),
//! but the paper's per-exposure cadence × depth × linking model (reproduced in
//! `p9-2022-des`) removes more of the Ω effect than our static footprint. We
//! pin isotropy for ϖ/ω and pin Ω as the labelled residual rather than
//! forcing it to isotropy. See `analysis::PAPER_CONCLUSION` for the reference
//! statement.

pub mod analysis;
pub mod des_sample;
pub mod selection;
pub mod statistics;
