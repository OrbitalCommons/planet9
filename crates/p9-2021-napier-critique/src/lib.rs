//! Reproduction of Napier et al. (2021)
//! "No Evidence for Orbital Clustering in the Extreme Trans-Neptunian
//! Objects" (PSJ 2, 59; arXiv:2102.05601).
//!
//! Headline result reproduced: the apparent clustering in the longitude of
//! perihelion ϖ of the ETNOs is a *selection artifact*, not evidence for
//! Planet Nine. The naive Rayleigh test — which compares the sample to a
//! pure-uniform null — reports significance, but that null is wrong: the
//! OSSOS + DES + S&T surveys cover the sky incompletely and to finite
//! depth, so the detected ϖ distribution is the intrinsic distribution
//! *times* a selection function. When the detections are debiased through
//! the composite survey selection function, the observed clustering is
//! *consistent with an isotropic intrinsic distribution* — the consistency
//! p-value lands in the published ~17%–94% range.
//!
//! Pipeline (`critique`):
//!   - ϖ sample from the vetted `p9_core::data::etno::BROWN_2017_SAMPLE`.
//!   - composite first+second-harmonic selection weight w(ϖ) for the
//!     OSSOS+DES+S&T union (amplitude/phase documented in
//!     `critique::SelectionFunction`).
//!   - seeded `StdRng` bias-aware Monte Carlo: consistency p that the
//!     observed mean resultant length R̄ arises under isotropy × selection.
//!   - the naive uniform-null Rayleigh p (`p9_core::analysis::circular`)
//!     for contrast.
//!
//! What is pinned in tests:
//!   1. naive Rayleigh p < 0.05 (clustered under the wrong null);
//!   2. debiased consistency p ∈ [0.17, 0.94] (consistent with isotropy);
//!   3. the debiased p ≫ the naive p (debiasing flips the conclusion).
//! The published 17%–94% band is kept as the labelled reference constant
//! `critique::NAPIER_2021_CONSISTENCY_BAND`.
//!
//! Residual: the selection function is a smooth two-harmonic stand-in for
//! the surveys' full pointing/depth histories (documented in
//! `SelectionFunction`), and the ϖ sample is the 10-object vetted
//! `BROWN_2017_SAMPLE` rather than Napier's exact OSSOS-debiased set. The
//! qualitative reproduction — the conclusion flips from "clustered" to
//! "consistent with isotropy" under debiasing — is robust; the precise
//! consistency p depends on the assumed amplitudes/phases.

pub mod critique;
