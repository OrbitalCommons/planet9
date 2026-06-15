//! Reproduction of Rowan-Robinson (2022),
//! "A search for Planet 9 in the IRAS data" (arXiv:2111.03831).
//!
//! Rowan-Robinson re-analysed the 1983 IRAS far-infrared (60/100 µm)
//! all-sky data, looking for the signature of a slow-moving solar-system
//! body: an unidentified 60 µm point source whose detections appear on some
//! of IRAS's hours-confirmed (HCON) passes but not all, plus an associated
//! single-HCON source from the Reject File offset by a few arcmin (the
//! body's motion between passes). After examining several hundred candidate
//! associations, **one** candidate survives the HCON detection/non-detection
//! pattern. A fitted orbit gives a heliocentric distance of **225 ± 15 AU**
//! and a mass of **3–5 Earth masses** — a *close*, *low-mass* Planet Nine,
//! quite different from the canonical 5–10 M⊕ at 400–800 AU.
//!
//! The candidate is explicitly *tentative*: it is a single survivor out of
//! hundreds of chance associations, near the IRAS detection limit, and the
//! paper calls for dynamical checks and a follow-up optical/NIR search in a
//! 2.5–4° annulus around the 1983 position before any claim of detection.
//!
//! # What this crate computes (no hard-coded answers)
//!
//! From the candidate's far-IR flux and an assumed effective temperature we
//! invert the blackbody + inverse-square law to recover the implied
//! heliocentric distance, reproducing the published ~225 AU. The forward
//! direction (flux from distance) reuses
//! [`p9_2025_iras_akari::thermal_model::P9ThermalParams`] verbatim — the
//! same Planck function, the same Neptune-anchored mass-radius relation from
//! [`p9_core::analysis::photometry::mass_radius_neptunian`] — so the two
//! sibling IRAS crates share one thermal model.
//!
//! - [`distance::implied_distance_au`] inverts `F_ν = π B_ν(T) (R/d)²` for
//!   `d`, given a flux, temperature and mass.
//! - [`distance::CandidateInversion`] bundles the implied distance with the
//!   cross-check flux a P9 at that distance would emit in each IRAS band.
//! - [`chance::chance_alignment_probability`] estimates the probability that
//!   a single-HCON association arises by chance given the IRAS source
//!   surface density.
//!
//! # Reference values
//!
//! The paper publishes the *fitted* distance (225 ± 15 AU) and mass
//! (3–5 M⊕) but not an explicit candidate flux in Jy. The reference flux in
//! [`REF_CANDIDATE`] is therefore the value *self-consistent* with the
//! published distance/mass under this crate's blackbody model at the assumed
//! temperature — it is labelled as such, and the regression tests pin the
//! round-trip (flux → distance → flux) and the published 225 AU, not a flux
//! the paper never stated.

pub mod chance;
pub mod distance;

/// Labelled reference quantities for the surviving IRAS candidate from
/// Rowan-Robinson (2022). The distance and mass are the paper's *fitted*
/// values; the flux and temperature are the self-consistent inputs that
/// reproduce them under the blackbody + inverse-square model (see module
/// docs — the paper does not publish an explicit flux).
#[derive(Debug, Clone, Copy)]
pub struct CandidateReference {
    /// Fitted heliocentric distance (AU). Published: 225 ± 15 AU.
    pub distance_au: f64,
    /// 1-σ distance uncertainty (AU). Published: ±15 AU.
    pub distance_err_au: f64,
    /// Fitted mass midpoint (Earth masses). Published range: 3–5 M⊕.
    pub mass_earth: f64,
    /// Lower / upper fitted mass (Earth masses).
    pub mass_lo_earth: f64,
    pub mass_hi_earth: f64,
    /// Assumed effective temperature (K) of the cold, internally-heated
    /// few-Earth-mass body. Not published; chosen in the cold-ice-giant
    /// range and consistent with the fitted distance.
    pub t_eff_k: f64,
    /// 60 µm flux density (Jy) consistent with the fitted distance/mass at
    /// `t_eff_k`. Back-derived, not published — see module docs.
    pub flux_60um_jy: f64,
    /// 100 µm flux density (Jy), likewise back-derived.
    pub flux_100um_jy: f64,
}

/// The surviving Rowan-Robinson (2022) candidate.
///
/// Distance and mass are the paper's headline fitted values. The flux/
/// temperature pair is the self-consistent input under this crate's model
/// (4 M⊕, 40 K → 225 AU at F₆₀ ≈ 0.35 Jy), labelled as back-derived.
pub const REF_CANDIDATE: CandidateReference = CandidateReference {
    distance_au: 225.0,
    distance_err_au: 15.0,
    mass_earth: 4.0,
    mass_lo_earth: 3.0,
    mass_hi_earth: 5.0,
    t_eff_k: 40.0,
    flux_60um_jy: 0.352,
    flux_100um_jy: 0.858,
};

/// IRAS 60 µm wavelength in metres.
pub const LAMBDA_60UM_M: f64 = 60.0e-6;
/// IRAS 100 µm wavelength in metres.
pub const LAMBDA_100UM_M: f64 = 100.0e-6;
