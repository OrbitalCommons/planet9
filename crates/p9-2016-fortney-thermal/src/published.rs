//! Labelled reference constants from Fortney et al. (2016), arXiv:1604.07424.
//!
//! These are the paper's stated values, kept as named constants so the model
//! outputs can be checked against them as residuals. They are NOT used as
//! inputs to the computation (the one model input is the internal-luminosity
//! normalization in the crate root, documented there).

/// Lower edge of the predicted effective-temperature band (K), for the
/// most massive / coldest interior cases at ~4.5 Gyr.
pub const TEFF_MIN_K: f64 = 35.0;

/// Upper edge of the predicted effective-temperature band (K).
pub const TEFF_MAX_K: f64 = 50.0;

/// Lower edge of the mass range the paper explores (M⊕).
pub const MASS_MIN_EARTH: f64 = 5.0;

/// Upper edge of the mass range (M⊕).
pub const MASS_MAX_EARTH: f64 = 20.0;

/// The paper's qualitative detectability conclusion: at 3–5 µm (the WISE/
/// Spitzer near-IR bands) the thermal flux of a cold (~40 K) Planet Nine, while
/// "~20 orders of magnitude larger than blackbody expectations" once a realistic
/// methane-depleted atmosphere is allowed, is still tiny on the Wien tail
/// compared with the far-IR peak. The brightest emission is in the far-IR/mm.
pub const SED_PEAK_REGIME: &str = "far-infrared (tens of microns)";
