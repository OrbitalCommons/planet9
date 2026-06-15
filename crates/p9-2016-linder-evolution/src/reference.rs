//! Published reference magnitudes from Linder & Mordasini (2016),
//! arXiv:1602.07465. These are the paper's *predictions* — labelled constants
//! used only as targets in the regression tests, never as computed answers.
//!
//! Nominal case: a 10 M⊕ Planet Nine at ~700 AU (heliocentric), effective
//! temperature ≈ 47 K. The optical/near-IR magnitudes are reflected sunlight;
//! the L/N/Q and WISE W1–W4 magnitudes are intrinsic thermal radiation.

/// Effective temperature of the nominal 10 M⊕ Planet Nine at 700 AU (K).
pub const T_EFF_10ME_700AU_K: f64 = 47.0;

/// Reflected-light V magnitude, 10 M⊕ at 700 AU.
pub const V_10ME_700AU: f64 = 21.7;
/// Reflected-light R magnitude, 10 M⊕ at 700 AU.
pub const R_10ME_700AU: f64 = 21.4;
/// Reflected-light I magnitude, 10 M⊕ at 700 AU.
pub const I_10ME_700AU: f64 = 21.0;

/// Thermal WISE W1 (3.4 µm) magnitude, 10 M⊕ at 700 AU.
pub const W1_10ME_700AU: f64 = 20.1;
/// Thermal WISE W2 (4.6 µm) magnitude, 10 M⊕ at 700 AU.
pub const W2_10ME_700AU: f64 = 20.1;
/// Thermal WISE W3 (12 µm) magnitude, 10 M⊕ at 700 AU.
pub const W3_10ME_700AU: f64 = 18.6;
/// Thermal WISE W4 (22 µm) magnitude, 10 M⊕ at 700 AU.
pub const W4_10ME_700AU: f64 = 10.2;

/// Thermal Q-band (~20 µm) magnitude, 10 M⊕ at 700 AU.
pub const Q_10ME_700AU: f64 = 10.7;
/// Thermal N-band (~10 µm) magnitude, 10 M⊕ at 700 AU.
pub const N_10ME_700AU: f64 = 19.9;
/// Thermal L-band (~3.5 µm) magnitude, 10 M⊕ at 700 AU.
pub const L_10ME_700AU: f64 = 20.1;

/// Reflected V magnitude at aphelion vs mass (Table in the paper): mass in M⊕,
/// V magnitude. The body brightens slowly with mass (R ∝ M^~0.27).
pub const V_VS_MASS_APHELION: [(f64, f64); 4] =
    [(5.0, 24.3), (10.0, 23.7), (20.0, 23.3), (50.0, 22.6)];

/// Intrinsic Q-band magnitude at aphelion vs mass: mass in M⊕, Q magnitude.
/// The thermal flux climbs steeply with mass (hotter, larger interior).
pub const Q_VS_MASS_APHELION: [(f64, f64); 4] =
    [(5.0, 14.6), (10.0, 11.7), (20.0, 9.2), (50.0, 5.8)];
