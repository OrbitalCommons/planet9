//! Labelled reference constants from Kohne & Batygin (2020), kept as named
//! values for cross-checks. These are the paper's CLAIMS, not derived
//! quantities; the computed values in `octupole`, `precession`, and
//! `libration` are pinned against the SCALINGS and qualitative structure these
//! describe, with residuals documented in the crate docs.

/// Planet Nine fiducial parameters used as the perturber in the paper's
/// secular maps (Batygin & Brown 2016 nominal, also used by Kohne & Batygin
/// 2020): a_9 = 700 AU, e_9 = 0.6, m_9 = 10 M_earth.
pub const P9_A_AU: f64 = 700.0;
pub const P9_E: f64 = 0.6;
pub const P9_MASS_EARTH: f64 = 10.0;

/// The octupole-order term is the LEADING apsidal-symmetry-breaking term; its
/// dimensionless strength relative to the quadrupole is
/// epsilon_oct = (15/4) alpha e_9 / (1 - e_9^2). For a clustered ETNO at
/// a ~ 250 AU (alpha ~ 0.36) and e_9 = 0.6, epsilon_oct ~ 1.4 -- i.e. the
/// octupole is NOT a small correction for the relevant orbits, which is the
/// paper's central point (octupole drives the clustering).
pub const PROTOTYPE_ALPHA: f64 = 250.0 / 700.0;

/// Qualitative claim: the quadrupole produces apsidal CIRCULATION, the octupole
/// introduces a dvarpi LIBRATION island (apsidal confinement). Reproduced (not
/// merely asserted) by the `libration` module's circulation-vs-libration
/// classification of the same initial condition under the two truncations.
pub const OCTUPOLE_PRODUCES_LIBRATION: bool = true;
