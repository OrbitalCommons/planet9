//! Reproduction of Karim Benakli (2026), "A Solar-System Window for Hidden
//! Stellar Companions" (arXiv:2606.25190).
//!
//! The paper asks whether the *nearest* stellar / sub-stellar object could be a
//! hidden-sector ("hidden-brane") companion bound to the Sun at hundreds to a
//! few thousand AU, rather than an ordinary star at ~1 pc. It builds a
//! phenomenological mass–distance map from an illustrative ephemeris tidal
//! envelope calibrated to Planet-Nine-like dynamical constraints, and reaches
//! two conclusions reproduced here.
//!
//! # 1. The tidal mass–distance envelope ([`envelope`])
//!
//! A distant perturber of mass `M` at distance `d` raises a tidal/quadrupole
//! acceleration on the inner Solar System that scales as `G·M/d³`. Holding that
//! amplitude fixed at the level the ephemerides still permit gives a maximum
//! allowed mass that grows as the cube of distance:
//!
//! ```text
//!     M_max(d) = M_anchor · (d / d_anchor)³
//! ```
//!
//! Calibrating the anchor to the Planet-Nine dynamical point —
//! [`ANCHOR_MASS_EARTH`] = 5 M⊕ at [`ANCHOR_DISTANCE_AU`] = 500 AU, the Brown &
//! Batygin (2019) revised best fit carried in
//! [`p9_core::types::P9Params::revised_2019`] — the envelope permits:
//!
//! | distance | M_max          | regime           |
//! |----------|----------------|------------------|
//! | 300 AU   | ≈ 1.1 M⊕       | Earth-scale      |
//! | 1000 AU  | ≈ 40 M⊕        | sub-Saturn       |
//! | 2000 AU  | ≈ 320 M⊕ ≈ 1 M_J | Jovian         |
//!
//! exactly the Earth → sub-Saturn → Jovian window the paper documents.
//!
//! # 2. A smooth dark-matter halo cannot do it ([`halo`])
//!
//! If the companion were instead a smooth galactic dark-matter over-density,
//! the enclosed mass is `M_DM(<r) = ρ_local·(4/3)π r³`. With the local density
//! ρ_local ≈ 0.01 M⊙/pc³ ([`RHO_LOCAL_DM_MSUN_PC3`], ≈ 0.4 GeV/cm³), the mass
//! within 1000 AU is *sub-Pluto* (≈ 0.73 [`M_PLUTO_KG`]). So a smooth halo
//! cannot supply a nearby planet-scale object; any such companion must be a
//! structured, gravitationally bound body.
//!
//! # p9-core reuse
//!
//! Unit conversions reuse [`p9_core::constants::PC_AU`] (parsec in AU) and
//! [`p9_core::units::SOLAR_MASS_KG`] / [`p9_core::units::EARTH_MASS_KG`] rather
//! than redefining them; the Planet-Nine anchor mirrors
//! [`p9_core::types::P9Params::revised_2019`].

pub mod envelope;
pub mod halo;

pub use envelope::mass_envelope_earth;
pub use halo::{dm_halo_mass_within_kg, dm_halo_mass_within_solar};

/// Planet-Nine-like anchor mass for the tidal envelope, in Earth masses.
///
/// Brown & Batygin (2019) revised best fit (see
/// [`p9_core::types::P9Params::revised_2019`], `mass_earth = 5.0`).
pub const ANCHOR_MASS_EARTH: f64 = 5.0;

/// Planet-Nine-like anchor distance for the tidal envelope, in AU.
///
/// Brown & Batygin (2019) revised semi-major axis (see
/// [`p9_core::types::P9Params::revised_2019`], `a = 500.0`).
pub const ANCHOR_DISTANCE_AU: f64 = 500.0;

/// Local dark-matter density in solar masses per cubic parsec.
///
/// ρ_local ≈ 0.01 M⊙/pc³ ≈ 0.4 GeV/cm³ (standard local DM density; e.g.
/// Read 2014, J. Phys. G 41, 063101). Used by the smooth-halo bound in
/// [`halo`].
pub const RHO_LOCAL_DM_MSUN_PC3: f64 = 0.01;

/// Pluto's mass in kilograms (1.303e22 kg; Stern et al. 2015, New Horizons).
/// The smooth dark-matter mass within 1000 AU is compared against this.
pub const M_PLUTO_KG: f64 = 1.303e22;
