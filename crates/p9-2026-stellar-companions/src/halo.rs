//! Smooth dark-matter halo bound on a hidden companion.
//!
//! If the nearby companion were instead a smooth over-density of galactic
//! dark matter, the most mass that could be assembled inside a radius `r` is
//! the local DM density integrated over the enclosed volume:
//!
//! ```text
//!     M_DM(<r) = ρ_local · (4/3)·π·r³
//! ```
//!
//! With the canonical local density ρ_local ≈ 0.01 M⊙/pc³ (≈ 0.4 GeV/cm³),
//! the mass inside 1000 AU is ≈ 9.5e21 kg ≈ 0.73 Pluto masses — sub-Pluto.
//! A smooth halo therefore cannot supply a nearby planet-scale object: any
//! such object must be a *structured*, gravitationally bound body (the paper's
//! "hidden-brane" companion).
//!
//! Unit chain: ρ_local is stored in M⊙/pc³; we convert the enclosed volume
//! (built in AU³) to pc³ via [`PC_AU`](p9_core::constants::PC_AU), multiply by
//! ρ_local to get solar masses, then to kilograms via
//! [`SOLAR_MASS_KG`](p9_core::units::SOLAR_MASS_KG).

use crate::RHO_LOCAL_DM_MSUN_PC3;
use p9_core::constants::PC_AU;
use p9_core::units::SOLAR_MASS_KG;
use std::f64::consts::PI;

/// Dark-matter mass (solar masses) enclosed within radius `r_au` (AU) for a
/// smooth halo of local density [`RHO_LOCAL_DM_MSUN_PC3`].
///
/// `M_DM(<r) = ρ_local · (4/3)·π·r³`.
pub fn dm_halo_mass_within_solar(r_au: f64) -> f64 {
    let r_pc = r_au / PC_AU;
    let volume_pc3 = (4.0 / 3.0) * PI * r_pc * r_pc * r_pc;
    RHO_LOCAL_DM_MSUN_PC3 * volume_pc3
}

/// Same enclosed dark-matter mass, expressed in kilograms.
pub fn dm_halo_mass_within_kg(r_au: f64) -> f64 {
    dm_halo_mass_within_solar(r_au) * SOLAR_MASS_KG
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::M_PLUTO_KG;
    use approx::assert_relative_eq;

    #[test]
    fn scales_as_radius_cubed() {
        let r = dm_halo_mass_within_solar(2000.0) / dm_halo_mass_within_solar(1000.0);
        assert_relative_eq!(r, 8.0, epsilon = 1e-9);
    }

    #[test]
    fn within_1000_au_is_sub_pluto() {
        let m = dm_halo_mass_within_kg(1000.0);
        assert!(
            m < M_PLUTO_KG,
            "smooth DM within 1000 AU ({m:e} kg) must be below Pluto ({M_PLUTO_KG:e} kg)"
        );
        // Sub-Pluto but the same order of magnitude (~0.73 M_Pluto).
        assert!(m > 0.5 * M_PLUTO_KG);
    }
}
