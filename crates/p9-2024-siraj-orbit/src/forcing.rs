//! Secular forcing strength of an exterior perturber on the clustered ETNOs.
//!
//! For an exterior perturber of mass `m` and semi-major axis `a_p` acting on
//! a test body of fixed semi-major axis `a`, the leading secular coupling
//! constant (the rate at which the perturber forces the test body's apsidal
//! line) scales as
//!
//!   ω̇_sec ∝ n (m_p/M☉) (a/a_p)³ ∝ m_p a^{3/2} / a_p³   (quadrupole,
//!   α = a/a_p ≪ 1 limit)
//!
//! At fixed test-particle semi-major axis `a`, the *mass-and-distance*
//! dependence of the forcing that the clustered population responds to is
//! commonly written in the reduced "secular strength" form
//!
//!   S(m, a_p) = m / a_p³,
//!
//! the leading m, a_p coupling once the (slowly varying, sample-fixed) a²/a_p²
//! geometric prefactor is absorbed into the normalization. This is the
//! scaling Siraj, Chyba & Tremaine (2024) use to translate the *observed*
//! apsidal-confinement strength into a perturber mass-distance relation, and
//! it produces a degeneracy ridge m ∝ a_p³ — a more distant perturber needs a
//! larger mass to force the same confinement.
//!
//! `m` is carried in Earth masses, `a_p` in AU. `S` therefore has units of
//! M⊕ / AU³; only ratios of `S` enter the likelihood, so the normalization is
//! immaterial.

/// Reduced secular forcing strength S(m, a_p) = m / a_p³.
///
/// `mass_earth`: perturber mass in Earth masses.
/// `a_p`: perturber semi-major axis in AU.
pub fn forcing_strength(mass_earth: f64, a_p: f64) -> f64 {
    mass_earth / a_p.powi(3)
}

/// The mass (Earth masses) that places the degeneracy ridge through a given
/// `a_p` at fixed forcing `s`: inverts S = m / a_p³ ⇒ m = S · a_p³.
pub fn mass_for_strength(s: f64, a_p: f64) -> f64 {
    s * a_p.powi(3)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_forcing_inverse_consistent() {
        let s = forcing_strength(4.4, 290.0);
        assert_relative_eq!(mass_for_strength(s, 290.0), 4.4, epsilon = 1e-9);
    }

    #[test]
    fn test_ridge_is_cubic_in_a() {
        // Doubling a_p along a fixed-S ridge requires 8x the mass.
        let s = forcing_strength(4.4, 290.0);
        let m1 = mass_for_strength(s, 290.0);
        let m2 = mass_for_strength(s, 580.0);
        assert_relative_eq!(m2 / m1, 8.0, epsilon = 1e-9);
    }
}
