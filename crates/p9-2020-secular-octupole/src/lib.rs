//! Reproduction of Kohne & Batygin (2020),
//! "Secular dynamics of distant Kuiper belt objects under Planet Nine"
//! (octupole-level secular theory).
//!
//! # Headline reproduced
//!
//! Extending the orbit-averaged secular theory of the ETNO-Planet Nine
//! interaction from QUADRUPOLE to OCTUPOLE order reproduces the apsidal
//! confinement (libration of `dvarpi = varpi - varpi_9`) that the quadrupole
//! cannot. At quadrupole order the coplanar doubly-averaged Hamiltonian is a
//! function of eccentricity ALONE, so the apsidal angle circulates uniformly
//! and there is no aligned/anti-aligned structure. The octupole term -- the
//! leading apsidal-symmetry-breaking term, scaling with the perturber
//! eccentricity `e_9` and with `alpha = a / a_9` -- introduces a `dvarpi`
//! libration island. For the relevant ETNO orbits (a ~ 250 AU, alpha ~ 0.36,
//! e_9 = 0.6) the dimensionless octupole strength
//! epsilon_oct = (15/4) alpha e_9 / (1 - e_9^2) is of order unity, so the
//! octupole is NOT a small correction: it is what drives the clustering.
//!
//! # What is computed (no hard-coded answers)
//!
//! - [`octupole`]: the coplanar secular Hamiltonian to quadrupole and octupole
//!   order. The quadrupole piece REUSES
//!   `p9_core::analysis::secular::coplanar_quadrupole`; the octupole piece adds
//!   the `+C_oct e e_9 cos(dvarpi)` term (positive sign, matching the exact
//!   ring average) with
//!   `C_oct = K_OCT (GM_9/a_9) alpha^3/(1-e_9^2)^{5/2}`, K_OCT = 1.05
//!   calibrated against p9-core's numerical ring average. The
//!   octupole-to-quadrupole strength ratio `epsilon_oct` and its scaling with
//!   `e_9` and `alpha` are pinned. The analytic octupole amplitude is
//!   cross-checked at small `alpha` against the cos(dvarpi) Fourier amplitude
//!   of `p9_core`'s EXACT numerical ring average
//!   (`numerical_secular_hamiltonian`), which contains the octupole and all
//!   higher multipoles automatically.
//!
//! - [`precession`]: the apsidal precession rate `d(dvarpi)/dt = dH/dGamma`.
//!   The quadrupole rate is `dvarpi`-independent (circulation); the octupole
//!   modulates it with the apsidal angle. The closed-form octupole modulation
//!   amplitude is pinned against the finite-difference aligned-vs-anti-aligned
//!   precession-rate spread.
//!
//! - [`libration`]: integrates the one-degree-of-freedom secular EOM in
//!   (Gamma, dvarpi) and classifies the apsidal motion as circulation or
//!   libration. The headline test shows the SAME initial orbit circulates
//!   under quadrupole-only but librates once the octupole is included.
//!
//! # Pinned results (see module tests)
//!
//! - Quadrupole coplanar Hamiltonian has no `dvarpi` dependence -> circulation.
//! - Octupole strength `epsilon_oct` scales LINEARLY with `e_9` (small e_9) and
//!   with `alpha`; the octupole coupling vanishes for a circular perturber.
//! - The exact ring-average cos(dvarpi) amplitude grows with both `e_9` and
//!   `alpha`, matching the analytic octupole scaling.
//! - The analytic octupole amplitude agrees with the ring-average amplitude at
//!   small `alpha` (convergent regime) within 5%.
//! - A test orbit that CIRCULATES under quadrupole-only LIBRATES once the
//!   octupole term is added (the apsidal-confinement island), with libration
//!   center at dvarpi = 0 and e* = epsilon_oct/3.
//!
//! # Residuals / honest caveats
//!
//! - **Quantitative validity of the analytic octupole.** The multipole
//!   expansion converges only for `alpha << 1`. At the `alpha ~ 0.3-0.9`
//!   relevant to Planet Nine it diverges, so the analytic `C_oct` term is only
//!   quantitatively accurate at small `alpha` (where it is cross-checked). For
//!   the prototype ETNO orbits the trustworthy reference is `p9_core`'s exact
//!   numerical ring average; the analytic SCALINGS (in `e_9`, `alpha`) and the
//!   qualitative libration-island structure are what the crate asserts at large
//!   `alpha`, not the analytic amplitude itself.
//!
//! - **Aligned vs anti-aligned libration center.** The coplanar octupole
//!   libration center comes out APSIDALLY ALIGNED (dvarpi = 0, e* =
//!   epsilon_oct/3); the observed ETNOs are ANTI-aligned with Planet Nine. The
//!   anti-alignment sign requires the inclined / perihelion-detached / resonant
//!   dynamics, not the coplanar secular fixed point. The crate reproduces the
//!   octupole-driven libration MECHANISM and the e*, epsilon_oct scalings, and
//!   documents the sign as out of reach of the coplanar reduction (the same
//!   honest caveat as p9-2019-selfgrav-disk).
//!
//! - **Eccentricity-inclination coupling.** The paper emphasizes that the
//!   octupole couples eccentricity and inclination (absent at the coplanar
//!   quadrupole level used here for the apsidal map). The libration
//!   demonstration is done in the coplanar limit, where the apsidal island is
//!   cleanest; the inclined e-i coupling is a property of the exact 3-D ring
//!   average (`p9_core`'s `numerical_secular_hamiltonian`, which is fully
//!   inclined) rather than asserted here as a separate scalar, to avoid
//!   over-claiming from the reduced coplanar model.
//!
//! - **Published scalars.** Kohne & Batygin (2020)'s specific libration
//!   amplitudes and periods depend on the full inclined map and the assumed P9
//!   orbit; they are kept as labelled reference constants in [`published`]
//!   rather than asserted, since the coplanar reduction here reproduces the
//!   MECHANISM (octupole -> libration island, octupole strength ~ e_9 alpha)
//!   not the exact island geometry.

pub mod libration;
pub mod octupole;
pub mod precession;
pub mod published;

#[cfg(test)]
mod integration_tests {
    use crate::libration::{classify_trajectory, ApsidalMotion, Order};
    use crate::octupole::{numerical_octupole_amplitude, octupole_strength_ratio};
    use crate::published;
    use p9_core::constants::{EARTH_MASS_SOLAR, GM_SUN};

    fn p9_gm() -> f64 {
        published::P9_MASS_EARTH * EARTH_MASS_SOLAR * GM_SUN
    }

    #[test]
    fn prototype_octupole_strength_is_order_unity() {
        // The crate's central claim: for the clustered-ETNO prototype the
        // octupole is NOT a small correction. epsilon_oct ~ O(1).
        let eps = octupole_strength_ratio(250.0, published::P9_A_AU, published::P9_E);
        assert!(
            (0.5..2.0).contains(&eps),
            "prototype epsilon_oct = {eps:.3} expected order unity"
        );
    }

    #[test]
    fn octupole_drives_libration_for_prototype_orbit() {
        // End-to-end: a prototype ETNO orbit circulates under quadrupole-only
        // and librates with the octupole on -- the OCTUPOLE_PRODUCES_LIBRATION
        // claim, computed rather than asserted as a constant.
        let gm = p9_gm();
        let (a, a_p, e_p) = (250.0, published::P9_A_AU, published::P9_E);
        let dt = 2.0e6;
        let (quad, _, _) =
            classify_trajectory(Order::Quadrupole, a, 0.3, 0.3, a_p, e_p, gm, dt, 400_000);
        let (oct, _, _) =
            classify_trajectory(Order::Octupole, a, 0.3, 0.3, a_p, e_p, gm, dt, 400_000);
        assert_eq!(quad, ApsidalMotion::Circulation);
        assert_eq!(oct, ApsidalMotion::Libration);
    }

    #[test]
    fn ring_average_apsidal_amplitude_nonzero_for_prototype() {
        // The exact ring average carries a real cos(dvarpi) term for the
        // prototype geometry -- the physical origin of the libration island.
        let gm = p9_gm();
        let a1 = numerical_octupole_amplitude(
            250.0,
            0.2,
            published::P9_A_AU,
            published::P9_E,
            gm,
            128,
            0.01 * published::P9_A_AU,
            48,
        );
        assert!(a1.abs() > 0.0, "expected nonzero apsidal amplitude");
    }
}
