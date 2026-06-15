//! Reproduction of Anderson & Kaib (2021),
//! "The effect of Planet Nine on the inclinations of the detached Kuiper belt"
//! (arXiv:2109.13307).
//!
//! HEADLINE: an inclined Planet Nine secularly excites and **broadens** the
//! inclination distribution of the detached Kuiper belt, producing more
//! moderately/highly inclined detached objects than a Planet-Nine-free Solar
//! System — a testable signature.
//!
//! This crate computes that signature from real secular theory on a seeded
//! synthetic detached population, reusing p9-core:
//!
//!   * [`forcing`] — linear Laplace-Lagrange forced inclination of a detached
//!     particle, built from p9-core's quadrupole secular coupling
//!     (`analysis::secular`) for Planet Nine and the giant-planet nodal
//!     precession rate (same closed form as p9-2025-clustering's
//!     `orbital_poles`).
//!   * [`population`] — a `StdRng`-seeded detached disk, its observed
//!     inclination distribution with and without Planet Nine, and the
//!     dispersion / high-i-tail statistics.
//!
//! Published claims are kept as labelled reference constants in [`reference`]
//! and the residuals against our reduced linear-secular model are documented
//! honestly in the integration tests below.

pub mod forcing;
pub mod population;

/// Labelled reference values from Anderson & Kaib (2021). These are *targets*
/// for the qualitative reproduction, not inputs to the calculation — the crate
/// computes its own inclination statistics and the tests compare against these.
pub mod reference {
    use p9_core::constants::DEG2RAD;

    /// Planet Nine inclination used by Anderson & Kaib (2021), ≈ 20° to the
    /// invariable plane (their adopted Batygin-Brown-class orbit).
    pub const P9_INCLINATION_DEG: f64 = 20.0;

    /// Planet Nine mass scale, ≈ 10 M⊕ (nominal Batygin & Brown class).
    pub const P9_MASS_EARTH: f64 = 10.0;

    /// The 8-planet (Kozai-Lidov) mechanism does *not* deposit low-inclination
    /// (i ≲ 20°) detached objects; an inclined Planet Nine does. So the key
    /// inclination scale separating the two scenarios is ~20°.
    pub const KOZAI_LIDOV_FLOOR_DEG: f64 = 20.0;

    /// Detached-disk semi-major-axis lower bound studied (AU).
    pub const DETACHED_A_MIN_AU: f64 = 250.0;

    /// Perihelion threshold for "detached" (decoupled from Neptune, AU).
    pub const DETACHED_Q_MIN_AU: f64 = 40.0;

    /// Helper: Planet Nine inclination in radians.
    pub fn p9_inclination_rad() -> f64 {
        P9_INCLINATION_DEG * DEG2RAD
    }
}

/// The nominal inclined Planet Nine adopted in this reproduction: a
/// Batygin-Brown-class orbit at the [`reference`] mass and inclination.
pub fn nominal_inclined_p9() -> p9_core::types::P9Params {
    p9_core::types::P9Params {
        i: reference::p9_inclination_rad(),
        mass_earth: reference::P9_MASS_EARTH,
        ..p9_core::types::P9Params::nominal_2016()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::{DEG2RAD, RAD2DEG};
    use population::*;

    #[test]
    fn test_end_to_end_broadening_signature() {
        // Full pipeline at the paper's nominal inclined Planet Nine: the
        // detached inclination distribution is broader and hotter with P9.
        let cfg = PopulationConfig::default();
        let pop = build_population(&cfg);
        let p9 = nominal_inclined_p9();

        let with = observed_inclinations(&pop, &p9);
        let without = observed_inclinations_no_p9(&pop);

        let sd_with = dispersion(&with) * RAD2DEG;
        let sd_without = dispersion(&without) * RAD2DEG;
        let mean_with = mean(&with) * RAD2DEG;
        let mean_without = mean(&without) * RAD2DEG;

        assert!(
            sd_with > sd_without,
            "broadening: σ_with = {sd_with:.2}° vs σ_without = {sd_without:.2}°"
        );
        assert!(
            mean_with > mean_without,
            "heating: ⟨i⟩_with = {mean_with:.2}° vs ⟨i⟩_without = {mean_without:.2}°"
        );

        // The without-P9 free disk sits below the Kozai-Lidov floor in the mean
        // (cold detached disk, σ ≈ 12°); P9 lifts a tail across it.
        let floor = reference::KOZAI_LIDOV_FLOOR_DEG * DEG2RAD;
        let f_with = fraction_above(&with, floor);
        let f_without = fraction_above(&without, floor);
        assert!(
            f_with > f_without,
            "i>{:.0}° fraction with P9 ({f_with:.3}) exceeds without ({f_without:.3})",
            reference::KOZAI_LIDOV_FLOOR_DEG
        );
    }

    #[test]
    fn test_forcing_strengthens_with_semimajor_axis_in_population() {
        // The effect strengthens at larger a (stronger α² secular coupling):
        // the outer half of the detached disk is forced more than the inner.
        let cfg = PopulationConfig::default();
        let pop = build_population(&cfg);
        let p9 = nominal_inclined_p9();
        let a_mid = 0.5 * (cfg.a_min + cfg.a_max);

        let inner: Vec<f64> = pop
            .iter()
            .filter(|p| p.a < a_mid)
            .map(|p| forcing::forced_inclination(p.a, p.e, p.i_free, &p9))
            .collect();
        let outer: Vec<f64> = pop
            .iter()
            .filter(|p| p.a >= a_mid)
            .map(|p| forcing::forced_inclination(p.a, p.e, p.i_free, &p9))
            .collect();

        let mean_inner = mean(&inner) * RAD2DEG;
        let mean_outer = mean(&outer) * RAD2DEG;
        assert!(
            mean_outer > mean_inner,
            "outer-disk forced i ({mean_outer:.2}°) must exceed inner ({mean_inner:.2}°)"
        );
    }

    #[test]
    fn test_reference_constants_are_self_consistent() {
        assert!((reference::p9_inclination_rad() - 20.0 * DEG2RAD).abs() < 1e-12);
        assert_eq!(nominal_inclined_p9().mass_earth, reference::P9_MASS_EARTH);
        assert!((nominal_inclined_p9().i * RAD2DEG - reference::P9_INCLINATION_DEG).abs() < 1e-9);
    }
}
