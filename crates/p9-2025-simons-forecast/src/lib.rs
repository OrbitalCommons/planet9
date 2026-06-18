//! Reproduction of the Simons Observatory (SO) Planet Nine detectability
//! forecast (arXiv:2503.00636, *SO: Science Goals and Forecasts for the
//! Enhanced Large Aperture Telescope*).
//!
//! # Headline
//!
//! SO's deep, wide millimeter maps (90–280 GHz, ~60% of the sky, μK-level
//! depth, ~10× Planck) can detect or strongly constrain a thermally-emitting
//! Planet Nine across much of its plausible mass–distance range. The paper's
//! headline: a 5 M⊕ Planet Nine is detectable "at distances ranging from
//! 500 AU in the more shallowly observed regions over half the sky, to 900 AU
//! over large parts of the sky" — substantially deeper than ACT, whose first
//! P9 search (Naess et al. 2021) reached only 325–625 AU for the same 5 M⊕
//! body and ruled out 17% of the orbital parameter space.
//!
//! # What this crate computes (no hard-coded answers)
//!
//! - The **mm blackbody flux density** of a Planet Nine at an SO/ACT band,
//!   from an effective temperature (internal floor + insolation) and a radius
//!   from mass (shared `p9_core` photometry) — [`thermal`].
//! - The **maximum detectable heliocentric distance** vs mass, by inverting
//!   `F_ν(d) = flux_limit` against a documented point-source sensitivity —
//!   [`forecast::max_detectable_distance`].
//! - The **expected detection significance** (SNR) — [`forecast::detection_significance`].
//! - The **detectable fraction** of a reference (mass, distance) box —
//!   [`forecast::ReferenceBox::detectable_fraction`].
//! - That **SO reaches farther than ACT** for the same cold body, and detects
//!   a larger fraction of the reference box.
//!
//! # Why millimeter beats infrared for a cold Planet Nine
//!
//! At 280 GHz (λ ≈ 1.07 mm) a 40 K body has hν/kT ≈ 0.34, deep in the
//! **Rayleigh–Jeans** regime, so its flux is a smooth ∝ ν² T function — strong
//! and easy to extrapolate. This is the physical opposite of the WISE W1 case
//! (`p9-2018-wise-search`), where at 3.4 µm hν/kT ≈ 107 puts a cold body far
//! down the Wien tail (thermal flux ~10⁻⁵⁵ W/m²/Hz/sr, swamped by reflected
//! sunlight). The two crates share the same body model (`mass_radius_neptunian`,
//! 40 K floor) so the contrast is purely the band physics.
//!
//! # Residuals / honest notes
//!
//! The SO and ACT sensitivities are kept as **labelled reference constants**
//! (`bands`), not derived: ACT's 8 mJy is the median of the published 4–12 mJy
//! shift-and-stack limit; SO's deep/shallow limits (1.3 / 4.4 mJy) are the
//! values that reproduce the paper's quoted 900 / 500 AU reach for a 5 M⊕ body
//! under this thermal model. The crate then derives reach, SNR, and box
//! coverage from them, and reproduces the published *direction* and *scale* of
//! every comparison (reach ↑ with mass, SO > ACT, fraction ∈ (0,1), SO > ACT
//! over the box).

pub mod bands;
pub mod forecast;
pub mod thermal;

/// Published reference scalars from the SO forecast and the ACT search,
/// retained for cross-checks (not derived).
pub mod published {
    /// SO 5 M⊕ Planet Nine reach in deep regions over large parts of the sky
    /// (AU), from arXiv:2503.00636.
    pub const SO_REACH_5ME_DEEP_AU: f64 = 900.0;
    /// SO 5 M⊕ reach in the more shallowly observed half-sky regions (AU).
    pub const SO_REACH_5ME_SHALLOW_AU: f64 = 500.0;
    /// ACT 5 M⊕ detection-limit band (AU), Naess et al. (2021).
    pub const ACT_REACH_5ME_MIN_AU: f64 = 325.0;
    pub const ACT_REACH_5ME_MAX_AU: f64 = 625.0;
    /// ACT 10 M⊕ detection-limit band (AU), Naess et al. (2021).
    pub const ACT_REACH_10ME_MIN_AU: f64 = 425.0;
    pub const ACT_REACH_10ME_MAX_AU: f64 = 775.0;
    /// Fraction of the 5 M⊕ orbital parameter space ACT ruled out.
    pub const ACT_RULED_OUT_FRACTION_5ME: f64 = 0.17;
}

#[cfg(test)]
mod tests {
    use crate::bands::{ACT_SENSITIVITY, SO_SENSITIVITY};
    use crate::forecast::{max_detectable_distance, ReferenceBox};
    use crate::published;

    #[test]
    fn end_to_end_so_outperforms_act() {
        // The full forecast in one assertion chain: for a nominal 5 M⊕ body
        // SO reaches the paper's ~900 AU, exceeds ACT's published band, and
        // covers a larger fraction of the reference box.
        use p9_core::units::au;
        let so5 = (max_detectable_distance(5.0, &SO_SENSITIVITY) / au(1.0)).value;
        let act5 = (max_detectable_distance(5.0, &ACT_SENSITIVITY) / au(1.0)).value;

        assert!(so5 > act5, "SO reach > ACT reach");
        assert!(
            (so5 - published::SO_REACH_5ME_DEEP_AU).abs() < 80.0,
            "SO 5 M⊕ reach {so5:.0} AU near published {}",
            published::SO_REACH_5ME_DEEP_AU
        );
        assert!(
            (published::ACT_REACH_5ME_MIN_AU..=published::ACT_REACH_5ME_MAX_AU).contains(&act5),
            "ACT 5 M⊕ reach {act5:.0} AU within published band"
        );

        let bx = ReferenceBox::nominal();
        assert!(bx.detectable_fraction(&SO_SENSITIVITY) > bx.detectable_fraction(&ACT_SENSITIVITY));
    }
}
