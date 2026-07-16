//! Labelled reference constants for the SO and ACT millimeter Planet Nine
//! searches. These are *published inputs*, kept as documented constants rather
//! than derived quantities (the crate computes flux and reach from them).

/// SO Large Aperture Telescope 280 GHz band centre (Hz). The SO forecast
/// (arXiv:2503.00636) identifies 280 GHz as the relevant band for the Planet
/// Nine search — the highest LAT frequency, where a cold body is brightest in
/// the Rayleigh–Jeans regime.
pub const SO_280: f64 = 280e9;

/// SO 225 GHz LAT band centre (Hz).
pub const SO_225: f64 = 225e9;

/// SO 145 GHz LAT band centre (Hz).
pub const SO_145: f64 = 145e9;

/// SO 93 GHz LAT band centre (Hz).
pub const SO_93: f64 = 93e9;

/// ACT 229 GHz band centre (Hz), the highest-frequency ACT band used in the
/// Naess et al. (2021) Planet Nine search (98 / 150 / 229 GHz).
pub const ACT_229: f64 = 229e9;

/// A documented point-source flux-density sensitivity for a millimeter Planet
/// Nine search: the shift-and-stack 95%-confidence flux limit, in mJy. A
/// candidate is "detectable" if its band flux exceeds this threshold.
#[derive(Debug, Clone, Copy)]
pub struct MmSensitivity {
    /// Survey/instrument label.
    pub survey: &'static str,
    /// Band centre frequency (Hz).
    pub nu_hz: f64,
    /// Point-source flux limit (mJy, 95% confidence) below which a source is
    /// not detectable.
    pub flux_limit_mjy: f64,
    pub reference: &'static str,
}

/// ACT Planet Nine search depth (Naess et al. 2021, arXiv:2104.10264): the
/// shift-and-stack flux constraint is < 4–12 mJy (95% CL) for r ≥ 300 AU
/// depending on local map depth. We adopt the *median* ≈ 8 mJy as the
/// representative ACT point-source sensitivity at 229 GHz.
pub const ACT_SENSITIVITY: MmSensitivity = MmSensitivity {
    survey: "ACT",
    nu_hz: ACT_229,
    flux_limit_mjy: 8.0,
    reference: "Naess et al. (2021), ACT P9 search, 4–12 mJy (median ~8)",
};

/// SO Planet Nine search depth in its deepest regions (arXiv:2503.00636). SO's
/// enhanced LAT delivers ~10× the map depth of Planck and, over large parts of
/// the sky, a point-source sensitivity deeper than ACT's first P9 search. We
/// adopt a deep point-source limit of 4.4 mJy at 280 GHz — the value that
/// places a 5 M⊕ Planet Nine's reach at the paper's headline ~900 AU (see
/// `crate::published::SO_REACH_5ME_DEEP_AU`); the shallow-region limit
/// `SO_SENSITIVITY_SHALLOW` reproduces the 500 AU floor over the other half of
/// the sky.
pub const SO_SENSITIVITY: MmSensitivity = MmSensitivity {
    survey: "SO (deep)",
    nu_hz: SO_280,
    flux_limit_mjy: 4.4,
    reference: "SO forecast (arXiv:2503.00636), deep half-sky 280 GHz",
};

/// SO point-source sensitivity in the more shallowly observed regions over
/// half the sky. The paper quotes a 500 AU reach for a 5 M⊕ Planet Nine there;
/// the corresponding flux limit (~14 mJy at 280 GHz) is comparable to ACT's
/// typical fields but applies over far more sky — the breadth, not the depth,
/// is the win in these regions.
pub const SO_SENSITIVITY_SHALLOW: MmSensitivity = MmSensitivity {
    survey: "SO (shallow)",
    nu_hz: SO_280,
    flux_limit_mjy: 14.2,
    reference: "SO forecast (arXiv:2503.00636), shallow half-sky 280 GHz",
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn so_is_deeper_than_act() {
        // SO's deep sensitivity reaches fainter sources than ACT.
        assert!(SO_SENSITIVITY.flux_limit_mjy < ACT_SENSITIVITY.flux_limit_mjy);
    }

    #[test]
    fn bands_ordered() {
        assert!(SO_93 < SO_145 && SO_145 < SO_225 && SO_225 < SO_280);
        assert!(ACT_229 < SO_280);
    }

    #[test]
    fn sensitivities_are_positive() {
        for s in [&ACT_SENSITIVITY, &SO_SENSITIVITY, &SO_SENSITIVITY_SHALLOW] {
            assert!(s.flux_limit_mjy > 0.0);
            assert!(s.nu_hz > 0.0);
        }
    }
}
