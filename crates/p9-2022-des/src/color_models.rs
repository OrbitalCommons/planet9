//! Surface color and albedo models from Table 1 of Belyakov, Bernardinelli
//! & Brown (2022), AJ 163, 216 (arXiv:2203.07642).
//!
//! Five surface-composition hypotheses for Planet Nine, each defined by a
//! geometric albedo and g-band-relative colors (g-r, g-i, g-z). The colors
//! and albedo act on detectability *through the band-dependent depths* of
//! the survey model: a red surface (Super-KBO) is fainter in g where DES is
//! deepest; an icy CH4 surface is bright everywhere. `paper_recovery` is the
//! published end-to-end recovery rate for each model; the crate's computed
//! counterpart lives in `recovery_analysis::compute_recovery_for_model`.

use serde::{Deserialize, Serialize};

use p9_core::analysis::photometry::planet_apparent_magnitude;

use crate::survey_model::DesBand;

/// Apparent magnitudes of an object in the DES bands.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct BandMagnitudes {
    pub g: f64,
    pub r: f64,
    pub i: f64,
    pub z: f64,
}

impl BandMagnitudes {
    /// Magnitude in a given DES band (Y is approximated by z; the detection
    /// model does not use Y).
    pub fn magnitude(&self, band: DesBand) -> f64 {
        match band {
            DesBand::G => self.g,
            DesBand::R => self.r,
            DesBand::I => self.i,
            DesBand::Z => self.z,
            DesBand::Y => self.z,
        }
    }
}

/// A surface color/albedo model for Planet Nine (Table 1).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ColorModel {
    /// Descriptive name of the model
    pub name: String,
    /// Geometric albedo (dimensionless)
    pub albedo: f64,
    /// g - r color index (magnitudes)
    pub g_minus_r: f64,
    /// g - i color index (magnitudes)
    pub g_minus_i: f64,
    /// g - z color index (magnitudes)
    pub g_minus_z: f64,
    /// Published recovery rate in DES for this model (Table 1)
    pub paper_recovery: f64,
}

impl ColorModel {
    /// Apparent magnitudes in the DES bands for a planet of `mass_earth`
    /// Earth masses at heliocentric distance `r_au` (opposition).
    ///
    /// The V magnitude comes from `p9_core::analysis::photometry`
    /// (Neptune-anchored mass-radius relation + reflected-sunlight
    /// 5 log10(r*Delta) law, replacing the stellar 5 log10(d/10) law that
    /// made P9 naked-eye bright). The V -> g conversion uses the Jester et
    /// al. (2005) photometric transformation V = g - 0.59(g-r) - 0.01,
    /// evaluated with the model's own g-r color, and the remaining bands
    /// follow from the model's g-relative colors.
    pub fn band_magnitudes(&self, mass_earth: f64, r_au: f64) -> BandMagnitudes {
        let v = planet_apparent_magnitude(mass_earth, self.albedo, r_au);
        let g = v + 0.59 * self.g_minus_r + 0.01;
        BandMagnitudes {
            g,
            r: g - self.g_minus_r,
            i: g - self.g_minus_i,
            z: g - self.g_minus_z,
        }
    }
}

/// Fiducial model of Brown & Batygin (2021): albedo 0.44, near-solar colors.
pub fn fiducial() -> ColorModel {
    ColorModel {
        name: "Fiducial (BB21)".to_string(),
        albedo: 0.44,
        g_minus_r: 0.53,
        g_minus_i: 0.55,
        g_minus_z: 0.55,
        paper_recovery: 0.870,
    }
}

/// Neptune-like model: albedo 0.5, strongly blue (bright in g, very faint
/// in i and z from methane absorption).
pub fn neptune_like() -> ColorModel {
    ColorModel {
        name: "Neptune-like".to_string(),
        albedo: 0.5,
        g_minus_r: -0.3,
        g_minus_i: -1.35,
        g_minus_z: -2.15,
        paper_recovery: 0.803,
    }
}

/// 40 K surface with 10% methane frost: high albedo 0.75, mildly red in
/// g-r, blue in i/z.
pub fn methane_40k() -> ColorModel {
    ColorModel {
        name: "40K, 10% CH4".to_string(),
        albedo: 0.75,
        g_minus_r: 0.6,
        g_minus_i: -0.3,
        g_minus_z: -0.2,
        paper_recovery: 0.881,
    }
}

/// Super-Ganymede model: albedo 0.43, moderately red.
pub fn super_ganymede() -> ColorModel {
    ColorModel {
        name: "Super-Ganymede".to_string(),
        albedo: 0.43,
        g_minus_r: 0.72,
        g_minus_i: 0.88,
        g_minus_z: 0.87,
        paper_recovery: 0.858,
    }
}

/// Super-KBO model: dark (albedo 0.1) and ultra-red.
pub fn super_kbo() -> ColorModel {
    ColorModel {
        name: "Super-KBO".to_string(),
        albedo: 0.1,
        g_minus_r: 1.0,
        g_minus_i: 1.25,
        g_minus_z: 1.5,
        paper_recovery: 0.768,
    }
}

/// Returns all five color models from Table 1.
pub fn all_models() -> Vec<ColorModel> {
    vec![
        fiducial(),
        neptune_like(),
        methane_40k(),
        super_ganymede(),
        super_kbo(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_models_has_five_entries() {
        assert_eq!(all_models().len(), 5);
    }

    #[test]
    fn fiducial_matches_paper_table1() {
        let m = fiducial();
        assert!((m.albedo - 0.44).abs() < 1e-10);
        assert!((m.g_minus_r - 0.53).abs() < 1e-10);
        assert!((m.paper_recovery - 0.870).abs() < 1e-10);
    }

    #[test]
    fn paper_recovery_rates_in_valid_range() {
        for model in all_models() {
            assert!(
                model.paper_recovery > 0.0 && model.paper_recovery <= 1.0,
                "{} recovery rate {} out of range",
                model.name,
                model.paper_recovery
            );
        }
    }

    #[test]
    fn albedos_in_valid_range() {
        for model in all_models() {
            assert!(
                model.albedo > 0.0 && model.albedo < 1.0,
                "{} albedo {} out of range",
                model.name,
                model.albedo
            );
        }
    }

    #[test]
    fn methane_highest_paper_recovery_super_kbo_lowest() {
        let models = all_models();
        let max_rate = models
            .iter()
            .map(|m| m.paper_recovery)
            .fold(f64::NEG_INFINITY, f64::max);
        let min_rate = models
            .iter()
            .map(|m| m.paper_recovery)
            .fold(f64::INFINITY, f64::min);
        assert!((methane_40k().paper_recovery - max_rate).abs() < 1e-10);
        assert!((super_kbo().paper_recovery - min_rate).abs() < 1e-10);
    }

    #[test]
    fn band_magnitudes_respect_colors_and_distance_law() {
        // Super-KBO is fainter in g than r by its g-r = 1.0.
        let mags = super_kbo().band_magnitudes(6.2, 500.0);
        assert!((mags.g - mags.r - 1.0).abs() < 1e-12);
        // Reflected light: doubling distance dims by ~3 mag.
        let far = super_kbo().band_magnitudes(6.2, 1000.0);
        assert!((far.g - mags.g - 3.01).abs() < 0.05);
        // Dark red surface is fainter than the bright fiducial in g.
        let fid = fiducial().band_magnitudes(6.2, 500.0);
        assert!(mags.g > fid.g);
    }
}
