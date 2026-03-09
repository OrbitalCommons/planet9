//! Color and albedo models from Table 1 of Belyakov et al. (2022).
//!
//! Five different surface composition assumptions for Planet Nine, each
//! with a geometric albedo, g-r color, and resulting recovery rate in DES.

use serde::{Deserialize, Serialize};

/// A surface color/albedo model for Planet Nine.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ColorModel {
    /// Descriptive name of the model
    pub name: String,
    /// Geometric albedo (dimensionless)
    pub albedo: f64,
    /// g - r color index (magnitudes)
    pub g_minus_r: f64,
    /// Recovery rate in DES for this color model (fraction)
    pub recovery_rate: f64,
}

/// Fiducial model: moderate albedo, neutral color.
/// Represents a generic icy/rocky surface at ~500 AU.
pub fn fiducial() -> ColorModel {
    ColorModel {
        name: "Fiducial".to_string(),
        albedo: 0.3,
        g_minus_r: 0.5,
        recovery_rate: 0.87,
    }
}

/// Neptune-like model: high albedo, blue color.
/// Assumes a hydrogen/helium-dominated atmosphere.
pub fn neptune_like() -> ColorModel {
    ColorModel {
        name: "Neptune-like".to_string(),
        albedo: 0.41,
        g_minus_r: 0.35,
        recovery_rate: 0.92,
    }
}

/// Methane at 40K model: moderate albedo, very red.
/// Assumes methane ice surface at distant equilibrium temperature.
pub fn methane_40k() -> ColorModel {
    ColorModel {
        name: "CH4 at 40K".to_string(),
        albedo: 0.25,
        g_minus_r: 0.8,
        recovery_rate: 0.83,
    }
}

/// Super-Ganymede model: low albedo, neutral color.
/// Based on scaling up Ganymede's surface properties.
pub fn super_ganymede() -> ColorModel {
    ColorModel {
        name: "Super-Ganymede".to_string(),
        albedo: 0.43,
        g_minus_r: 0.45,
        recovery_rate: 0.90,
    }
}

/// Super-KBO model: very low albedo, ultra-red.
/// Assumes surface similar to distant Kuiper Belt objects.
pub fn super_kbo() -> ColorModel {
    ColorModel {
        name: "Super-KBO".to_string(),
        albedo: 0.12,
        g_minus_r: 1.0,
        recovery_rate: 0.78,
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
    fn fiducial_matches_paper() {
        let m = fiducial();
        assert!((m.albedo - 0.3).abs() < 1e-10);
        assert!((m.g_minus_r - 0.5).abs() < 1e-10);
        assert!((m.recovery_rate - 0.87).abs() < 1e-10);
    }

    #[test]
    fn recovery_rates_in_valid_range() {
        for model in all_models() {
            assert!(
                model.recovery_rate > 0.0 && model.recovery_rate <= 1.0,
                "{} recovery rate {} out of range",
                model.name,
                model.recovery_rate
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
    fn neptune_like_highest_recovery() {
        let models = all_models();
        let max_rate = models
            .iter()
            .map(|m| m.recovery_rate)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((neptune_like().recovery_rate - max_rate).abs() < 1e-10);
    }
}
